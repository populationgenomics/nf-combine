#!/usr/bin/env python3

"""Writes a bash script of `gcloud storage objects compose` commands that join VCF fragments in GCS into one VCF,
without downloading any of it.

Ported from cpg-flow-seqr-loader's scripts/write_gcloud_compose_script.py, with only the standard library, so that it
runs in the gcloud container. The manifest lists the header first, then the shards in genomic order: either names in
--vcf_dir (a separate-header export), or whole gs:// paths. Because the header and shards are BGZF blocks, joining the
objects in that order gives a valid bgzipped VCF.

A compose takes at most 32 objects, so the fragments are composed in rounds of 32 into temp objects, until 32 or
fewer remain. The last round composes straight to the output, which must be in the same bucket as the fragments. If
--final_output is given, the output is then copied there, which may be in any bucket.
"""

import argparse
import logging
import os

logger = logging.getLogger(__name__)

COMPOSE = 'gcloud storage objects compose'
COPY = 'gcloud storage cp'

# the most source objects one compose accepts
COMPOSE_LIMIT = 32

# an object's line in `gcloud storage ls -l`: size, timestamp, gs:// path
LISTING_FIELDS = 3

# a BGZF file containing only the end-of-file marker block, i.e. a fragment from an empty partition. These must be
# skipped: htsjdk (GATK) reads an EOF block mid-file as the end of the file, silently truncating the VCF
BGZF_EOF_BYTES = 28


def read_fragment_sizes(listing_file: str) -> dict[str, int]:
    """Parses `gcloud storage ls -l` output into fragment sizes.

    Args:
        listing_file: The saved listing. Its summary line and any non-object lines are ignored.

    Returns:
        {object basename: size in bytes}
    """
    sizes = {}
    with open(listing_file) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) == LISTING_FIELDS and fields[0].isdigit() and fields[2].startswith('gs://'):
                sizes[os.path.basename(fields[2])] = int(fields[0])
    return sizes


def compose_commands(fragments: list[str], output: str, temp_dir: str) -> list[str]:
    """The compose commands that join fragments, in order, into output.

    Args:
        fragments: gs:// paths, in the order they're to be joined.
        output: Where the joined object is written.
        temp_dir: Where intermediate rounds write, e.g. gs://bucket/compose. Overwritten on a rerun.

    Returns:
        One command per compose, in the order they must run.
    """
    commands = []
    merge_round = 1
    while len(fragments) > COMPOSE_LIMIT:
        merged = []
        for index in range(0, len(fragments), COMPOSE_LIMIT):
            target = f'{temp_dir}/round_{merge_round}/chunk_{index // COMPOSE_LIMIT}.vcf.bgz'
            commands.append(f'{COMPOSE} {" ".join(fragments[index : index + COMPOSE_LIMIT])} {target}')
            merged.append(target)
        logger.info(f'Round {merge_round}: {len(fragments)} fragments -> {len(merged)}')
        fragments = merged
        merge_round += 1
    commands.append(f'{COMPOSE} {" ".join(fragments)} {output}')
    return commands


def main(argv: list[str] | None = None) -> None:
    """Writes the compose script for one set of fragments.

    Args:
        argv: The command line, without the program name. None: sys.argv.

    Raises:
        ValueError: If a manifest entry isn't in the listing, or no fragment holds any variants.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--manifest', required=True, help='the fragments to join, header first')
    parser.add_argument('--vcf_dir', help='the directory the manifest names are in. Unset: the names are gs:// paths')
    parser.add_argument(
        '--sizes',
        help='gcloud storage ls -l output covering the fragments: those of 28 bytes or less (a BGZF EOF block at most) '
        'are skipped. Unset: none are',
    )
    parser.add_argument('--output', required=True, help='the joined VCF, in the same bucket as the fragments')
    parser.add_argument('--temp_dir', required=True, help='where intermediate compose rounds write')
    parser.add_argument('--final_output', help='where the joined VCF is copied to, if anywhere')
    parser.add_argument('--script', required=True, help='the bash script to write')
    args = parser.parse_args(argv)

    with open(args.manifest) as handle:
        names = [line.strip() for line in handle if line.strip()]
    fragments = [f'{args.vcf_dir.rstrip("/")}/{name}' for name in names] if args.vcf_dir else names

    if args.sizes:
        # keyed by basename, as gcloud lists them
        sizes = read_fragment_sizes(args.sizes)
        if missing := [path for path in fragments if os.path.basename(path) not in sizes]:
            raise ValueError(f'{len(missing)} fragments in the manifest are not in the listing: {missing[:10]}')

        if empty := [path for path in fragments if sizes[os.path.basename(path)] <= BGZF_EOF_BYTES]:
            logger.info(f'Skipping {len(empty)} empty fragments')

        fragments = [path for path in fragments if sizes[os.path.basename(path)] > BGZF_EOF_BYTES]
        # the header is always one fragment: anything less than two means there are no variants
        if len(fragments) < 2:  # noqa: PLR2004
            raise ValueError(f'No fragment in {args.manifest} holds any variants')

    commands = compose_commands(fragments, args.output, args.temp_dir.rstrip('/'))
    if args.final_output:
        commands.append(f'{COPY} {args.output} {args.final_output}')

    with open(args.script, 'w') as handle:
        handle.write('#!/bin/bash\nset -euo pipefail\n\n')
        for command in commands:
            handle.write(f'{command}\n')


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()
