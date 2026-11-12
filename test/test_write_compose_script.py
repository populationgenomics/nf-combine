"""bin/write_compose_script.py: empty shards are skipped, and composes stay within GCS's 32-object limit."""

from pathlib import Path

import pytest
from write_compose_script import BGZF_EOF_BYTES, COMPOSE_LIMIT, compose_commands, main

VCF_DIR = 'gs://bucket/densify/x_separate.vcf.bgz'


def write_export(tmp_path: Path, sizes: dict[str, int]) -> tuple[str, str]:
    """A manifest of the shards in order, and a `gcloud storage ls -l` listing of them."""
    manifest = tmp_path / 'manifest.txt'
    manifest.write_text(''.join(f'{name}\n' for name in sizes))
    listing = tmp_path / 'sizes.txt'
    lines = [f'{size}  2026-11-06T00:00:00Z  {VCF_DIR}/{name}' for name, size in sizes.items()]
    listing.write_text('\n'.join([*lines, f'TOTAL: {len(sizes)} objects, 0 bytes (0B)']) + '\n')
    return str(manifest), str(listing)


def compose_args(manifest: str, listing: str, script: Path) -> list[str]:
    return [
        '--manifest', manifest,
        '--vcf_dir', VCF_DIR,
        '--sizes', listing,
        '--output', 'gs://bucket/out.vcf.bgz',
        '--temp_dir', 'gs://bucket/tmp',
        '--script', str(script),
    ]  # fmt: skip


def run(tmp_path: Path, sizes: dict[str, int]) -> list[str]:
    manifest, listing = write_export(tmp_path, sizes)
    script = tmp_path / 'compose.sh'
    main(compose_args(manifest, listing, script))
    return [line for line in script.read_text().splitlines() if line.startswith('gcloud')]


def test_few_fragments_compose_straight_to_the_output(tmp_path):
    commands = run(tmp_path, {'header.bgz': 500, 'part-0.bgz': 100, 'part-1.bgz': 100})
    assert commands == [
        f'gcloud storage objects compose {VCF_DIR}/header.bgz {VCF_DIR}/part-0.bgz {VCF_DIR}/part-1.bgz '
        'gs://bucket/out.vcf.bgz',
    ]


def test_empty_fragments_are_skipped(tmp_path):
    commands = run(tmp_path, {'header.bgz': 500, 'part-0.bgz': BGZF_EOF_BYTES, 'part-1.bgz': 100})
    assert 'part-0.bgz' not in commands[0]
    assert 'part-1.bgz' in commands[0]


def test_only_empty_fragments_fail(tmp_path):
    with pytest.raises(ValueError, match='holds any variants'):
        run(tmp_path, {'header.bgz': 500, 'part-0.bgz': BGZF_EOF_BYTES})


def test_fragments_missing_from_the_listing_fail(tmp_path):
    manifest, listing = write_export(tmp_path, {'header.bgz': 500})
    Path(manifest).write_text('header.bgz\npart-0.bgz\n')
    with pytest.raises(ValueError, match='not in'):
        main(compose_args(manifest, listing, tmp_path / 'compose.sh'))


def test_many_fragments_compose_in_rounds_in_order():
    fragments = [f'gs://b/{i}' for i in range(COMPOSE_LIMIT**2 + 1)]
    commands = compose_commands(fragments, 'gs://b/out', 'gs://b/tmp')
    for command in commands:
        assert len(command.split()) - 5 <= COMPOSE_LIMIT  # 4 words of command, 1 output
    # 1025 -> 33 -> 2 -> the output
    assert commands[-1].split()[4:] == [
        'gs://b/tmp/round_2/chunk_0.vcf.bgz',
        'gs://b/tmp/round_2/chunk_1.vcf.bgz',
        'gs://b/out',
    ]
    # the first round keeps the manifest order
    assert commands[0].split()[4:6] == ['gs://b/0', 'gs://b/1']


def test_a_manifest_of_paths_composes_every_fragment_then_copies_to_the_final_output(tmp_path):
    """Without --vcf_dir or --sizes: the manifest's gs:// paths are composed as they are, empty or not."""
    paths = ['gs://work/a/header.vcf.bgz', 'gs://work/b/part-0.subset.vcf.bgz', 'gs://work/a/eof.bgz']
    manifest = tmp_path / 'fragments.txt'
    manifest.write_text(''.join(f'{path}\n' for path in paths))
    script = tmp_path / 'compose.sh'
    main([
        '--manifest', str(manifest),
        '--output', 'gs://work/a/compose/out.vcf.bgz',
        '--temp_dir', 'gs://work/a/compose/intermediates',
        '--final_output', 'gs://dataset/joint_vcf/out.vcf.bgz',
        '--script', str(script),
    ])  # fmt: skip
    commands = [line for line in script.read_text().splitlines() if line.startswith('gcloud')]
    assert commands == [
        f'gcloud storage objects compose {" ".join(paths)} gs://work/a/compose/out.vcf.bgz',
        'gcloud storage cp gs://work/a/compose/out.vcf.bgz gs://dataset/joint_vcf/out.vcf.bgz',
    ]


def test_without_a_final_output_nothing_is_copied(tmp_path):
    assert not any(command.startswith('gcloud storage cp') for command in run(tmp_path, {'h': 500, 'p': 100}))
