#!/usr/bin/env python3

"""Dataproc lifecycle for one Hail job: create a Hail cluster, run a Hail script on it, and tear it down.

Run by the COMBINE and DENSIFY processes (see DESIGN.md). The submitter's own options come first, then
`--`, then the Hail script's arguments, which are passed on unchanged:

    submit_hail.py --cluster-name ... --script combiner.py --pyfiles nf_combine_hail \\
        --complete gs://.../x.vds/variant_data/_SUCCESS ... -- --out gs://.../x.vds ...

- --pyfiles: package directories or .py files, uploaded with HailDataprocCluster.upload_pyfiles. The script can
  import them on the cluster.
- --files: local files, localised by basename into the job's working directory on the cluster.
- --complete: the paths that exist once the job's outputs are complete. If they all exist already, as on a
  retry after the job succeeded but a later step failed, no cluster is created. If any is missing once the job
  finishes, the run fails.
- --fetch REMOTE LOCAL: copies a file the job wrote into the working directory, for the process's outputs.

Everything is uploaded before the cluster starts, so a bad path fails without the cost of starting a cluster.
Writes --output: a JSON record of the cluster and job (null if no cluster was needed), and any --record values.
"""

import argparse
import json
import logging
import sys

from cpg_utils import to_path
from cpg_utils.dataproc_runner import HailDataprocCluster, install_sigterm_cleanup, parse_label_kvs

logger = logging.getLogger(__name__)


def split_argv(argv: list[str]) -> tuple[list[str], list[str]]:
    """Splits the command line at the first `--` into the submitter's arguments and the Hail script's.

    Args:
        argv: The command line, without the program name.

    Returns:
        The submitter's arguments, and the Hail script's (empty if there's no `--`).

    >>> split_argv(['--script', 's.py', '--', '--out', 'x', '--', 'y'])
    (['--script', 's.py'], ['--out', 'x', '--', 'y'])
    >>> split_argv(['--script', 's.py'])
    (['--script', 's.py'], [])
    """
    if '--' not in argv:
        return argv, []
    index = argv.index('--')
    return argv[:index], argv[index + 1 :]


def parse_records(items: list[str]) -> dict[str, str | None]:
    """Parses key=value pairs for the output record. An empty value is recorded as null.

    Args:
        items: The key=value strings.

    Returns:
        The values, by key.

    Raises:
        ValueError: If an item has no '='.

    >>> parse_records(['vds=gs://b/x.vds', 'base_vds='])
    {'vds': 'gs://b/x.vds', 'base_vds': None}
    """
    if bad := [item for item in items if '=' not in item]:
        raise ValueError(f'--record values must be key=value, got {bad}')
    return {key: value or None for key, value in (item.split('=', 1) for item in items)}


def is_complete(paths: list[str]) -> bool:
    """Checks whether every completion marker exists.

    Args:
        paths: The completion markers.

    Returns:
        True if they all exist.
    """
    return all(to_path(path).exists() for path in paths)


def parse_args(argv: list[str]) -> argparse.Namespace:
    """Parses the submitter's own arguments, as built by the COMBINE and DENSIFY processes.

    Args:
        argv: The arguments before `--`.

    Returns:
        The parsed arguments.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--cluster-name', required=True, help='Prefix for the cluster name; must start with a letter')
    parser.add_argument('--region', required=True, help='GCP region')
    parser.add_argument('--project', required=True, help='GCP project the cluster runs in')
    parser.add_argument('--staging-bucket', required=True, help='GCS bucket the script and its files are staged to')
    parser.add_argument('--temp-bucket', required=True, help='GCS bucket for Dataproc temp files')
    parser.add_argument('--autoscaling-policy', default=None, help='Autoscaling policy ID or full resource URI')
    parser.add_argument('--labels', nargs='*', default=[], help='Cluster labels as key=value pairs')

    # cluster sizing and lifecycle, forwarded to HailDataprocCluster
    parser.add_argument('--master-type', required=True)
    parser.add_argument('--worker-type', required=True)
    parser.add_argument('--num-workers', type=int, required=True)
    parser.add_argument('--num-secondary-workers', type=int, required=True)
    parser.add_argument('--preemptible-secondary-workers', action='store_true')
    parser.add_argument('--max-age', type=int, required=True, help='Cluster auto-delete TTL (seconds)')
    parser.add_argument('--max-idle', type=int, required=True, help='Cluster idle-delete TTL (seconds)')
    parser.add_argument('--boot-disk-size', type=int, required=True, help='Per-node boot disk (GB)')
    parser.add_argument('--num-local-ssds', type=int, required=True, help='Local SSDs per node, for shuffle')
    parser.add_argument('--service-account', required=True, help='Identity the cluster VMs run as')

    # the job
    parser.add_argument('--script', required=True, help='Local path to the Hail script')
    parser.add_argument('--pyfiles', nargs='*', default=[], help='Package directories or .py files to import')
    parser.add_argument('--files', nargs='*', default=[], help="Local files for the job's working directory")
    parser.add_argument('--complete', nargs='+', required=True, help='Paths that exist once the outputs are complete')
    parser.add_argument(
        '--fetch',
        nargs=2,
        action='append',
        default=[],
        metavar=('REMOTE', 'LOCAL'),
        help='A file the job writes, to copy into the working directory afterwards. Repeatable',
    )
    parser.add_argument('--output', required=True, help='Local path for the JSON record of the cluster and job')
    parser.add_argument('--record', nargs='*', default=[], help='key=value pairs added to the --output record')
    return parser.parse_args(argv)


def run_on_cluster(args: argparse.Namespace, hail_args: list[str]) -> dict:
    """Runs the Hail script on a new cluster, and tears the cluster down.

    Args:
        args: The submitter's arguments.
        hail_args: The Hail script's arguments.

    Returns:
        The cluster and job, for the output record.
    """
    cluster = HailDataprocCluster(
        project=args.project,
        region=args.region,
        autoscaling_policy=args.autoscaling_policy,
        cluster_name_prefix=args.cluster_name,
        staging_bucket=args.staging_bucket,
        temp_bucket=args.temp_bucket,
        labels=parse_label_kvs(args.labels) if args.labels else None,
        master_type=args.master_type,
        worker_type=args.worker_type,
        num_workers=args.num_workers,
        num_secondary_workers=args.num_secondary_workers,
        max_age_seconds=args.max_age,
        max_idle_seconds=args.max_idle,
        boot_disk_size_gb=args.boot_disk_size,
        preemptible_workers=args.preemptible_secondary_workers,
        num_local_ssds=args.num_local_ssds,
        service_account=args.service_account,
        internal_ip_only=False,
    )

    # none of these start the cluster
    script_uri = cluster.upload(args.script)
    file_uris = [cluster.upload(path) for path in args.files]
    pyfile_uris = cluster.upload_pyfiles(args.pyfiles) if args.pyfiles else []

    # tear the cluster down promptly if Nextflow cancels the task (SIGTERM)
    install_sigterm_cleanup(cluster)
    with cluster:
        job = cluster.run_job(script_uri, args=hail_args, file_uris=file_uris, python_file_uris=pyfile_uris)
        return {
            'name': cluster.name,
            'job_id': job.job_uuid,
            'region': args.region,
            'project': args.project,
            'hail_version': cluster.hail_version,
            'autoscaling_policy_uri': cluster.autoscaling_policy_uri,
            'labels': cluster.labels,
            'pyfiles': pyfile_uris,
        }


def main(argv: list[str] | None = None) -> None:
    """Runs one Hail job on a Dataproc cluster, unless its outputs are already complete, and writes its records.

    Args:
        argv: The command line, without the program name. Defaults to sys.argv[1:].

    Raises:
        RuntimeError: If the outputs are still incomplete once the job finishes.
    """
    own_argv, hail_args = split_argv(sys.argv[1:] if argv is None else argv)
    args = parse_args(own_argv)
    record = parse_records(args.record)

    if is_complete(args.complete):
        logger.info('The outputs are already complete, so no cluster is needed')
        record['cluster'] = None
    else:
        record['cluster'] = run_on_cluster(args, hail_args)

    if not is_complete(args.complete):
        missing = [path for path in args.complete if not to_path(path).exists()]
        raise RuntimeError(f'The job finished, but these outputs are missing: {missing}')

    for remote, local in args.fetch:
        with open(local, 'wb') as handle:
            handle.write(to_path(remote).read_bytes())
    with open(args.output, 'w') as handle:
        json.dump(record, handle, indent=2)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')
    main()
