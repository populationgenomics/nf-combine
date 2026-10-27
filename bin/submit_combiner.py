#!/usr/bin/env python3

"""Dataproc lifecycle for one combine: create a Hail cluster, run combiner.py on it, and tear it down.

Run by the COMBINE process (see DESIGN.md, ADR 0003). If the VDS is already complete, as on a retry whose
first attempt combined but failed to register, no cluster is created. Writes, in the working directory:
- output.json: the cluster and job, for the run records
- combine_report.json: the base VDS used, and the SGs added, removed and refreshed (written by combiner.py)
"""

import argparse
import json
import logging

from cpg_utils import to_path
from cpg_utils.dataproc_runner import HailDataprocCluster, install_sigterm_cleanup, parse_label_kvs

logger = logging.getLogger(__name__)


def vds_is_complete(vds: str) -> bool:
    """Checks whether a VDS is fully written.

    Hail marks each of a VDS's two matrix tables with _SUCCESS once it's fully written.

    Args:
        vds: The VDS path.

    Returns:
        True if both matrix tables are marked complete.
    """
    return all(to_path(f'{vds}/{part}/_SUCCESS').exists() for part in ['reference_data', 'variant_data'])


def parse_args() -> argparse.Namespace:
    """Parses the command line, as built by the COMBINE process.

    Returns:
        The parsed arguments.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cluster-name', required=True, help='Prefix for the cluster name; must start with a letter')
    parser.add_argument('--region', required=True, help='GCP region')
    parser.add_argument('--project', required=True, help='GCP project the cluster runs in')
    parser.add_argument('--script', required=True, help='Local path to the Hail script (combiner.py)')
    parser.add_argument('--samplesheet', required=True, help='Sample sheet TSV: sg_id, gvcf, dataset')
    parser.add_argument('--staging-bucket', required=True, help='GCS bucket the script and sheet are staged to')
    parser.add_argument('--temp-bucket', required=True, help='GCS bucket for Dataproc and combiner temp files')
    parser.add_argument('--vds', required=True, help='VDS path to write')
    parser.add_argument('--intervals', required=True, help='Regions BED; a gs:// URI read directly by Hail')
    parser.add_argument('--base-vds', default=None, help='VDS to build on')
    parser.add_argument('--refresh', nargs='*', default=[], help='SG IDs to refresh from their gVCFs')
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

    # combiner tuning, forwarded to combiner.py
    parser.add_argument('--branch-factor', type=int, required=True)
    parser.add_argument('--target-records', type=int, required=True)
    parser.add_argument('--gvcf-batch-size', type=int, required=True)
    parser.add_argument('--force', action='store_true', help='Discard any saved combiner plan')
    return parser.parse_args()


def main() -> None:
    """Runs one combine on a Dataproc cluster, unless the VDS is already complete, and writes its records.

    Raises:
        ValueError: If --intervals isn't a gs:// path.
        RuntimeError: If the VDS is still incomplete once the combine finishes.
    """
    args = parse_args()
    if not args.intervals.startswith('gs://'):
        raise ValueError(f'--intervals must be a gs:// path, got {args.intervals}')

    # keyed by the VDS name, so a retry finds the plan and report of its earlier attempt
    name = to_path(args.vds).stem
    temp = f'{args.temp_bucket.rstrip("/")}/nf-combine'
    report = f'{temp}/combine_reports/{name}.json'

    output: dict = {'vds': args.vds, 'base_vds': args.base_vds}
    if vds_is_complete(args.vds):
        logger.info(f'{args.vds} is already complete, so no cluster is needed')
        output['cluster'] = None
    else:
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
        # tear the cluster down promptly if Nextflow cancels the task (SIGTERM)
        install_sigterm_cleanup(cluster)

        with cluster:
            script_uri = cluster.upload(args.script)
            samplesheet = cluster.upload(args.samplesheet)
            hail_args = [
                '--samplesheet', to_path(samplesheet).name,
                '--out', args.vds,
                '--temp_dir', temp,
                '--report', report,
                '--intervals', args.intervals,
                '--branch', str(args.branch_factor),
                '--target_records', str(args.target_records),
                '--gvcf_batch', str(args.gvcf_batch_size),
            ]  # fmt: skip
            if args.base_vds:
                hail_args += ['--base_vds', args.base_vds]
            if args.refresh:
                hail_args += ['--refresh', *args.refresh]
            if args.force:
                hail_args.append('--force')
            job = cluster.run_job(script_uri, args=hail_args, file_uris=[samplesheet])
            output['cluster'] = {
                'name': cluster.name,
                'job_id': job.job_uuid,
                'region': args.region,
                'project': args.project,
                'hail_version': cluster.hail_version,
                'autoscaling_policy_uri': cluster.autoscaling_policy_uri,
                'labels': cluster.labels,
            }

    if not vds_is_complete(args.vds):
        raise RuntimeError(f'The combine finished, but {args.vds} is not complete')
    with open('output.json', 'w') as handle:
        json.dump(output, handle, indent=2)
    with open('combine_report.json', 'w') as handle:
        handle.write(to_path(report).read_text())


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')
    main()
