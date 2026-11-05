# Hail runs on an ephemeral Dataproc cluster, driven by a non-spot submitter task

A small Nextflow task creates a Dataproc cluster with `cpg_utils.dataproc_runner`, runs the combiner on it,
and tears it down. This is the approach already proven in `cohort_demo`, and it scales to cohort-sized
combines. The submitter lives as long as the combine does, which can be many hours. If it were preempted,
its SIGTERM handler would tear the cluster down. So it runs on a regular VM. Nextflow sets spot only for a
whole pipeline (`google.batch.spot`), not per process, so this pipeline's `standard` profile turns spot off
altogether, unlike `nf-align`'s. The only other task is the small `BUILD_SAMPLESHEET`. It gets one retry, which resumes from the combiner's saved plan.

## Considered Options

- **Hail Query-on-Batch**: no cluster to manage, but not what the existing combiner code has been run on.
- **Spark in local mode on one large Batch VM**: simple, but only viable for small cohorts.

## Consequences

The Hail version comes from the cluster image `HailDataprocCluster` selects, not from `cpg_container`.
Registration happens inside the submitter once the job succeeds. The VDS is written directly to its final
path, so there's no asynchronous Nextflow publish to wait for.
