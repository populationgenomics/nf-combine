# One output dataset per run, but its cohorts may span datasets

A run is launched for one `dataset` and one `sequencing_type`, and it writes and registers its single VDS
in that dataset. Its cohorts are looked up in that dataset's Metamist project, but their SGs may belong to
other datasets. Joint VDSs routinely cross projects, and the operating service account has read, write and
Metamist access across them. This is a deliberate departure from nf-align's ADR 0001, which fails a run on
any SG from another project.

## Consequences

Each sample sheet row records the SG's own dataset. This stage doesn't use it. Later stages (not yet built)
are expected to fan outputs back out to each SG's dataset bucket, which the Dataproc cluster can write to
directly.
