# One output dataset per run, but its cohorts may span datasets

A run is launched for one `dataset` and one `sequencing_type`, and it writes and registers its single VDS
in that dataset. Cohort lookup is **not** scoped to the same dataset, so SGs may belong to other datasets. 
Joint VDSs routinely cross projects, so the operating service account requires read and Metamist access across 
all datasets referenced in input cohorts, but only requires write access to the one dataset (for now). 

## Consequences

Each sample sheet row records the SG's own dataset. The combiner stage doesn't use it. Later stages (not yet built)
are expected to fan outputs back out to each SG's dataset bucket, which the Dataproc cluster can write to directly.
