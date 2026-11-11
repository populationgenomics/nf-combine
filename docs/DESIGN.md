# nf-combine: design brief

Combines the single-sample gVCFs of a set of Metamist cohorts into one Hail VDS. A run normally builds on
the newest earlier VDS for the dataset: it removes withdrawn SGs, refreshes changed ones, and adds new ones.
The VDS is then densified into VCF exports, filtered with VQSR, and split into one VCF per dataset.
It's a Nextflow (≥ 26.04, strict syntax) pipeline. Small tasks run on Google Batch, and the combine and densify
each run on their own ephemeral Dataproc cluster. It can be tested offline with stub runs.

The pipeline is an adaptation of `cpg-seqera-trial-workflows/workflows/cohort_demo`, which provides the 
Hail script and the Dataproc submitter. Densify, VQSR and the compose steps are ported from
`cpg-flow-seqr-loader`. Everything around them follows the conventions of `nf-align`.

Terms are defined in [GLOSSARY.md](../GLOSSARY.md), [RISKS.md](RISKS.md) lists what may bite later.
[PROCESS.md](PROCESS.md) walks through every step in plain language, with the reasons behind the complex ones, and
[WORKFLOW.md](WORKFLOW.md) has diagrams of the same flow.

## Scope of a run

A run writes and registers one VDS in one dataset, for one sequencing type. Its cohorts may contain SGs from
multiple datasets, so long as the service account running the pipeline has permission to query metamist for 
all the included cohorts, and read access to the relevant buckets. Each of those datasets also gets its own VCF,
written to its own bucket, so the service account needs write access to those buckets as well.

| param                                                   | meaning                                                                                                                                                                                                                                                                     |
|---------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `dataset`                                               | Metamist project name, without the `-test` suffix. The VDS is written to and registered in this dataset                                                                                                                                                                     |
| `sequencing_type`                                       | `genome` or `exome`. Also selects which interval files are used                                                                                                                                                                                                             |
| `accessLevel`                                           | `full` / `standard` (equivalent), or `test`. Selects `<dataset>` or `<dataset>-test`, and the `-main` or `-test` bucket                                                                                                                                                     |
| `arGuid`                                                | analysis-runner GUID, used in resource labels, the run-record directory and the registered analysis. Mandatory                                                                                                                                                              |
| `cohorts`                                               | the cohorts whose SGs are combined. Mandatory, and the only input: there is no `samplesheet` or `sgids` input                                                                                                                                                               |
| `base_vds`                                              | optional. An explicit base VDS, overriding discovery                                                                                                                                                                                                                        |
| `force`                                                 | default `false`. Makes a fresh build: no base VDS is used and any saved combiner plan is discarded. Can't be combined with `base_vds`                                                                                                                                       |
| `refresh`                                               | optional list of SG IDs to refresh. Refreshing is the removal of gVCF content, followed by re-adding from the same source path. Typically done when the SG should be included, but a gVCF was truncated/malformed                                                           |
| `genome_combiner_intervals`, `exome_combiner_intervals` | the regions BED the combiner partitions on, chosen by `sequencing_type`. Both default to variant-balanced BEDs under `gs://cpg-common-main/rd_references/` (1000 intervals for genome, 500 for exome). The run refuses to start if the one for its sequencing type is unset |
| `genome_densify_intervals`, `exome_densify_intervals`   | the coarser BEDs (500 for genome, 100 for exome) the VDS is read on when densifying. They set the dense MT's partitions, and so the number of shards in each VCF export. Must be `gs://` paths                                                                              |
| `vqsr`                                                  | default `true`. If `false`, the run stops after densifying                                                                                                                                                                                                                  |
| `vqsr_fragments_per_job`, `annotate_fragments_per_job`  | how many shards each VQSR task, and each task copying VQSR's results onto the whole fragments, handles. Both default to 100                                                                                                                                                 |
| `snp_filter_level`, `indel_filter_level`, `vqsr_*`      | VQSR's truth-sensitivity filter levels (99.8 and 99.0), and its resource VCFs, as in `cpg-flow-seqr-loader`                                                                                                                                                                |
| `target_records`, `branch_factor`, `gvcf_batch_size`    | combiner tuning, with the same defaults for both sequencing types                                                                                                                                                                                                           |
| `register`                                              | default `true`. If `false`, the VDS is written but not registered                                                                                                                                                                                                           |
| `dataproc_*`                                            | cluster shape and lifecycle, as in `cohort_demo`. `dataproc_project` was removed as each cluster can only be started in the service account's own project                                                                                                                   |
| `dataproc_autoscaling_policy`                           | defaults to `cpg-dataproc-medium`, also available are `-light` and `-heavy`. Same names in every project                                                                                                                                                                    |
| `cpg_container`                                         | the docker container with required libraries/packages to run the python components & interact with dataproc API                                                                                                                                                             |
| `dataproc_service_account`                              | in this implementation always has the pattern `seqera-{accessLevel}@{dataset}.iam.gserviceaccount.com`                                                                                                                                                                      |
| `dataproc_staging_bucket`                               | staging and Temp bucket are the same: `gs://cpg-<dataset>-<main\|test>-tmp`                                                                                                                                                                                                 |

## Flow

```
--cohorts ─ BUILD_SAMPLESHEET ─┐
                               ▼
              sheet, base VDS, registered VDSs ─ planCombine ─┬─ VDS exists & registered ─▶ no combine ─┐
                                                              ├─ any check fails ─▶ fail, no cluster    │
                                                              └─ COMBINE (Dataproc, then register) ─────┤
                                                                                                        ▼
                                  run records ◀─ DENSIFY (Dataproc) ◀─ planDensify ◀─ VDS, new or existing
                                                     │
                                                     ▼
                         dense MT, sites-only VCF shards (two layouts), whole-VCF fragments
                                                     │ --vqsr
                                                     ▼
                         VQSR ─▶ VQSR sites-only VCF, annotated whole fragments ─ DATASET_VCFS ─▶ one VCF per dataset
```

Everything that can fail cheaply fails before any cluster starts.

### 1. Sample sheet (`BUILD_SAMPLESHEET`, `bin/build_samplesheet.py`)

- One GraphQL query, scoped to `project(name: <dataset>[-test])`. It fetches the cohorts' status, their SGs
  with each SG's active `gvcf` analyses, and the project's active `vds` analyses. Analysis status is ignored.
- Cohorts that are missing from the project, inactive or empty fail the run.
- An SG in several cohorts is kept once, and the overlap is logged. SGs of another sequencing type are
  skipped with a warning.
- Each SG must have exactly one distinct active gVCF path. Several analyses of the same path count as one.
  An SG with none, or with several different paths, fails the run, and every offender is listed in one
  error.
- Each row records the SG's participant's primary external ID (`participant`).
- Each row records the SG's own dataset (its project, without `-test`). It splits the VQSR-annotated VCF into
  one VCF per dataset (`DATASET_VCFS`).
- The sheet is named for the VDS it makes:
  `<hash>_<count>.tsv`. The hash is the first 38 hex characters of the SHA-256 of the sorted SG IDs, joined,
  then `refresh:` and the sorted refresh IDs if there are any. The hash reaches `main.nf` through the sheet's filename.
- Alongside the sheet it writes the discovered base VDS (the newest active `vds` analysis with the same
  `sequencing_type`, whatever its cohorts) and the paths of all registered VDSs.
- This stage executes a directory-cached existence check on all gVCF files and indexes. Could be implemented in Groovy
  but the functools caching is a massive speed improvement.

### 2. Pre-flight (`planCombine` in `main.nf`)

- Params: types and required values are checked by the `params` block in `main.nf`; `validateParams` checks
  allowed values, not both `force` and `base_vds`, and combiner intervals for the sequencing type.
- Each `refresh` ID must be in the sheet.
- The base VDS: none with `force`, else `base_vds`, else the discovered one. It must exist.
- The output: `<dataset bucket, or output_root>/vds/<sequencing_type>/<sheet name>.vds`.
- The output must not exist, unless it's registered:

  | exists | outcome                                                     |
  |--------|-------------------------------------------------------------|
  | no     | combine                                                     |
  | yes    | no combine, logged, no cluster; the VDS is still densified  |

### 3. Combine (`COMBINE`, `bin/submit_hail.py`, `bin/combiner.py`)

A submitter task creates the Dataproc cluster (named `combine-<ar-guid>-<random>`) and runs the Hail script
on it. When the job finishes, or the task receives SIGTERM, it tears the cluster down. The pipeline doesn't
use spot VMs at all, because Nextflow can't exempt a single process from spot. The task gets one retry, which
resumes from the combiner plan saved at `<temp bucket>/nf-combine/combiner_plan/<name>.json`, where `<name>`
is the VDS filename without `.vds`.

If the VDS is already complete (both `reference_data/_SUCCESS` and `variant_data/_SUCCESS` exist), the
submitter skips the cluster. That's the retry after a combine that succeeded and a registration that
failed. The Hail script itself refuses to start over a partial VDS.

The Hail script:

- removes SGs that are in the base VDS but not in the sheet
- removes refreshed SGs that are in the base VDS
- adds the sheet's SGs that the base lacks, which includes the refreshed ones, from their gVCFs
- writes a report of the base VDS used and the SGs added, removed and refreshed, to
  `<temp bucket>/nf-combine/combine_reports/<name>.json`

Refresh IDs that aren't in the base VDS, and all refresh IDs in a fresh build, produce a warning and are
treated as ordinary additions.

Removal comes first, and drops any alleles only the removed SGs carried. If there's nothing to add, the
result is written straight to the VDS path. Otherwise it's written to
`<temp bucket>/nf-combine/combiner_removal/<name>.vds` and becomes the combiner's base. Because the base has
changed, any saved plan is discarded. A retry that finds that VDS complete skips the removal and resumes the
plan.

### 4. Registration

Registration happens in the `COMBINE` task, once the submitter succeeds. The VDS is written directly to its final
path, so there's no asynchronous publish to wait for. It creates a `vds` analysis in
`<dataset>[-test]`, linked to the cohorts. Its `meta` holds `ar-guid`, `sequencing_type`, `base_vds`, and
`refresh` when that list isn't empty.

### 5. Densify (`DENSIFY`, `bin/submit_hail.py`, `bin/densify.py`)

Every VDS is densified, whether this run combined it or it already existed. This is a port of
`cpg-flow-seqr-loader`'s `CreateDenseMtFromVdsWithHail` stage. The task has the same shape as `COMBINE`. A submitter
task creates a cluster (named `densify-<ar-guid>-<random>`), with the same `dataproc_*` settings as the combine,
runs the Hail script on it, and tears it down. The task gets one retry. The helpers the script imports
(`hail_lib/nf_combine_hail`) are shipped to the cluster as pyfiles.

The Hail script:

- reads the VDS on the densify intervals, densifies it, keeps only sites where at least one SG is non-ref, and
  checkpoints the result
- computes GATK-style INFO fields and splits multiallelics. Each split row's allele-specific INFO is cut down to
  its own alt. Hail's split copies the whole site's arrays onto every row, and GATK reads the first value, so
  without this every row would get the first alt's values. No INFO field is pipe-delimited, as VQSR can't parse
  them.
- writes the dense MT, then exports it three times, each as a directory of shards in partition order:
  - sites-only, a header per shard: each shard is a valid VCF on its own, for VQSR's per-shard steps
  - sites-only, separate header: the header and the shards join into one VCF, to train VQSR's models
  - whole fragments, a header per shard, with genotypes: VQSR's results are copied onto these, and they're split
    into the dataset VCFs

Everything is written under `<temp bucket>/nf-combine/densify/<name>/`. Because it's named for the VDS, a retry,
or a later run on the same VDS, reuses whatever is complete. Each step is skipped if its output is complete. The
MT counts as complete once it has `_SUCCESS`, and an export once it has `shard-manifest.txt`, which is written last.
An export directory without a manifest is a partial write, and is removed before exporting again. If the MT and all
three manifests already exist, the submitter skips the cluster.

The dense MT isn't registered or published. Nothing downstream reads it: it's kept only to resume from, and the
temp bucket's lifecycle rules remove it. The VCF exports are the product. The task copies each export's manifest
into its work directory, so the VQSR steps can split their work across the shards.

Sites-only shards and whole fragments come from the same partitions. The nth whole fragment holds the same records
as the nth sites-only shard. Hail gives each export's shards different names, so later steps pair them by position.

### 6. VQSR and dataset VCFs

With `vqsr` on, and unless the VDS's VQSR VCF already exists, the `VQSR` subworkflow filters the sites-only
exports and copies the results onto the whole fragments. The `DATASET_VCFS` subworkflow then writes one VCF per
dataset in the sheet, to `<dataset bucket>/vcf/<sequencing_type>/<name>.vcf.bgz`, skipping any that already exist.
See [PROCESS.md](PROCESS.md) for the steps and the reasons behind them.

### 7. Run records

`combiner_samplesheets/<ar-guid>/` holds:

- `samplesheet.tsv`: as in `nf-align`
- `output.json`: the VDS, the base VDS, and the cluster and job (or `null` if no cluster was needed)
- `combine_report.json`: the Hail script's report
- `densify_output.json`: the VDS, the MT and the three export directories, with the cluster and job (or `null`
  if no cluster was needed)

## Testing

- `test/test_build_samplesheet.py`: pytest for `bin/build_samplesheet.py`, with the query validated against
  the schema bundled with the metamist client.
- `test/test_nextflow_stub.py`: `nextflow run -stub -profile test` checks of the wiring and the pre-flight
  checks. `BUILD_SAMPLESHEET`'s stub copies its outputs from `task.ext.fixtures`, set from config. One real run swaps fakes in on `PATH` to check the submitter and registration commands as rendered.
  The tests check that the VDS path follows the sheet's name; `vds_name` itself is unit tested, including
  against `cohort_demo`'s hash.

## Deliberately out of scope (for now)

- Fanning outputs other than the dataset VCFs out to each SG's dataset.
- Guarding against two concurrent runs on the same SG set.
- Cleaning up temp files and saved plans. The bucket's lifecycle rules handle them.
- Building a VDS from individual `sgids`, or from a prepared sample sheet.
