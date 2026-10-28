# nf-combine: design brief

Combines the single-sample gVCFs of a set of Metamist cohorts into one Hail VDS. A run normally builds on
the newest earlier VDS for the dataset: it removes withdrawn SGs, refreshes changed ones, and adds new ones.
It's a Nextflow (≥ 25.10, strict syntax) pipeline. Small tasks run on Google Batch, and the combine itself
runs on an ephemeral Dataproc cluster. It can be tested offline with stub runs.

The pipeline adapts `cpg-seqera-trial-workflows/workflows/cohort_demo`, which provides the Hail script and
the Dataproc submitter. Everything around them follows the conventions of `nf-align`.

Terms are defined in [GLOSSARY.md](../GLOSSARY.md), which is shared across CPG workflows. Why things are
the way they are is recorded in [adr/](adr/). [RISKS.md](RISKS.md) lists what may bite later.

## Scope of a run

A run writes and registers one VDS in one dataset, for one sequencing type. Its cohorts may contain SGs from
other datasets ([ADR 0001](adr/0001-one-output-dataset-cohorts-may-span-datasets.md)).

| param | meaning |
|---|---|
| `dataset` | Metamist project name, without the `-test` suffix. The VDS is written to and registered in this dataset |
| `sequencing_type` | `genome` or `exome`. Also selects which interval files are used |
| `accessLevel` | `full` / `standard` (equivalent), or `test`. Selects `<dataset>` or `<dataset>-test`, and the `-main` or `-test` bucket |
| `arGuid` | analysis-runner GUID, used in resource labels, the run-record directory and the registered analysis. Mandatory |
| `cohorts` | the cohorts whose SGs are combined. Mandatory, and the only input: there is no `samplesheet` or `sgids` input |
| `base_vds` | optional. An explicit base VDS, overriding discovery |
| `force` | default `false`. Makes a fresh build: no base VDS is used and any saved combiner plan is discarded. Can't be combined with `base_vds` |
| `refresh` | optional list of SG IDs to refresh |
| `genome_combiner_intervals`, `exome_combiner_intervals` | the regions BED the combiner partitions on, chosen by `sequencing_type`. Both default to variant-balanced BEDs under `gs://cpg-common-main/rd_references/` (5000 intervals for genome, 500 for exome). The run refuses to start if the one for its sequencing type is unset |
| `genome_densify_intervals`, `exome_densify_intervals` | the coarser BEDs (2000 for genome, 100 for exome) for densifying the VDS. Declared with defaults, but not used yet: no stage densifies |
| `target_records`, `branch_factor`, `gvcf_batch_size` | combiner tuning, with the same defaults for both sequencing types |
| `register` | default `true`. If `false`, the VDS is written but not registered |
| `dataproc_*` | cluster shape and lifecycle, as in `cohort_demo`. `dataproc_project` defaults to the dataset |
| `dataproc_service_account`, `dataproc_staging_bucket`, `dataproc_autoscaling_policy`, `cpg_container` | no defaults. The run refuses to start without them; the `test` profile sets placeholders |

The Dataproc temp bucket is `gs://cpg-<dataset>-<main|test>-tmp`.

## Flow

```
--cohorts ─ BUILD_SAMPLESHEET ─┐
                               ▼
              sheet, base VDS, registered VDSs ─ planCombine ─┬─ VDS exists & registered ─▶ no-op, exit 0
                                                              ├─ any check fails ─▶ fail, no cluster
                                                              └─ COMBINE (Dataproc, then register) ─▶ run records
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
  error ([ADR 0002](adr/0002-fail-run-on-missing-or-ambiguous-gvcf.md)).
- Each row records the SG's own dataset (its project, without `-test`). Nothing uses it yet; later stages
  will.
- The sheet is named for the VDS it makes ([ADR 0004](adr/0004-vds-named-by-content-never-overwritten.md)):
  `<hash>_<count>.tsv`. The hash is the first 38 hex characters of the SHA-256 of the sorted SG IDs, joined,
  then `refresh:` and the sorted refresh IDs if there are any. Without a refresh this is exactly
  `cohort_demo`'s naming, and as in `cohort_demo` the name reaches `main.nf` through the sheet's filename.
- Alongside the sheet it writes the discovered base VDS (the newest active `vds` analysis with the same
  `sequencing_type`, whatever its cohorts) and the paths of all registered VDSs.

### 2. Pre-flight (`planCombine` in `main.nf`)

- Params: `cohorts` given, not both `force` and `base_vds`, and combiner intervals for the sequencing type.
- Each `refresh` ID must be in the sheet.
- Each gVCF and its `.tbi` must exist. Every missing file is listed in one error.
- The base VDS: none with `force`, else `base_vds`, else the discovered one. It must exist.
- The output: `<outputDir>/vds/<sequencing_type>/<sheet name>.vds`.
- The output must not exist, unless it's registered:

  | exists | registered | outcome |
  |---|---|---|
  | no | no | combine |
  | yes | yes | successful no-op, logged, no cluster |
  | yes | no | fail: "would create X, but X already exists and is not registered" |
  | no | yes | fail: registered in Metamist but missing |

### 3. Combine (`COMBINE`, `bin/submit_combiner.py`, `bin/combiner.py`)

A submitter task creates the Dataproc cluster (named `combine-<ar-guid>-<random>`) and runs the Hail script
on it. When the job finishes, or the task receives SIGTERM, it tears the cluster down. The pipeline doesn't
use spot VMs at all, because Nextflow can't exempt a single process from spot
([ADR 0003](adr/0003-hail-on-ephemeral-dataproc-via-non-spot-submitter.md)). The task gets one retry, which
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

### 5. Run records

`combiner_samplesheets/<ar-guid>/` holds:

- `samplesheet.tsv`: as in `nf-align`
- `output.json`: the VDS, the base VDS, and the cluster and job (or `null` if no cluster was needed)
- `combine_report.json`: the Hail script's report

## Testing

- `test/test_build_samplesheet.py`: pytest for `bin/build_samplesheet.py`, with the query validated against
  the schema bundled with the metamist client.
- `test/test_nextflow_stub.py`: `nextflow run -stub -profile test` checks of the wiring and the pre-flight
  checks. `BUILD_SAMPLESHEET`'s stub copies its outputs from `task.ext.fixtures`, set from config. One real run swaps fakes in on `PATH` to check the submitter and registration commands as rendered.
  The tests check that the VDS path follows the sheet's name; `vds_name` itself is unit tested, including
  against `cohort_demo`'s hash.
- A local Hail integration test on tiny gVCFs comes later, once the removal and refresh code exists.

## Deliberately out of scope (for now)

- Fanning outputs out to each SG's dataset.
- Guarding against two concurrent runs on the same SG set.
- Cleaning up temp files and saved plans. The bucket's lifecycle rules handle them.
- Building a VDS from individual `sgids`, or from a prepared sample sheet.
