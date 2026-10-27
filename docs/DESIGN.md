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
| `sequencing_type` | `genome` or `exome`. Also selects the default `intervals` file |
| `accessLevel` | `full` / `standard` (equivalent), or `test`. Selects `<dataset>` or `<dataset>-test`, and the `-main` or `-test` bucket |
| `arGuid` | analysis-runner GUID, used in resource labels, the run-record directory and the registered analysis. Mandatory |
| `cohorts` **or** `samplesheet` + `output_vds` | exactly one input. A sample sheet is an alternative input with its own output: `output_vds` is required with `samplesheet` and refused without it. There is no `sgids` input |
| `base_vds` | optional. An explicit base VDS, overriding discovery |
| `force` | default `false`. Makes a fresh build: no base VDS is used and any saved combiner plan is discarded. Can't be combined with `base_vds` |
| `refresh` | optional list of SG IDs to refresh |
| `genome_intervals`, `exome_intervals` | the regions BED the combiner partitions on, chosen by `sequencing_type`. `exome_intervals` has no default yet |
| `target_records`, `branch_factor`, `gvcf_batch_size` | combiner tuning, with the same defaults for both sequencing types |
| `register` | default `true`. Only cohort runs register, whatever this says |
| `dataproc_*` | cluster shape and lifecycle, as in `cohort_demo`. `dataproc_project` defaults to the dataset |
| `dataproc_service_account`, `dataproc_staging_bucket`, `dataproc_autoscaling_policy`, `cpg_container` | no defaults. The run refuses to start without them; the `test` profile sets placeholders |

The Dataproc temp bucket is `gs://cpg-<dataset>-<main|test>-tmp`.

## Flow

```
--cohorts ─ BUILD_SAMPLESHEET ─┐            (or --samplesheet)
                               ▼
              sheet, base VDS, registered VDSs ─ planCombine ─┬─ VDS exists & registered ─▶ no-op, exit 0
                                                              ├─ any check fails ─▶ fail, no cluster
                                                              └─ COMBINE (Dataproc, then register) ─▶ run records
```

Everything that can fail cheaply fails before any cluster starts.

### 1. Sample sheet (`BUILD_SAMPLESHEET`, `bin/build_samplesheet.py`, cohort runs only)

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

A `--samplesheet` run skips this step, and writes to the VDS path given by `--output_vds` instead. The sheet
needs the columns `sg_id`, `gvcf` and `dataset`, and its filename doesn't matter. Metamist isn't consulted at
all, so there's no discovered base VDS and no VDS counts as registered. Everything from the pre-flight on is
the same for both inputs.

### 2. Pre-flight (`planCombine` in `main.nf`)

- Params: exactly one of `cohorts` / `samplesheet`, `output_vds` (ending `.vds`) exactly when
  `samplesheet` is given, not both `force` and `base_vds`, and intervals for the sequencing type.
- The sheet: the three columns, no empty values, no repeated SG, at least one row.
- Each `refresh` ID must be in the sheet.
- Each gVCF and its `.tbi` must exist. Every missing file is listed in one error.
- The base VDS: none with `force`, else `base_vds`, else the discovered one. It must exist.
- The output: `<outputDir>/vds/<sequencing_type>/<sheet name>.vds` for a cohort run, or `output_vds`.
- The output must not exist, unless it's registered:

  | exists | registered | outcome |
  |---|---|---|
  | no | no | combine |
  | yes | yes | successful no-op, logged, no cluster |
  | yes | no | fail: "would create X, but X already exists and is not registered" |
  | no | yes | fail: registered in Metamist but missing |

  A sample-sheet run never sees a VDS as registered, so it fails on any existing output.

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
treated as ordinary additions. Removing SGs from a VDS (`remove_sgs` in `bin/combiner.py`) is not written
yet: any run that needs it fails with `NotImplementedError` and a list of the SGs.

### 4. Registration

Registration happens in the `COMBINE` task, once the submitter succeeds. The VDS is written directly to its final
path, so there's no asynchronous publish to wait for. A cohort run creates a `vds` analysis in
`<dataset>[-test]`, linked to the cohorts. Its `meta` holds `ar-guid`, `sequencing_type`, `base_vds`, and
`refresh` when that list isn't empty. Sample-sheet runs don't register.

### 5. Run records

`combiner_samplesheets/<ar-guid>/` holds:

- `samplesheet.tsv`: cohort runs only, as in `nf-align`; a supplied sheet isn't republished
- `output.json`: the VDS, the base VDS, and the cluster and job (or `null` if no cluster was needed)
- `combine_report.json`: the Hail script's report

## Testing

- `test/test_build_samplesheet.py`: pytest for `bin/build_samplesheet.py`, with the query validated against
  the schema bundled with the metamist client.
- `test/test_nextflow_stub.py`: `nextflow run -stub -profile test` checks of the wiring and the pre-flight
  checks. For cohort runs, `BUILD_SAMPLESHEET`'s stub copies its outputs from `task.ext.fixtures`, set from
  config. One real run swaps fakes in on `PATH` to check the submitter and registration commands as rendered.
  Cohort runs check that the VDS path follows the sheet's name; `vds_name` itself is unit tested, including
  against `cohort_demo`'s hash.
- A local Hail integration test on tiny gVCFs comes later, once the removal and refresh code exists.

## Deliberately out of scope (for now)

- Fanning outputs out to each SG's dataset.
- Guarding against two concurrent runs on the same SG set.
- Cleaning up temp files and saved plans. The bucket's lifecycle rules handle them.
- Building a VDS from individual `sgids`.
- Registering sample-sheet runs.
