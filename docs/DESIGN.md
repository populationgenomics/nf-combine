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
| `cohorts` **or** `samplesheet` | exactly one. There is no `sgids` input |
| `base_vds` | optional. An explicit base VDS, overriding discovery |
| `force` | default `false`. Makes a fresh build: no base VDS is used and any saved combiner plan is discarded. Can't be combined with `base_vds` |
| `refresh` | optional list of SG IDs to refresh |
| `intervals` | defaults by `sequencing_type` |
| `target_records`, `branch_factor`, `gvcf_batch_size` | combiner tuning, with the same defaults for both sequencing types |
| `dataproc_*` | cluster shape and lifecycle, as in `cohort_demo` |
| `dataproc_service_account`, `dataproc_staging_bucket`, `dataproc_autoscaling_policy`, `cpg_container` | no defaults. The `standard` profile refuses to run without them |

The Dataproc temp bucket is `gs://cpg-<dataset>-<main|test>-tmp`.

## Flow

```
--cohorts | --samplesheet ─ BUILD_SAMPLESHEET ─┬─ VDS exists & registered ─▶ no-op, exit 0
                                               ├─ VDS exists, unregistered ─▶ fail
                                               └─ COMBINE (Dataproc) ─ register ─▶ run records
```

### 1. Sample sheet and pre-flight (`BUILD_SAMPLESHEET`)

Everything that can fail cheaply fails here, before any cluster starts.

- Check params: exactly one of `cohorts` / `samplesheet`, and not both `force` and `base_vds`.
- `--cohorts`: one GraphQL query, scoped to `project(name: <dataset>[-test])`, fetches the cohorts' SGs and
  their active `gvcf` analyses. Analysis status is ignored. An SG in several cohorts is kept once, and the
  overlap is logged.
- Each SG must have exactly one distinct active gVCF path. Several analyses of the same path count as one.
  An SG with none, or with several different paths, fails the run, and every offender is listed in one
  error ([ADR 0002](adr/0002-fail-run-on-missing-or-ambiguous-gvcf.md)).
- `--samplesheet`: a supplied sheet with the columns `sg_id`, `gvcf` and `dataset`. Metamist is not queried.
- Each gVCF and its `.tbi` must exist. Every missing file is listed in one error.
- Each `refresh` ID must be in the requested set, or the run fails.
- Choose the base VDS:
  - `force`: none.
  - `base_vds`: that VDS.
  - Cohort run: the newest active `vds` analysis in the dataset with the same `sequencing_type`.
  - Sample-sheet run: none.
- Name the output ([ADR 0004](adr/0004-vds-named-by-content-never-overwritten.md)). It's
  `gs://cpg-<dataset>-<main|test>/vds/<sequencing_type>/<hash>_<count>.vds`, where the hash is taken over
  the sorted SG IDs and the sorted refresh IDs.
- If that path exists and an active `vds` analysis points at it, end the run as a successful no-op. If it
  exists but isn't registered, fail with "would create X, but X already exists and is not registered".

The sample sheet's `dataset` column is recorded only. Later stages will use it.

### 2. Combine (`COMBINE`)

A submitter task on a non-spot VM creates the Dataproc cluster and runs the Hail script on it. When the job
finishes, or the task receives SIGTERM, it tears the cluster down
([ADR 0003](adr/0003-hail-on-ephemeral-dataproc-via-non-spot-submitter.md)). It gets one retry, which
resumes from the combiner plan saved under the temp bucket and keyed by the output hash.

The Hail script:

- removes SGs that are in the base VDS but not in the request
- removes refreshed SGs
- adds requested SGs that aren't in the base VDS, plus the refreshed SGs, from their gVCFs

Refresh IDs that aren't in the base VDS, and all refresh IDs in a fresh build, produce a warning and are
treated as ordinary additions. The removal and refresh logic is still to be written.

### 3. Registration

Registration happens inside the submitter, once the job succeeds. The VDS is written directly to its final
path, so there's no asynchronous publish to wait for. A cohort run creates a `vds` analysis in
`<dataset>[-test]`, linked to the cohorts. Its `meta` holds `ar-guid`, `sequencing_type`, `base_vds`, and
`refresh` when that list isn't empty. Sample-sheet runs don't register.

### 4. Run records

`combiner_samplesheets/<ar-guid>/` holds the sample sheet, the base VDS used, the SGs added and removed, and
the cluster's `output.json`.

## Testing

- Stub runs of the Nextflow wiring, under the `test` profile, which turns registration off.
- pytest for the pure Python logic: turning query results into a sample sheet, de-duplicating SGs, the gVCF
  rules, choosing the base VDS, the output hash, and checking refresh IDs.
- A local Hail integration test on tiny gVCFs comes later, once the removal and refresh code exists.

## Deliberately out of scope (for now)

- Fanning outputs out to each SG's dataset.
- Guarding against two concurrent runs on the same SG set.
- Cleaning up temp files and saved plans. The bucket's lifecycle rules handle them.
- Building a VDS from individual `sgids`.
- Registering sample-sheet runs.
