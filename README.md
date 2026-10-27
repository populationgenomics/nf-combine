# nf-combine

Nextflow/Seqera workflow for combining gVCFs into VDSs.

Combines the single-sample gVCFs of a set of Metamist cohorts into one Hail VDS. The combine runs on an
ephemeral Dataproc cluster, and builds on the dataset's newest earlier VDS. See [docs/DESIGN.md](docs/DESIGN.md).

```bash
nextflow run main.nf \
    --dataset <dataset> --accessLevel standard --sequencing_type genome --arGuid <guid> \
    --cohorts COH1,COH2 \
    --cpg_container <image> \
    --dataproc_service_account <sa> --dataproc_staging_bucket gs://<bucket> --dataproc_autoscaling_policy <policy>
```

For a test or one-off run, replace `--cohorts` with a prepared sheet (`sg_id`, `gvcf`, `dataset`) and the VDS
path to write: `--samplesheet sheet.tsv --output_vds gs://<bucket>/<path>.vds`. Nothing is registered.

Development:

```bash
uv sync
uv run pytest -rs        # stub-run tests need nextflow on PATH
uv run ruff check . && uv run ruff format --check .
nextflow lint main.nf modules
```
