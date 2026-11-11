# How the nf-combine workflow runs

A run takes a set of Metamist cohorts and produces two things: one combined Hail VDS, and one VQSR-filtered VCF for
each dataset that has samples in those cohorts. Each step checks whether its output
already exists and skips the work if it does. A failed or repeated run therefore picks up where the last one stopped.
Output locations are generated based on the combined hash of all sequencing groups, so workflow resumption is reliable. 

Parameters and design decisions are in [DESIGN.md](DESIGN.md), and terms are defined in [GLOSSARY.md](../GLOSSARY.md).

---

## 1. Check the parameters

The run checks its settings before doing any work:

- The sequencing type is genome or exome.
- The access level is full, standard or test.
- `--force` and `--base_vds` aren't both set.
- The interval files for the sequencing type are set, and they are `gs://` paths.

**Why:** Dataproc clusters are the expensive part. Anything that can fail cheaply fails before a cluster starts.

## 2. Build the sample sheet (`BUILD_SAMPLESHEET`)

One Metamist query fetches the cohorts, their sequencing groups (SGs) and each SG's gVCF. The result is a table with
one row per SG: its ID, its gVCF path, the dataset it belongs to and its participant's primary external ID.

- Missing, inactive or empty cohorts stop the run.
- An SG that appears in several cohorts is kept once.
- SGs of the wrong sequencing type are skipped with a warning.
- Every SG must have exactly one gVCF. All SGs that don't are listed together in one error.
- Every gVCF and its index are checked to exist.

It also finds the base VDS: the newest registered VDS of the same sequencing type in the dataset.

**Why the sheet is named with a hash:** the filename is a hash of the sorted SG IDs (plus any refresh IDs), followed
by the SG count. The VDS takes the same name. The same set of SGs therefore always maps to the same VDS path, which is
how a later run knows the work is already done.

## 3. Plan the combine (`planCombine`)

This step works out the output path (`<dataset bucket>/vds/<seq type>/<name>.vds`) and which base VDS to build on.
With `--force` there is no base, a supplied `--base_vds` is used next, and the discovered one otherwise.

- If the VDS is already complete, the combine is skipped and the run goes straight to densifying.
- Otherwise, the run goes on to combine.

It also checks whether the VDS's VQSR VCF (step 6) and every dataset's VCF of the VDS (step 7) already exist. If
the dataset VCFs all do, the VDS isn't densified, and nothing after the combine runs for it except the seqr VCFs
(step 8). If the VQSR VCF does, the VDS isn't densified either, and the missing dataset VCFs are cut from it.

**Why it checks the `_SUCCESS` marker rather than the folder:** a half-written VDS still has its folder. Hail writes
the marker last, so the marker means the VDS is complete.

## 4. Combine (`COMBINE`, on a temporary Dataproc cluster)

A small task starts a Dataproc cluster, runs the Hail combiner (`bin/combiner.py`) on it, and shuts the cluster down.
The cluster is also shut down if the task is cancelled. Starting from the base VDS, the combiner:

1. removes SGs that are no longer in the cohorts
2. removes SGs listed in `--refresh`, so they can be re-added from their current gVCF
3. adds every SG the base doesn't have, including the refreshed ones

It writes a report of what was added, removed and refreshed. If registration is on, it then records the VDS in
Metamist.

**Why it's built this way:**

- **Removals first.** Removing SGs also drops alleles that only those SGs carried. The cleaned-up VDS then becomes the
  base for the additions.
- **Plans saved by name.** Hail's combiner saves a plan as it works, stored under the VDS name. The single retry
  resumes from that plan instead of starting again.
- **No cluster if the VDS exists.** If the VDS is already complete, the cluster isn't started at all. This covers the
  case where the combine succeeded but registration failed.
- **No spot VMs anywhere.** Nextflow can't exclude a single process from spot instances, so the whole pipeline avoids
  them.
- **Uploads before the cluster starts.** Everything the job needs is uploaded first, so a wrong path fails before any
  cluster cost.

## 5. Densify (`DENSIFY`, on a second temporary cluster)

Every VDS is densified, whether this run made it or it already existed, unless its VQSR VCF, or every dataset's VCF
of it, already exists.

*Why the check is here, not at VQSR:* densifying and VQSR exist only to make the VQSR VCF, and it exists only to make
the dataset VCFs. Once either exists there's nothing left for densifying or VQSR to do. Checking here also avoids
starting a cluster to re-densify a VDS whose densified outputs have been removed from the temp bucket.

**What is a VDS?** A compact format that stores reference blocks and variant calls separately. Most tools can't read
it directly.

This step (`bin/densify.py`):

1. Converts the VDS into a regular sample-by-site matrix (a "dense" MatrixTable).
2. Keeps only sites where at least one sample has a non-reference allele.
3. Computes GATK-style quality annotations, which VQSR needs.
4. Splits sites with several alternate alleles into one row per allele. Each row keeps only its own allele's
   allele-specific values.
5. Exports three sets of VCF files:
   - sites-only shards, each a complete VCF with its own header, which can be processed independently
   - sites-only shards with a single separate header, designed to be joined into one file
   - "whole" shards that include every sample's genotypes, each with its own header

Each export writes a list of its shards (a manifest) last. Later steps use the manifests to decide how to split up
their work.

**Why it's built this way:**

- **Trimming per-allele values.** Hail's multi-allelic split copies all of a site's per-allele values onto every split
  row. GATK reads the first value in each row, so without trimming, every split row would get the first allele's
  numbers.
- **No pipe-delimited fields.** VQSR can't parse them, so they're written in another form.
- **Resumable exports.** An export folder without a manifest is treated as half-written, deleted, and redone.
  Re-running continues from the last finished output.
- **Intermediate files in the temp bucket.** They're kept only so a run can resume, and the bucket's lifecycle rules
  clean them up.

## 6. VQSR filtering (`VQSR` subworkflow)

This step runs on every densified VDS when `--vqsr` is set. VQSR is GATK's
machine-learning filter. It learns what real variants look like from trusted reference sets, then scores and filters
every site. Its result is the **VQSR VCF**: one genome-wide VCF of every sample in the run, with VQSR's filters and
annotations, at `<temp bucket>/nf-combine/vqsr/<name>.vcf.bgz` (`<output_root>/tmp/nf-combine/...` for local
runs).

*Why it's kept, and why in the temp bucket:* it's the point a run resumes from. A later run that finds it skips
densifying and VQSR, and only cuts the missing dataset VCFs from it. It's large, so it's left to the temp bucket's
lifecycle rules. Once they remove it, a missing dataset VCF means densifying again.

1. **Join the sites-only shards (`COMPOSE_SITES_ONLY`).** The separate-header shards are joined into one genome-wide
   file using Google Cloud Storage's compose operation, so nothing is downloaded.
   *Why:* the models have to be trained on the whole callset. Compose joins at most 32 files at a time, so it runs in
   rounds. This works because each shard is already a complete compressed block.
2. **Train the models.** One model is trained for indels and one for SNPs, each on the joined file.
3. **Score SNPs in chunks (`TRAIN_VQSR_SNP_TRANCHES`).** The header-per-shard files are grouped into chunks of
   `vqsr_fragments_per_job`, and each chunk is scored against the SNP model.
   *Why chunks:* with one task per shard, starting the container and copying the reference files cost almost all of
   the time. Chunking spreads that cost over many shards.
4. **Gather the SNP results (`GATHER_VQSR_SNP_TRANCHES`).** Every chunk's scoring thresholds are merged into one
   genome-wide set.
5. **Apply the filters (`APPLY_VQSR`).** For each chunk, the SNP and then the indel filters are applied shard by
   shard.
   *Why per shard works:* the indel model covers the whole genome, and GATK looks up results by genomic position.
   Applying it to one shard at a time gives the same result as applying it to the whole file.
6. **Copy the results onto the genotypes (`ANNOTATE_FRAGMENTS`).** Each whole shard receives the FILTER values and
   VQSR annotations from the matching sites-only shard.
   *Why it's more involved than it sounds:*
   - The whole shards are grouped into chunks of their own size (`annotate_fragments_per_job`). Each shard's filtered
     file is regrouped to match.
   - Shards are paired by their position in the export, not by name, because the names don't match (Hail incorporates a hash value into each fragment name). The shard counts
     must be equal, and each pair's record counts are compared before annotating.
   - Each whole shard is annotated from its own small filtered file rather than one joined across shards. Without a
     region, bcftools would read a joined file from the start of the chromosome every time.
   - Each annotated shard is written without its header, ready to join. The header is written once, as its own small
     file. An empty shard is an empty file.
7. **Join the VQSR VCF (`COMPOSE_FRAGMENTS`).** The header, the annotated shards in genomic order and an end-of-file
   marker are joined in cloud storage, so nothing is downloaded. The result is then copied to the VQSR VCF's path.
   *Why the end-of-file markers are removed:* bcftools ends every compressed file with an end-of-file marker. GATK and
   htsjdk stop reading when they reach one, even in the middle of a file. The marker is removed from every piece, then
   added once at the very end.

## 7. One VCF per dataset (`SUBSET_DATASET_VCF`)

A plan is made for each dataset in the sample sheet. It covers that dataset's SGs and the output path
`<dataset bucket>/joint_vcf/<seq type>/<name>.vcf.bgz`. Datasets whose VCF already exists are skipped. Each missing
one is cut from the VQSR VCF, whether this run made it or an earlier one did. Without `--vqsr`, no dataset VCFs are
made.

1. **Cut down the VQSR VCF (`SUBSET_DATASET_VCF`).** The task copies the VQSR VCF to a scratch directory on its own
   disk, keeps only the dataset's samples and the sites where at least one of them has a non-reference allele, checks
   the samples are exactly the dataset's, indexes the result, and uploads the index and then the VCF. The scratch
   directory is deleted when the task ends, whether it succeeds or not.
   *Why the counts aren't recalculated:* AC and AN are left as computed across the whole callset (bcftools'
   `--no-update` option). Because each row holds a single alternate allele, a row is simply kept or dropped, and no
   allele trimming is needed.
   *Why the VCF goes last:* each upload appears in one step, so if the VCF exists, its index does too. A dataset VCF
   from before they were indexed still counts as done.
2. **Register the dataset's VCF (`REGISTER_DATASET_VCF`).** If registration is on, the VCF is recorded in the
   dataset's own Metamist project as a `dataset_vcf` analysis, linked to every SG in it, with the sequencing type and
   ar-guid as metadata. A failed registration is retried once.
   *Why only VCFs this run made:* a VCF that already exists was registered by the run that made it, so registering it
   again would duplicate the record.

## 8. seqr VCFs (`REHEADER_VCF`)

Only for the datasets listed in `seqr_vcf`. Each one's VCF is written a second time, beside the first, as
`<dataset bucket>/joint_vcf/<seq type>/<name>_seqr.vcf.bgz`. It has the same rows, but each sample is named by its
participant's primary external ID (the sample sheet's `participant` column) rather than its SG ID. It's built from
the dataset's VCF whether this run made it or an earlier one did. Datasets whose seqr VCF index already exists are
skipped.

The task copies the dataset's VCF to a scratch directory on its own disk, renames the samples with
`bcftools reheader`, checks that every sample was renamed, indexes the result, and uploads the VCF and then its index.
The scratch directory is deleted when the task ends, whether it succeeds or not.

*Why scratch, not the work directory:* anything in the work directory is kept, so a copy of every dataset VCF would
pile up there. *Why the index goes last:* if it exists, the VCF upload before it finished.

*Why it's checked before combining:* a seqr VCF needs every SG in the dataset to have a participant ID, and no two
SGs to share one, or two samples would have the same name. The run checks this as soon as the sample sheet is built,
before any cluster starts. A dataset listed in `seqr_vcf` with no SGs in the sheet only gets a warning.

---

## Run records

Each run writes `combiner_samplesheets/<ar-guid>/` in the dataset bucket. It holds:

- the sample sheet
- the cluster and job details
- the combine report (what was added, removed and refreshed)
- the densify record
