# Workflow diagrams

How a run moves from Metamist cohorts to one VCF per dataset. Boxes are processes, rounded boxes are files or
records, and diamonds are the checks that let a re-run skip finished work. See [DESIGN.md](DESIGN.md) for the reasons
behind each step.

## Overview

```mermaid
flowchart TD
    params[/"--dataset, --cohorts, --refresh,<br/>--force, --base_vds"/]
    metamist[("Metamist")]

    params --> BUILD
    metamist --> BUILD
    BUILD["BUILD_SAMPLESHEET<br/><i>bin/build_samplesheet.py</i>"]
    BUILD --> sheet(["#lt;hash#gt;_#lt;n SGs#gt;.tsv<br/>SG, gVCF, dataset"])
    BUILD --> base(["base_vds.txt<br/>newest earlier VDS"])

    sheet --> planC{{"planCombine()<br/>VDS _SUCCESS exists?"}}
    base --> planC
    planC -- "no" --> COMBINE
    planC -- "yes: already combined" --> vds

    subgraph dp1 ["Ephemeral Dataproc cluster"]
        COMBINE["COMBINE<br/><i>bin/combiner.py</i><br/>gVCFs + base VDS → VDS"]
    end
    COMBINE --> vds[("#lt;dataset bucket#gt;/vds/#lt;seq type#gt;/#lt;name#gt;.vds")]
    COMBINE -. "--register" .-> metamist

    vds --> planD{{"planDensify()"}}
    planD --> DENSIFY
    subgraph dp2 ["Ephemeral Dataproc cluster"]
        DENSIFY["DENSIFY<br/><i>bin/densify.py</i><br/>VDS → dense MT → VCF exports"]
    end
    DENSIFY --> so(["sites-only,<br/>header per shard"])
    DENSIFY --> sh(["sites-only,<br/>separate header"])
    DENSIFY --> wf(["whole fragments<br/>(with genotypes)"])

    so --> planV{{"planVqsr()<br/>--vqsr set, and<br/>VQSR VCF .tbi missing?"}}
    sh --> planV
    wf --> planV
    planV -- "yes" --> VQSR[["VQSR subworkflow"]]
    planV -- "no" --> skipV(["skipped"])

    VQSR --> vqsrvcf(["#lt;name#gt;_vqsr.vcf.gz<br/>genome-wide sites-only"])
    VQSR --> annot(["annotated whole fragments<br/>+ shared header"])

    sheet --> planDS{{"planDatasetVcfs()<br/>one plan per dataset in the sheet;<br/>final VCF missing?"}}
    annot --> DSV[["DATASET_VCFS subworkflow"]]
    planDS -- "yes" --> DSV
    DSV --> dsvcf[("#lt;each dataset bucket#gt;/vcf/#lt;seq type#gt;/#lt;name#gt;.vcf.bgz")]
```

## VQSR subworkflow

Two chunkings of the same partitions run side by side. Sites-only shards are chunked by `vqsr_fragments_per_job`
and whole fragments by `annotate_fragments_per_job`. They are paired by their position in the export, not by name.

```mermaid
flowchart TD
    sh(["separate-header sites-only shards"]) --> COMPOSE["COMPOSE_SITES_ONLY<br/>gcloud storage objects compose"]
    COMPOSE --> composed(["one genome-wide sites-only VCF"])

    composed --> TIM["TRAIN_VQSR_INDEL_MODEL"]
    composed --> TSM["TRAIN_VQSR_SNP_MODEL"]
    TIM --> indel(["indel recal + tranches"])
    TSM --> snpmodel(["snp.model"])

    so(["header-per-shard sites-only shards"]) --> chunks["chunk by<br/>vqsr_fragments_per_job"]
    chunks --> TST["TRAIN_VQSR_SNP_TRANCHES<br/>× chunk"]
    snpmodel --> TST
    TST -- "per-shard tranches" --> GATHER["GATHER_VQSR_SNP_TRANCHES"]
    TST -- "per-shard recal" --> APPLY
    GATHER --> APPLY
    indel --> APPLY
    chunks --> APPLY["APPLY_VQSR<br/>× chunk<br/>SNP then indel"]

    APPLY -- "chunk VCFs, in genomic order" --> CONCAT["CONCAT_VCFS"]
    CONCAT --> vqsrvcf(["#lt;name#gt;_vqsr.vcf.gz"])

    APPLY -- "per-shard VQSR VCFs" --> regroup["regroup by position into<br/>annotate_fragments_per_job chunks"]
    wf(["whole fragments"]) --> wchunks["chunk by<br/>annotate_fragments_per_job<br/>(count must match shards)"]
    wchunks --> ANNOT["ANNOTATE_FRAGMENTS<br/>× chunk<br/>bcftools annotate FILTER + VQSR INFO"]
    regroup --> ANNOT
    ANNOT --> frags(["annotated whole fragments"])
    ANNOT -- "chunk 0's header" --> header(["header.vcf"])
```

## DATASET_VCFS subworkflow

Every annotated chunk is subset for every dataset as soon as it's ready. Each dataset's VCF is then joined in GCS
without being downloaded.

```mermaid
flowchart TD
    header(["annotated header"]) --> SH["SUBSET_HEADER<br/>× dataset"]
    frags(["annotated whole fragments, by chunk"]) --> SF["SUBSET_FRAGMENTS<br/>× chunk × dataset<br/>bcftools view -s SGs --min-ac 1"]
    plans(["dataset plans<br/>(SGs per dataset)"]) --> SH
    plans --> SF

    SH --> dh(["dataset header.vcf.bgz + eof.bgz"])
    SF --> sub(["headless subset fragments"])

    dh --> order["order: header, fragments by chunk then<br/>part name, EOF block"]
    sub --> order
    order --> CF["COMPOSE_FRAGMENTS<br/>gcloud storage objects compose"]
    CF --> out[("#lt;dataset bucket#gt;/vcf/#lt;seq type#gt;/#lt;name#gt;.vcf.bgz")]
```
