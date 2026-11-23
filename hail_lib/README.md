# hail_lib

`nf_combine_hail` is the Python package that the Hail scripts in `bin/` import on the Dataproc cluster.
Processes take it as a `path` input. `bin/submit_hail.py` then zips it and passes it to
`HailDataprocCluster.upload_pyfiles`, and Dataproc puts it on the job's `PYTHONPATH` (`python_file_uris`).

| module         | contents                                                                                  |
|----------------|-------------------------------------------------------------------------------------------|
| `intervals`    | `read_bed_file_as_intervals`, used by the combiner and by densifying                      |
| `sparse_mt`    | `default_compute_info`: GATK-style site and allele-specific INFO, ACs and lowqual flags    |
| `annotations`  | FS, SOR, allele-balance, lowqual and adj expressions used by `sparse_mt`                  |
| `vcf`          | `adjust_vcf_incompatible_types`: casts and pipe-delimits INFO fields for VCF export       |

`sparse_mt`, `annotations` and `vcf` are copied from
[cpg-flow-seqr-loader](https://github.com/populationgenomics/cpg-flow-seqr-loader)'s `hail_scripts`, which had
in turn copied them from [gnomad_methods](https://github.com/broadinstitute/gnomad_methods). Only what
densifying uses was kept. The logic is unchanged, and the docstrings were rewritten in Google style.

Hail isn't installed locally, so nothing here is unit tested yet.
