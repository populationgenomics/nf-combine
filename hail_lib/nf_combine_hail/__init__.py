"""Hail helpers imported by the Hail scripts in bin/ when they run on a Dataproc cluster.

The submitter (bin/submit_hail.py) zips this package and ships it with each job as a pyfile, so the
scripts can `import nf_combine_hail.<module>`. Only the Spark driver imports it: Hail builds its
expressions on the driver, so executors never need this code.

annotations, sparse_mt and vcf are copied from cpg-flow-seqr-loader's hail_scripts, which had in turn
copied them from gnomad_methods. Only what densifying a VDS uses was kept. See hail_lib/README.md.
"""
