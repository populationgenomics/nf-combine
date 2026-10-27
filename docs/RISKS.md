# Risks: things that may come back to bite

Each entry gives the risk, why it matters, and what currently guards against it. None of these has been
verified by experiment yet. They're inferred from the design session and from `cohort_demo`.

## Metamist

1. **gVCF analyses for SGs from other datasets.** The query is scoped to the run's project, but cohort SGs
   can belong to other projects. Nobody has yet checked that `sequencingGroups { analyses }` returns those
   SGs' gVCF analyses, rather than filtering them by the queried project. If it filters them, those SGs
   look as if they have no gVCF, and the run fails loudly
   ([ADR 0002](adr/0002-fail-run-on-missing-or-ambiguous-gvcf.md)). It doesn't produce a wrong VDS.
2. **Base VDS discovery ignores cohorts.** The newest active `vds` of the right sequencing type in the
   dataset is used, whatever cohorts it's linked to. A dataset with several unrelated VDSs could pick a
   distant base. This is still correct, because the Hail script removes SGs that weren't requested, but it
   could be slow. Pass `base_vds` explicitly to choose a different one.
3. **A gVCF re-called at the same SG ID is invisible.** The VDS doesn't record which gVCF each SG came
   from. The only remedies are `refresh` and `force`.

## Runs

4. **Concurrent runs on the same SG set.** Both runs pass the pre-flight check, share the saved combiner
   plan and temp path, and race to write the same VDS. Nothing guards against this. It's considered
   unlikely, because runs are launched by hand.
5. **A VDS written but not registered.** If the combine succeeds and registration fails on both attempts,
   the VDS exists without an analysis. Every later run on the same SG set then fails at pre-flight by
   design ([ADR 0004](adr/0004-vds-named-by-content-never-overwritten.md)). A person has to delete the VDS
   or register it, carrying the run's metadata by hand.
6. **Preemption of the submitter.** Spot is off for the whole pipeline
   ([ADR 0003](adr/0003-hail-on-ephemeral-dataproc-via-non-spot-submitter.md)). If a launch config turns
   `google.batch.spot` back on, a preemption tears the cluster down. The retry resumes from the saved plan,
   so the cost is time rather than correctness.
7. **Sample-sheet runs share plans by filename.** The combiner plan and temp files are keyed by the VDS
   filename only. Two `--output_vds` paths with the same filename in different directories, such as
   `a/test.vds` and `b/test.vds`, would resume each other's plan. Content-named cohort VDSs can't collide this
   way. Give one-off VDSs distinct names, or use `--force`.
8. **Combiner tuning is genome-derived.** The defaults for `target_records`, `branch_factor` and
   `gvcf_batch_size` are shared by exome and genome, and haven't been tuned on exome data.

## Implementation assumptions (not yet checked against GCS or Hail)

9. **`file(<vds>).exists()` on GCS.** The pre-flight check treats a VDS path as existing if Nextflow says
   so. A VDS is a prefix with no object of its own, so this depends on Nextflow's GCS filesystem treating a
   prefix as a directory. If it doesn't, an existing VDS looks absent. A registered one then fails the run
   as "registered but missing", which is safe but wrong. An unregistered complete one goes on to `COMBINE`,
   where the submitter skips the cluster and registers it, against ADR 0004. Check this on a real bucket
   before the first production run.
10. **`_SUCCESS` markers.** The submitter takes a VDS as complete when both `reference_data/_SUCCESS` and
    `variant_data/_SUCCESS` exist. That relies on Hail writing `_SUCCESS` at the root of each matrix table.
11. **gVCF existence checks run one by one on the head node.** For thousands of SGs that's thousands of
    sequential GCS requests during pre-flight. It's slow but bounded. Make it parallel if it hurts.
12. **Cross-dataset gVCFs (risk 1).** In the bundled Metamist schema, `GraphQLSequencingGroup.analyses` has
    an optional `project` filter that defaults to null. That suggests other projects' analyses are returned,
    but it hasn't been checked against a live query.

## Inherited from `cohort_demo`

13. `dataproc_autoscaling_policy` is used but was never declared in `cohort_demo`'s config, and its
    `COMPARE_VDS` uses an undeclared `params.gcp_project`. Don't copy these across unchanged.
