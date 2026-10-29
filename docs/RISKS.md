# Risks: things that may come back to bite

Each entry gives the risk, why it matters, and what currently guards against it. None of these has been
verified by experiment yet.

## Metamist

1. **Base VDS discovery ignores cohorts.** The newest active `vds` of the right sequencing type in the
   dataset is used, whatever cohorts it's linked to. A dataset with several unrelated VDSs could pick a
   distant base. This is still correct, because the Hail script removes SGs that weren't requested, but it
   could be slow. Pass `base_vds` explicitly to choose a different one.
2. **A gVCF re-called at the same SG ID is invisible.** The VDS doesn't record which gVCF each SG came
   from. The only remedies are `refresh` and `force`. More explicitly in

## Runs

3. **Concurrent runs on the same SG set.** Both runs pass the pre-flight check, share the saved combiner
   plan and temp path, and race to write the same VDS. Nothing guards against this. It's considered
   unlikely, because runs are launched by hand.
4. **A VDS written but not registered.** If the combine succeeds and registration fails on both attempts,
   the VDS exists without an analysis. Every later run on the same SG set then fails at pre-flight by
   design ([ADR 0004](adr/0004-vds-named-by-content-never-overwritten.md)). A person has to delete the VDS
   or register it, carrying the run's metadata by hand.
5. **Preemption of the submitter.** Spot is off for the whole pipeline
   ([ADR 0003](adr/0003-hail-on-ephemeral-dataproc-via-non-spot-submitter.md)). If a launch config turns
   `google.batch.spot` back on, a preemption tears the cluster down. The retry resumes from the saved plan,
   so the cost is time rather than correctness.
6. **Combiner tuning is genome-derived.** The defaults for `target_records`, `branch_factor` and
   `gvcf_batch_size` are shared by exome and genome, and haven't been tuned on exome data.

## Implementation assumptions (not yet checked against GCS or Hail)

7. **`_SUCCESS` markers.** The submitter takes a VDS as complete when both `reference_data/_SUCCESS` and
    `variant_data/_SUCCESS` exist. That relies on Hail writing `_SUCCESS` at the root of each matrix table.
8. **Cross-dataset gVCFs (risk 1).** In the bundled Metamist schema, `GraphQLSequencingGroup.analyses` has
    an optional `project` filter that defaults to null. That suggests other projects' analyses are returned,
    but it hasn't been checked against a live query.
