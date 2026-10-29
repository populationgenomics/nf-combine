# A VDS is named by its content and never overwritten

The output is `<hash>_<count>.vds`. The hash is taken over the sorted SG IDs plus the sorted refresh IDs,
and the count is the number of SGs. Refresh IDs are part of the hash because a refresh changes a VDS
without changing its SG set. Without them, a refresh would get the base VDS's name and become a no-op.

As in `cohort_demo`, the name is made where the SGs are fetched, in `build_samplesheet.py`, and reaches the
workflow as the sample sheet's filename. That keeps it in Python, where it's unit tested. Every run starts from
cohorts, so every VDS is named this way: there is no input for a prepared sample sheet or an explicit output path.

An existing VDS is never overwritten, even with `force`. If the target path already exists, the first
stage stops before any cluster starts:

- If the VDS is registered in Metamist, the run ends as a successful no-op.
- If a path to an existing VDS is generated (i.e. by checking `VDS/../_SUCCESS`) markers, the planCombine
  method routes the process to the branch 'done', meaning the Combiner process is never run, so the registration
  call inside the process never fires. Registering the VDS is to be done by hand. Deleting an unwanted VDS is
  also an irregular behaviour - Hail cannot overwrite a partially written VDS, and the Hail team have confirmed
  this is a feature not a bug - in the event of partial writes, deletion and registration is by hand.
