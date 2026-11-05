# A VDS is named by its content and never overwritten

The output is `<hash>_<count>.vds`. The hash is taken over the sorted SG IDs plus the sorted refresh IDs,
and the count is the number of SGs. Refresh IDs are part of the hash because a refresh changes a VDS
without changing its SG set. Without them, a refresh would get the base VDS's name and become a no-op.

An existing VDS is never overwritten, even with `force`. If the target path already exists, the first
stage stops before any cluster starts:

- If the VDS is registered in Metamist, the run ends as a successful no-op.
- If it isn't registered, the run fails with "would create X, but X already exists and is not registered".
  It doesn't register the VDS itself. A registration made after the fact would lack the run's metadata
  (cohorts, base VDS, refresh list), so the record would be inaccurate. Deleting a VDS is abnormal and is
  left to a person.
