# A missing or ambiguous gVCF fails the whole run

Before any cluster starts, every requested SG must have exactly one distinct active gVCF. Several active
analyses pointing at the same file count as one. If any SG has none, or has two different files, the run
fails and every offender is listed in one error. A VDS is a single artefact covering the whole set: one
that silently leaves out SGs, or picks an arbitrary gVCF, is worse than one that is never built. Because
the VDS is built in one job, there's also no per-SG output to salvage.

Analysis status is not checked. Any registered gVCF is treated as complete, because `completed` is a
Metamist status that CPG doesn't use reliably. Active/archived SGs are checked as part of the query, and
inactive cohorts cause the workflow to fail.
