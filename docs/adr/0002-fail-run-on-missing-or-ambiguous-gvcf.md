# A missing or ambiguous gVCF fails the whole run

Before any cluster starts, every requested SG must have exactly one distinct active gVCF. Several active
analyses pointing at the same file count as one. If any SG has none, or has two different files, the run
fails and every offender is listed in one error. A VDS is a single artefact covering the whole set: one
that silently leaves out SGs, or picks an arbitrary gVCF, is worse than one that is never built. Because
the VDS is built in one job, there's also no per-SG output to salvage. This is the opposite of nf-align's
ADR 0007, which drops failed SGs and carries on.

Analysis status is not checked. Any registered gVCF is treated as complete, because `completed` is a
Metamist status that CPG doesn't use reliably.
