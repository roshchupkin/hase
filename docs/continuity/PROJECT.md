# HASE reference pilot

HASE computes high-dimensional association regressions. This branch adds a small regression-kernel reference and a bounded Python 3.13
compatibility repair. Regression formulas and singular-matrix policy are preserved.

Source baseline: `32b6e0a838a8691b41b7d4db9b13af0fe6cb0d5d`.
Branch: `continuity/hase-python3`; draft candidate, not released.
Pilot initiator: Gennady Roshchupkin. Operational maintainer, scientific reviewer
and release authority are not yet confirmed for the pilot.

The original ExampleStudy and CLI ran in the recovered Python 2.7 Docker image.
The new probe runs the regression matrix functions on 16 synthetic subjects,
three variants, two phenotypes, an intercept and two covariates (12 residual df).
The public [Docker reference command](RUN_REFERENCE.md) runs the probe and independent
comparison without private Kit access or host scientific dependencies. The archived
pilot also recorded checks through Kit; see [public evidence](evidence/README.md).

The upstream README/LICENSE indicate Apache-2.0; the CLI banner says GPLv3. This
existing discrepancy remains unresolved. Publication citation: Roshchupkin et al.,
Scientific Reports 6, 36076 (2016), doi:10.1038/srep36076.

Read [VALIDATION.md](VALIDATION.md) for the proposed numerical contract. The reviewer remains
unassigned; a passed execution/comparison is not scientific approval. Other modes,
identifier alignment, chunk invariance and degenerate inputs remain untested.
