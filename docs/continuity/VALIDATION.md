# HASE small regression reference contract

Status: provisional engineering reference. Scientific reviewer: unassigned.
Scientific approval: pending. Tolerances are proposals, not approved thresholds.

## Model and input

`checks/continuity/fixture.json` contains 16 explicit synthetic integer observations.
Rows have declared subject order; genotype rows are variants, phenotype columns
are traits. Each variant/phenotype pair fits ordinary least squares with intercept
and both declared covariates. The four-column design has full rank; residual df=12.
All arithmetic is deterministic. No random seed or participant data is involved.

`hase_probe.py` calls the original A_covariates, A_tests, A_inverse, B_covariates,
B4, C_matrix and HASE functions. HASE returns t statistics and standard errors;
the probe reconstructs the tested coefficient as t times SE. Actual output is JSON
with input identity, dimensions, kernel hash and all six result cells.

## Independent expectation

`independent_ols.py` does not import HASE, NumPy or SciPy. It uses fractions for
Gram-Schmidt orthogonalisation, a triangular coefficient solution and residual
sum of squares. For the final genotype coefficient, variance is residual variance
divided by the squared norm of its orthogonalised regressor. Fractions remain
exact until the final square root and division for SE/t.

The independent calculation first checks a separate hand-calculable five-point
anchor: x=[-2,-1,0,1,2], y=[-3,-4,2,0,5], beta=2, SSE=14, df=3,
coefficient variance=7/15. This catches a basic reference-implementation error.

Compare all 18 beta/SE/t values. Provisional absolute tolerance is 1e-12 and
relative tolerance 1e-10: a rounding allowance for this small full-rank integer
fixture against an exact calculation. These thresholds do not apply to arbitrary
GWAS, different precision or ill-conditioned data without review. Exact rational
expectations and observed errors are exported for the reviewer's inspection.

## Reproduction and archived Kit integration

For a fresh checkout use the [standalone Docker command](RUN_REFERENCE.md). It
produces observed.json and comparison.json directly using the same comparator;
it does not execute Kit or produce a Kit report. The following describes the
archived local integration, with [selected public exports](evidence/README.md).

The local Docker controller runs the legacy probe, exports its observation and
then creates a fresh Python 3 container. A temporary copy of this project receives
the observed JSON at `.continuity/observations/observed.json`. The Kit executes
the `hase-ols-reference` check and saves its normal JSON evidence. The comparator
also exports detailed numerical comparisons beneath `.continuity/comparisons/`.

Declared reference hashes cover fixture, contract, independent implementation,
comparator, probe, HASE kernel and observed input. Container Git context remains
unknown in the Kit report; the controller separately records the host worktree's
base SHA, dirty state, file hashes, patch hash, image IDs and execution outcomes.

A negative-control run in a separate temporary project adds 0.01 to one observed
coefficient. The same check must fail without changing references or tolerances.
This is labelled negative-control evidence, not a failed HASE scientific run.

## Limits

This exercises the matrix kernel, not the file-reading CLI, subject matching,
MAF filtering, chunking, p-values, meta-analysis, singular designs, constant
phenotypes, missing data or scientific validity of the original method. The [bounded Python 3 port](PYTHON3.md) has scoped engineering evidence;
actual reviewer acceptance remains separate work. No scientific
implementation is repaired or numerically redefined in this increment.
