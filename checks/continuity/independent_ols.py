"""Exact rational Gram-Schmidt OLS reference; does not import HASE or NumPy."""
from fractions import Fraction
from math import sqrt


def dot(left, right):
    return sum((a * b for a, b in zip(left, right)), Fraction(0))


def reference(fixture):
    n = len(fixture["subject_ids"])
    covariates = fixture["covariates"]
    phenotypes = fixture["phenotypes"]
    if len(covariates) != n or len(phenotypes) != n:
        raise ValueError("Input subject counts differ")
    cov_count = len(fixture["covariate_ids"])
    phenotype_count = len(fixture["phenotype_ids"])
    if any(len(row) != cov_count for row in covariates):
        raise ValueError("Covariate dimensions differ")
    if any(len(row) != phenotype_count for row in phenotypes):
        raise ValueError("Phenotype dimensions differ")
    if len(fixture["genotypes"]) != len(fixture["variant_ids"]):
        raise ValueError("Variant dimensions differ")
    cells = []
    for genotype in fixture["genotypes"]:
        if len(genotype) != n:
            raise ValueError("Genotype subject count differs")
        columns = [[Fraction(1)] * n]
        columns += [[Fraction(row[j]) for row in covariates] for j in range(cov_count)]
        columns += [[Fraction(value) for value in genotype]]
        k = len(columns)
        df = n - k
        if df <= 0:
            raise ValueError("Non-positive residual degrees of freedom")
        # X = U R, with orthogonal unnormalised columns U and unit-diagonal R.
        # Fractions avoid rounding in projection, fitting and residual variance.
        orthogonal = []
        triangular = [[Fraction(0) for _ in range(k)] for _ in range(k)]
        for j, column in enumerate(columns):
            residual = list(column)
            for i, previous in enumerate(orthogonal):
                coefficient = dot(previous, column) / dot(previous, previous)
                triangular[i][j] = coefficient
                residual = [a - coefficient * b for a, b in zip(residual, previous)]
            if dot(residual, residual) == 0:
                raise ValueError("Rank-deficient design is outside this fixture's contract")
            orthogonal.append(residual)
            triangular[j][j] = Fraction(1)
        row = []
        for phenotype in range(phenotype_count):
            y = [Fraction(values[phenotype]) for values in phenotypes]
            gamma = [dot(u, y) / dot(u, u) for u in orthogonal]
            beta = [Fraction(0)] * k
            for j in reversed(range(k)):
                beta[j] = gamma[j] - sum(triangular[j][l] * beta[l] for l in range(j + 1, k))
            errors = [y[i] - sum(columns[j][i] * beta[j] for j in range(k)) for i in range(n)]
            sse = dot(errors, errors)
            # The last regressor's coefficient variance is sigma^2 / ||U_last||^2.
            variance = (sse / df) / dot(orthogonal[-1], orthogonal[-1])
            if variance <= 0:
                raise ValueError("Non-positive coefficient variance is outside this contract")
            se = sqrt(float(variance))
            row.append({"beta": float(beta[-1]), "standard_error": se,
                        "t_statistic": float(beta[-1]) / se,
                        "beta_exact": str(beta[-1]), "variance_exact": str(variance),
                        "residual_sum_squares_exact": str(sse)})
        cells.append(row)
    return {"degrees_of_freedom": df, "cells": cells}


def verify_hand_calculated_anchor():
    """A five-point example: x'x=10, SSE=14, df=3, beta=2, var(beta)=7/15."""
    anchor = {"subject_ids": list(range(5)), "covariate_ids": [],
              "covariates": [[] for _ in range(5)], "variant_ids": ["x"],
              "genotypes": [[-2, -1, 0, 1, 2]], "phenotype_ids": ["y"],
              "phenotypes": [[-3], [-4], [2], [0], [5]]}
    result = reference(anchor)
    cell = result["cells"][0][0]
    if (result["degrees_of_freedom"], cell["beta_exact"], cell["variance_exact"],
            cell["residual_sum_squares_exact"]) != (3, "2", "7/15", "14"):
        raise ValueError("Independent OLS implementation failed its hand-calculated anchor")
