"""Extract one tiny regression observation from original HASE; Python 2 compatible."""
from __future__ import print_function
import hashlib
import json
import os
import sys

import numpy as np

# The scientific source is mounted separately; never edit it to import the probe.
sys.path.insert(0, "/source")
from hdgwas import hdregression as hase


def main():
    fixture_path, output_path = sys.argv[1:]
    with open(fixture_path, "rb") as stream:
        raw = stream.read()
    fixture = json.loads(raw)
    cov = np.asarray(fixture["covariates"], dtype=np.float64)
    gen = np.asarray(fixture["genotypes"], dtype=np.float64)
    phen = np.asarray(fixture["phenotypes"], dtype=np.float64)
    acov = hase.A_covariates(cov, intercept=True)
    atest = hase.A_tests(cov, gen, intercept=True)
    inverse = hase.A_inverse(acov, atest)
    n_constant = inverse.shape[1] - 1
    df = phen.shape[0] - inverse.shape[1]
    t_stat, se = hase.HASE(hase.B4(phen, gen), inverse,
                          hase.B_covariates(cov, phen, intercept=True),
                          hase.C_matrix(phen), n_constant, df)
    beta = t_stat * se
    cells = [[{"beta": float(beta[v, 0, p]),
               "standard_error": float(se[v, 0, p]),
               "t_statistic": float(t_stat[v, 0, p])}
              for p in range(phen.shape[1])] for v in range(gen.shape[0])]
    with open(hase.__file__.replace(".pyc", ".py"), "rb") as stream:
        kernel_hash = hashlib.sha256(stream.read()).hexdigest()
    result = {"schema_version": 1, "fixture_sha256": hashlib.sha256(raw).hexdigest(),
              "kernel_sha256": kernel_hash, "degrees_of_freedom": int(df),
              "subject_ids": fixture["subject_ids"], "variant_ids": fixture["variant_ids"],
              "phenotype_ids": fixture["phenotype_ids"], "covariate_ids": fixture["covariate_ids"],
              "intercept": True, "cells": cells,
              "scope": "Regression matrix kernel only; no CLI, file alignment, conversion or meta-analysis claim."}
    with open(output_path, "w") as stream:
        json.dump(result, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write("\n")
    print("Saved HASE observations to " + output_path)


if __name__ == "__main__":
    main()
