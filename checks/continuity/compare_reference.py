"""Kit entry point: compare an observed HASE kernel result with independent OLS."""
import argparse
import hashlib
import json
import math
from pathlib import Path

from independent_ols import reference, verify_hand_calculated_anchor

HERE = Path(__file__).resolve().parent


def compare(observed_path):
    verify_hand_calculated_anchor()
    raw = (HERE / "fixture.json").read_bytes()
    fixture = json.loads(raw)
    contract = json.loads((HERE / "contract.json").read_text())
    observed = json.loads(observed_path.read_text())
    expected = reference(fixture)
    failures = []
    if observed.get("fixture_sha256") != hashlib.sha256(raw).hexdigest():
        failures.append("Observed fixture hash differs")
    kernel = HERE.parents[1] / "hdgwas/hdregression.py"
    if observed.get("kernel_sha256") != hashlib.sha256(kernel.read_bytes()).hexdigest():
        failures.append("Observed HASE kernel hash differs")
    for key in ("subject_ids", "variant_ids", "phenotype_ids", "covariate_ids"):
        if observed.get(key) != fixture[key]:
            failures.append(key + " differ")
    if observed.get("intercept") is not True:
        failures.append("Intercept setting differs")
    if observed.get("degrees_of_freedom") != expected["degrees_of_freedom"]:
        failures.append("Residual degrees of freedom differ")
    cells = observed.get("cells")
    rows = len(fixture["variant_ids"])
    columns = len(fixture["phenotype_ids"])
    if not isinstance(cells, list) or len(cells) != rows or any(
            not isinstance(row, list) or len(row) != columns for row in cells):
        raise ValueError("Observed result dimensions differ")
    comparisons = []
    for v, variant in enumerate(fixture["variant_ids"]):
        for p, phenotype in enumerate(fixture["phenotype_ids"]):
            for metric in ("beta", "standard_error", "t_statistic"):
                value = cells[v][p].get(metric)
                target = expected["cells"][v][p][metric]
                passed = (type(value) in (int, float) and math.isfinite(value) and
                          math.isclose(value, target, rel_tol=contract["relative_tolerance"],
                                       abs_tol=contract["absolute_tolerance"]))
                comparisons.append({"variant": variant, "phenotype": phenotype,
                                    "metric": metric, "actual": value, "expected": target,
                                    "absolute_error": abs(value - target) if type(value) in (int, float) else None,
                                    "passed": passed})
                if not passed:
                    failures.append(variant + "/" + phenotype + "/" + metric + " differs")
    return {"schema_version": 1, "status": "failed" if failures else "passed",
            "scientific_review": "not_assessed", "contract_status": contract["status"],
            "fixture_sha256": hashlib.sha256(raw).hexdigest(),
            "observed_sha256": hashlib.sha256(observed_path.read_bytes()).hexdigest(),
            "expected_exact": expected, "comparisons": comparisons, "failures": failures,
            "scope": "Three full-rank variants, two phenotypes, intercept and two covariates; provisional tolerances."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--observed", type=Path, default=Path(".continuity/observations/observed.json"))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        result = compare(args.observed)
    except (ValueError, KeyError, TypeError, OSError) as error:
        result = {"status": "failed", "scientific_review": "not_assessed", "error": str(error)}
    text = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        # Preserve earlier comparison evidence, just as Kit preserves reports.
        with args.output.open("x") as stream:
            stream.write(text)
    print(text, end="")
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
