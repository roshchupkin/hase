"""Compare only hash-identified outputs from this pilot's trusted local Docker runs."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def result_manifest(path):
    run = json.loads(path.read_text())
    if run.get("status") != "passed" or run.get("exit_code") != 0:
        raise ValueError("An input ExampleStudy execution did not pass")
    artifacts = {Path(a["path"]).name: a for a in run["artifacts"]
                 if a["path"].startswith("results/") and a["path"].endswith(".npy")}
    if not artifacts:
        raise ValueError("No saved result artifacts")
    return run, artifacts


def compare(legacy, modern, legacy_manifest, modern_manifest):
    contract = json.loads((HERE/"contract.json").read_text())
    atol, rtol = contract["absolute_tolerance"], contract["relative_tolerance"]
    old_run, old_files = result_manifest(legacy_manifest)
    new_run, new_files = result_manifest(modern_manifest)
    if set(old_files) != set(new_files):
        raise ValueError("Result filenames/coverage differ")
    old_inputs, new_inputs = old_run["source_files_sha256"], new_run["source_files_sha256"]
    inputs = [p for p in old_inputs if p.startswith("test/ExampleStudy/")]
    if not inputs or any(old_inputs[p] != new_inputs.get(p) for p in inputs):
        raise ValueError("ExampleStudy input file hashes differ")
    files, failures, total = [], [], 0
    for name in sorted(old_files):
        paths = (legacy/name, modern/name)
        hashes = [digest(p) for p in paths]
        if hashes != [old_files[name]["sha256"], new_files[name]["sha256"]]:
            raise ValueError("Result artifact hash differs from its recorded run: " + name)
        # These dictionaries were created by our own recorded jobs, are mounted
        # read-only, and were hashed above. This is NOT an arbitrary pickle loader.
        old, new = [np.load(p, allow_pickle=True, encoding="latin1").item() for p in paths]
        expected_keys = {"index", "phenotype", "MAF", "SE", "t-stat"}
        if set(old) != expected_keys or set(new) != expected_keys:
            raise ValueError("Result dictionary schema differs: " + name)
        shape = np.asarray(old["t-stat"]).shape
        if len(shape) != 1 or any(np.asarray(d[k]).shape != shape for d in (old, new) for k in expected_keys):
            raise ValueError("Result dimensions differ: " + name)
        total += shape[0]
        file = {"name": name, "associations": shape[0], "sha256": hashes,
                "identity_mismatches": 0, "metrics": {}}
        for metric in ("MAF", "SE", "t-stat", "beta"):
            stat = {"compared": 0, "nonfinite": 0, "outside_tolerance": 0,
                    "max_absolute_error": 0.0, "max_error_to_allowance": 0.0}
            for start in range(0, shape[0], 100000):
                selection = slice(start, min(start+100000, shape[0]))
                if metric == "MAF":
                    for key in ("index", "phenotype"):
                        left = np.asarray(old[key][selection]).astype(str)
                        right = np.asarray(new[key][selection]).astype(str)
                        file["identity_mismatches"] += int(np.count_nonzero(left != right))
                left = old["t-stat"][selection] * old["SE"][selection] if metric == "beta" else old[metric][selection]
                right = new["t-stat"][selection] * new["SE"][selection] if metric == "beta" else new[metric][selection]
                finite = np.isfinite(left) & np.isfinite(right)
                difference = np.abs(left-right)
                allowance = atol + rtol*np.abs(left)
                stat["compared"] += int(left.size)
                stat["nonfinite"] += int(np.count_nonzero(~finite))
                stat["outside_tolerance"] += int(np.count_nonzero(finite & (difference > allowance)))
                if np.any(finite):
                    stat["max_absolute_error"] = max(stat["max_absolute_error"], float(np.max(difference[finite])))
                    stat["max_error_to_allowance"] = max(stat["max_error_to_allowance"], float(np.max((difference/allowance)[finite])))
            if stat["nonfinite"] or stat["outside_tolerance"]:
                failures.append(name + "/" + metric + " differs")
            file["metrics"][metric] = stat
        if file["identity_mismatches"]:
            failures.append(name + "/association identities differ")
        files.append(file)
    if total != 20000000:
        failures.append("Expected 20,000,000 saved associations for this bundled fixture")
    return {"schema_version": 1, "status": "failed" if failures else "passed",
            "scientific_review": "not_assessed", "associations": total,
            "absolute_tolerance": atol, "relative_tolerance": rtol,
            "legacy_manifest_sha256": digest(legacy_manifest), "modern_manifest_sha256": digest(modern_manifest),
            "input_files_matched": len(inputs), "files": files, "failures": failures,
            "scope": "Full saved ExampleStudy output consistency with the archived legacy run, using unchanged provisional numerical tolerances. Not an independent validation of the large study."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy", type=Path, required=True)
    parser.add_argument("--modern", type=Path, required=True)
    parser.add_argument("--legacy-manifest", type=Path, required=True)
    parser.add_argument("--modern-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = compare(args.legacy, args.modern, args.legacy_manifest, args.modern_manifest)
    except (ValueError, KeyError, TypeError, OSError) as error:
        result = {"status": "failed", "scientific_review": "not_assessed", "error": str(error)}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"status": result["status"], "output": str(args.output)}))
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
