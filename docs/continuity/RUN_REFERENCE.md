# Run the tiny independent reference

From a fresh checkout of this PR branch, on Linux with Docker and a POSIX shell:

```sh
sh checks/continuity/run_reference.sh
```

Run from the checkout root. No Python, scientific packages, private Kit repository,
archived pilot data or workspace controller is required on the host. Building
requires network access to the public Python base image and package index. The
recipe pins an immutable base and resolved dependency versions. This pilot was
checked on Linux x86_64; other architectures have not been verified. The image ID
identifies the local build; rebuilding pinned versions need not yield identical
image bytes. Downloads are version pinned, not package-hash locked.

The script builds `continuity-hase-python3:pilot-01`, then uses a disposable,
network-disabled container with read-only source/root and bounded resources. It
runs the synthetic HASE probe and exact-rational OLS comparator. Exit code 0 and
`"status": "passed"` in `continuity-reference-output/comparison.json` mean all 18
beta/SE/t comparisons passed. `observed.json` and `image-id.txt` are saved alongside
it. This direct command does not run Kit; archived Kit evidence is [here](evidence/README.md).

The output directory must be new. To retain another run or choose a location:

```sh
sh checks/continuity/run_reference.sh /path/to/new-output-directory
```

Use an existing parent directory and a path without commas (Docker bind-mount
syntax). On failure, retain any partial output and the command log. A failed
build or comparison exits nonzero; inspect the output before trying a new directory.

The [contract](VALIDATION.md) remains provisional: 16 synthetic subjects, three
variants, two phenotypes, intercept and two covariates; absolute tolerance 1e-12,
relative tolerance 1e-10. A pass is an engineering observation. Model/tolerance
approval needs a named scientific reviewer. This command does not run ExampleStudy,
converters or meta-analysis, or validate degenerate designs or the scientific method.
