# Python 3 compatibility verification

Date: 18 October 2026. Engineering scope: passed. Scientific approval: pending.

Original base SHA: 32b6e0a838a8691b41b7d4db9b13af0fe6cb0d5d.
Candidate branch: continuity/hase-python3. No release or scientific approval.
Modern image: sha256:c006446d193bfa1a817daed5a4bcb651d90592337ee2c398ea3a6fe4ad659b08.

Completed in fresh, resource-limited, network-disabled Docker containers:

- CLI help passed.
- All 18 small-fixture beta/SE/t comparisons matched independent exact OLS.
- Deliberately adding 0.01 to one coefficient produced a failed Kit report.
- The bundled ExampleStudy completed and wrote all four result files.
- Kit compared all 20,000,000 saved associations against the archived Python 2
  baseline. File coverage, shapes, association identities, input hashes and
  finite MAF/SE/t/reconstructed beta values matched within unchanged provisional
  abs 1e-12 / rel 1e-10 tolerances. No identity or numerical failures occurred.

Maximum absolute differences against the large legacy outputs: MAF=0,
SE=1.63e-14, t=8.40e-14, reconstructed coefficient=8.58e-14 (rounded upward).
Large-output consistency is distinct from independent validation of the study.

Evidence run IDs (UTC):

- CLI: 20261018T120208Z-help-c473bda2.
- Reference producer: 20261018T120210Z-python3-reference-produce-77fe77fe.
- Kit reference/negative control: 20261018T120212Z-python3-reference-verify-06338b1e.
- Final ExampleStudy: 20261018T120252Z-example-b6fb52ac.
- Kit legacy comparison: 20261018T120346Z-python3-example-comparison-96edd05a.

The final ExampleStudy also exercises the narrow legacy-HDF5 ID reader and its
imports. After fixing a missing PyTables import, only the affected example and
new output comparison were executed; the kernel, reference checker, fixture,
contract and locked environment were unchanged from their passing reference run.
Two earlier failed examples (byte metadata, then missing import) and their logs
are retained in the local pilot evidence. No passing failure state was invented.

The tested final runtime patch SHA-256 was
ed54a16ea844370aa6effd7d0a50096675de4388269e9065b21cd94abf9faf49.
This verification file and the status paragraph were added afterwards; all runtime
and checking code hashes still identify the executed files. Outer Docker records
retain source/patch/input/output hashes, image identities, logs and cleanup.
Kit reports retain scientific_review=not_assessed and their missing container Git
context; outer source provenance supplies separate observations without rewriting
those Kit reports.

The original and reference checkouts remain preserved. No host scientific package
installation, numerical policy change, complete-mode certification, SPACE adapter,
Platform deployment or scientific release was performed. Review of the provisional
model/tolerances, licence discrepancy and release ownership remains outstanding.

Selected [public evidence records](evidence/README.md) resolve within this branch.
The [standalone Docker guide](RUN_REFERENCE.md) reproduces the tiny independent
reference from a fresh checkout without the private integration workspace.
