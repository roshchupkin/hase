# Public evidence for the bounded Python 3 repair

18 October 2026. Engineering checks passed for the declared regression scope;
scientific review is **not assessed**. This is a selected export from the retained
local pilot archive. [Export manifest](exports.json) gives original and delivered
SHA-256 hashes and archive identifiers. It is not the complete execution archive.
Original records were not modified. Archive paths in JSON identify provenance;
they are not public hyperlinks.

| Observation | Record |
| --- | --- |
| Environment, source base, original run-manifest hashes and exact large-output maxima | [Sanitized summary](summary.json) |
| Independent OLS: 18 beta/SE/t comparisons passed | [Comparison and exact rational expectations](positive-comparison.json), [synthetic observation](positive-observed.json), [Kit export](positive-kit-export.json) |
| Negative control: one coefficient +0.01, exactly one failure | [Comparison](negative-control-comparison.json), [altered synthetic observation](negative-control-observed.json), [failed Kit export](negative-control-kit-export.json) |
| All 20,000,000 saved ExampleStudy associations matched legacy identities and provisional thresholds | [Per-file comparison](example-comparison.json), [Kit export](example-kit-export.json) |

Comparison/observation files are byte-identical selected records containing only
synthetic values or aggregate comparison statistics/hashes. Kit exports are labelled
sanitized derivatives: command/environment details and unrelated fields are omitted;
original hashes are embedded. The summary is also a labelled derivative. Kit's
missing Git context (`source.status=unavailable`) and unassessed scientific review
remain explicit. Separate outer manifest hashes do not fill in Kit's missing fields.
No private Kit implementation or Foundation documents are included.

The archived Python 3.13 image was
`sha256:c006446d193bfa1a817daed5a4bcb651d90592337ee2c398ea3a6fe4ad659b08`.
Source base: `32b6e0a838a8691b41b7d4db9b13af0fe6cb0d5d`.
The original CLI help and full ExampleStudy passed in that image. The two failed
modern example attempts (legacy HDF byte metadata, then missing PyTables import)
remain archived and are identified in the summary; the final successful example
exercised the corrected reader. Full logs/manifests and large arrays remain outside
this public bundle. Large-study data provenance/redistribution still needs review.

Absolute 1e-12 / relative 1e-10 provisional thresholds were not widened. Largest
absolute legacy differences: SE 1.6209256159527285e-14, t 8.393286066166183e-14,
reconstructed beta 8.570921750106208e-14, MAF 0. The large comparison measures
consistency; the tiny exact-rational OLS reference is the independent check.

[Run the tiny reference](../RUN_REFERENCE.md) without Kit or the original workspace.
Reproducing the large archived comparison requires its original trusted run
manifests and arrays, which are intentionally not included here.

Remaining gates: named scientific reviewer, model/tolerance approval, operational
and release authority, Apache-2.0 LICENSE versus GPLv3 CLI-banner discrepancy,
and stewardship. Converter/meta-analysis runtime and other HASE modes are unverified.
This evidence grants no release, global Kit/Runner acceptance or deployment status.
