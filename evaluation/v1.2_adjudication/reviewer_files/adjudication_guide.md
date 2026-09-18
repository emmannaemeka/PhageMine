# Frozen v1.2 functional adjudication guide

Review only the unresolved cases in `blinded_adjudication_cases.xlsx`. The seven genomes and all prediction files are evaluation data. Do not use these cases to tune PhageMine v1.2, alter the frozen benchmark, or infer which tool produced Prediction_A or Prediction_B. Evidence fields are copied from the frozen archive and may contain conflicting or incomplete database support. No new searches are part of this review.

For each target prediction, complete `final_judgment` (also recorded as `adjudication_class`), evidence basis, uncertainty, evidence support, strength, specificity, conflict, preferred product name, error category, and a short rationale (`reviewer_notes`). `NOT_EVALUABLE` is reserved for cases excluded by the frozen rules; `UNRESOLVABLE` means the case is eligible but the supplied evidence cannot distinguish plausible functions. A blank decision is missing data, never an incorrect call.

Use the controlled definitions below:

* `CORRECT`: same biological function or accepted equivalent.
* `PARTIALLY_CORRECT`: meaningful functional information, but incomplete or partly mismatched.
* `TOO_GENERAL`: compatible but materially less informative than the supported reference.
* `UNSUPPORTED_SPECIFIC`: a specific claim exceeds the supplied evidence.
* `INCORRECT`: conflicts with the supported function.
* `UNRESOLVABLE`: evidence is insufficient to distinguish plausible interpretations.
* `NOT_EVALUABLE`: cannot contribute under the frozen reference-evaluability rules.

`CORRECT` is the primary strict positive outcome. `PARTIALLY_CORRECT` is retained separately and receives half credit only in the prespecified utility sensitivity analysis. `TOO_GENERAL`, unsupported specificity, and incorrect calls do not receive strict credit. Do not treat a database hit as independent truth; record conflict and source limitations.

Return the completed workbook or CSV as a new raw submission with reviewer ID and timestamp. Never overwrite the blank template, the frozen archive, or the blinding key. Multiple reviewers must submit independently; disagreements remain separate until a documented consensus step.
