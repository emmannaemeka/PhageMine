# Frozen functional-utility definition

Frozen on 18 September 2026 at repository commit `4933d8fa532a0157c24c3c01dbf1d8adfa787ed7`, before any genuine adjudication outcome was imported. For genome `g`, `U_g = sum(u_i) / N_g`, where `N_g` is the number of frozen evaluable reference loci for that genome. Weights are `CORRECT=1.00`, `PARTIALLY_CORRECT=0.50`, `TOO_GENERAL=0.25`, `UNSUPPORTED_SPECIFIC=0.00`, and `INCORRECT=0.00`. `UNRESOLVABLE` is missing and is reported with counts and identification bounds; it is never silently converted to zero. `NOT_EVALUABLE` is excluded from `N_g`.

The strict secondary endpoint is `CORRECT / N_g`. The inclusive sensitivity endpoint is `(CORRECT + PARTIALLY_CORRECT) / N_g`. Unsupported specificity remains a separate safety endpoint. The primary paired effect is the seven-genome mean of `D_g = U_g(PhageMine)-U_g(Pharokka)`, with frozen margin `+0.05`; a lower 95% confidence bound must exceed +0.05 to meet the prespecified criterion. Exact sign permutation and genome-level bootstrap are reported as small-panel sensitivities. Loci are never treated as independent biological replicates for the primary comparison.

This definition is immutable for v1.2 adjudication. It does not alter annotation behavior, frozen inputs, matching, synonym rules, or the pre-adjudication benchmark.
