"""Current label scoring policy, separate from frozen historical benchmarks."""
import re
from .functional_benchmark import normalize, informative as historical_informative, PENDING


def informative(product):
    text = normalize(product)
    return (historical_informative(product)
            and text not in {'na', 'n a', 'nan', 'none', 'null'}
            and not re.search(r'\b(duf\d*|upf\d*)\b', text))


def classify(prediction, reference, synonyms):
    if not informative(reference):
        return 'NON_EVALUABLE_REFERENCE'
    if not informative(prediction):
        return 'ABSTENTION'
    if normalize(prediction) == normalize(reference):
        return 'EXACT_PRODUCT_AGREEMENT'
    if (normalize(prediction), normalize(reference)) in synonyms:
        return 'EQUIVALENT_FUNCTION'
    return PENDING
