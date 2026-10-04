"""Current label scoring policy, separate from frozen historical benchmarks."""
import re
from .functional_benchmark import normalize, informative as historical_informative, PENDING


def identifier_only_product(product):
    """Locus/entry designations do not describe protein function."""
    text = re.sub(r'\s+', ' ', str(product or '').strip().lower())
    text = re.sub(r'^(?:(?:putative|probable|predicted|conserved)\s+)+', '', text)
    return bool(re.fullmatch(
        r'(?:(?:protein|gene product)\s+(?:[a-z]|(?:gp|orf)?\d+(?:\.\d+)?[a-z]?|'
        r'[a-z]+\d+[a-z0-9.-]*|[a-z0-9]+[-.][a-z0-9]+)|'
        r'(?:gp|orf)\s*\d+(?:\.\d+)?[a-z]?(?:\s+protein)?)', text))


def informative(product):
    text = normalize(product)
    return (historical_informative(product) and not identifier_only_product(product)
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
