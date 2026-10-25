"""Review a phage genome package against NCBI GenBank requirements before submission.

Each rule mirrors a published NCBI or INSDC requirement and cites it. This is a
pre-submission review only: NCBI's own validator (table2asn) and GenBank staff
remain the authority. Findings are ERROR (expected to be flagged or returned by
NCBI) or WARNING (review recommended).

Sources (consulted October 2026):
  DISC    https://www.ncbi.nlm.nih.gov/genbank/asndisc.examples/
  VAL     https://www.ncbi.nlm.nih.gov/genbank/genome_validation/
  NOMEN   https://www.ncbi.nlm.nih.gov/genbank/internatprot_nomenguide/
  ANNOT   https://www.ncbi.nlm.nih.gov/genbank/genomesubmit_annotation/
  SOURCE  https://www.ncbi.nlm.nih.gov/WebSub/html/help/genbank-source-table.html
  COUNTRY https://www.ncbi.nlm.nih.gov/genbank/collab/country/
  MISSING https://www.insdc.org/technical-specifications/missing-value-reporting/
  SEQ     https://www.ncbi.nlm.nih.gov/genbank/sequencecheck/virus/
"""
from __future__ import annotations

import re
from datetime import date
from pathlib import Path

SOURCES = {
    "DISC": "https://www.ncbi.nlm.nih.gov/genbank/asndisc.examples/",
    "VAL": "https://www.ncbi.nlm.nih.gov/genbank/genome_validation/",
    "NOMEN": "https://www.ncbi.nlm.nih.gov/genbank/internatprot_nomenguide/",
    "ANNOT": "https://www.ncbi.nlm.nih.gov/genbank/genomesubmit_annotation/",
    "SOURCE": "https://www.ncbi.nlm.nih.gov/WebSub/html/help/genbank-source-table.html",
    "COUNTRY": "https://www.ncbi.nlm.nih.gov/genbank/collab/country/",
    "MISSING": "https://www.insdc.org/technical-specifications/missing-value-reporting/",
    "SEQ": "https://www.ncbi.nlm.nih.gov/genbank/sequencecheck/virus/",
}

# INSDC geo_loc_name list (COUNTRY), retrieved October 2026.
APPROVED_GEO_NAMES = frozenset("""Afghanistan|Albania|Algeria|American Samoa|Andorra|Angola|Anguilla|Antarctica|Antigua and Barbuda|Arctic Ocean|Argentina|Armenia|Aruba|Ashmore and Cartier Islands|Atlantic Ocean|Australia|Austria|Azerbaijan|Bahamas|Bahrain|Baltic Sea|Baker Island|Bangladesh|Barbados|Bassas da India|Belarus|Belgium|Belize|Benin|Bermuda|Bhutan|Bolivia|Borneo|Bosnia and Herzegovina|Botswana|Bouvet Island|Brazil|British Virgin Islands|Brunei|Bulgaria|Burkina Faso|Burundi|Cambodia|Cameroon|Canada|Cape Verde|Cayman Islands|Central African Republic|Chad|Chile|China|Christmas Island|Clipperton Island|Cocos Islands|Colombia|Comoros|Cook Islands|Coral Sea Islands|Costa Rica|Cote d'Ivoire|Croatia|Cuba|Curacao|Cyprus|Czechia|Democratic Republic of the Congo|Denmark|Djibouti|Dominica|Dominican Republic|Ecuador|Egypt|El Salvador|Equatorial Guinea|Eritrea|Estonia|Eswatini|Ethiopia|Europa Island|Falkland Islands (Islas Malvinas)|Faroe Islands|Fiji|Finland|France|French Guiana|French Polynesia|French Southern and Antarctic Lands|Gabon|Gambia|Gaza Strip|Georgia|Germany|Ghana|Gibraltar|Glorioso Islands|Greece|Greenland|Grenada|Guadeloupe|Guam|Guatemala|Guernsey|Guinea|Guinea-Bissau|Guyana|Haiti|Heard Island and McDonald Islands|Honduras|Hong Kong|Howland Island|Hungary|Iceland|India|Indian Ocean|Indonesia|Iran|Iraq|Ireland|Isle of Man|Israel|Italy|Jamaica|Jan Mayen|Japan|Jarvis Island|Jersey|Johnston Atoll|Jordan|Juan de Nova Island|Kazakhstan|Kenya|Kerguelen Archipelago|Kingman Reef|Kiribati|Kosovo|Kuwait|Kyrgyzstan|Laos|Latvia|Lebanon|Lesotho|Liberia|Libya|Liechtenstein|Line Islands|Lithuania|Luxembourg|Macau|Madagascar|Malawi|Malaysia|Maldives|Mali|Malta|Marshall Islands|Martinique|Mauritania|Mauritius|Mayotte|Mediterranean Sea|Mexico|Micronesia, Federated States of|Midway Islands|Moldova|Monaco|Mongolia|Montenegro|Montserrat|Morocco|Mozambique|Myanmar|Namibia|Nauru|Navassa Island|Nepal|Netherlands|New Caledonia|New Zealand|Nicaragua|Niger|Nigeria|Niue|Norfolk Island|North Korea|North Macedonia|North Sea|Northern Mariana Islands|Norway|Oman|Pacific Ocean|Pakistan|Palau|Palmyra Atoll|Panama|Papua New Guinea|Paracel Islands|Paraguay|Peru|Philippines|Pitcairn Islands|Poland|Portugal|Puerto Rico|Qatar|Republic of the Congo|Reunion|Romania|Ross Sea|Russia|Rwanda|Saint Barthelemy|Saint Helena|Saint Kitts and Nevis|Saint Lucia|Saint Martin|Saint Pierre and Miquelon|Saint Vincent and the Grenadines|Samoa|San Marino|Sao Tome and Principe|Saudi Arabia|Senegal|Serbia|Seychelles|Sierra Leone|Singapore|Sint Maarten|Slovakia|Slovenia|Solomon Islands|Somalia|South Africa|South Georgia and the South Sandwich Islands|South Korea|South Sudan|Southern Ocean|Spain|Spratly Islands|Sri Lanka|State of Palestine|Sudan|Suriname|Svalbard|Sweden|Switzerland|Syria|Taiwan|Tajikistan|Tanzania|Tasman Sea|Thailand|Timor-Leste|Togo|Tokelau|Tonga|Trinidad and Tobago|Tromelin Island|Tunisia|Turkey|Turkmenistan|Turks and Caicos Islands|Tuvalu|Uganda|Ukraine|United Arab Emirates|United Kingdom|Uruguay|USA|Uzbekistan|Vanuatu|Venezuela|Viet Nam|Virgin Islands|Wake Island|Wallis and Futuna|West Bank|Western Sahara|Yemen|Zambia|Zimbabwe""".split("|"))

# INSDC missing-value terms (MISSING); mandatory fields use the "missing: <term>" form.
_MISSING_TERMS = ("not applicable", "not collected", "not provided", "restricted access", "control sample",
                  "sample group", "synthetic construct", "lab stock", "third party data",
                  "data agreement established pre-2023", "endangered species", "human-identifiable")
MISSING_VALUES = frozenset({"missing", *_MISSING_TERMS, *(f"missing: {term}" for term in _MISSING_TERMS)})

# Translation table 11 (bacteria and their phages): permitted initiation codons.
TABLE11_STARTS = frozenset({"ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA"})
MONTHS = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
IUPAC_DNA = set("ACGTRYSWKMBDHVN")


def finding(code: str, severity: str, message: str, action: str, source: str, ncbi_name: str | None = None, **where) -> dict:
    return {"code": code, "severity": severity, "message": message, "recommended_action": action,
            "ncbi_check": ncbi_name or code, "ncbi_reference": SOURCES[source], **where}


# ------------------------------------------------------------------ dates ----

def _parse_date(text: str) -> tuple[int, int | None, int | None] | None:
    """Accepts DD-Mon-YYYY, Mon-YYYY, YYYY, YYYY-MM, YYYY-MM-DD (SOURCE, VAL)."""
    text = text.strip()
    if m := re.fullmatch(r"(\d{2})-([A-Z][a-z]{2})-(\d{4})", text):
        if m.group(2) in MONTHS:
            return int(m.group(3)), MONTHS.index(m.group(2)) + 1, int(m.group(1))
    elif m := re.fullmatch(r"([A-Z][a-z]{2})-(\d{4})", text):
        if m.group(1) in MONTHS:
            return int(m.group(2)), MONTHS.index(m.group(1)) + 1, None
    elif m := re.fullmatch(r"(\d{4})(?:-(\d{2})(?:-(\d{2}))?)?", text):
        return int(m.group(1)), int(m.group(2)) if m.group(2) else None, int(m.group(3)) if m.group(3) else None
    return None


def check_collection_date(value: str | None, today: date | None = None) -> list[dict]:
    if not value or value.strip().lower() in MISSING_VALUES:
        return []
    today = today or date.today()
    parsed = _parse_date(value)
    if parsed is None:
        return [finding("BAD_COLLECTION_DATE", "ERROR", f"collection_date '{value}' is not in an accepted format.",
                        "Use DD-Mon-YYYY (14-Mar-2025), Mon-YYYY, YYYY, or ISO 8601 (2025-03-14), or an INSDC missing-value term.",
                        "VAL", "SEQ_DESCR_BadCollectionDate")]
    year, month, day = parsed
    try:
        when = date(year, month or 1, day or 1)
    except ValueError:
        return [finding("BAD_COLLECTION_DATE", "ERROR", f"collection_date '{value}' is not a real date.",
                        "Correct the day or month.", "VAL", "SEQ_DESCR_BadCollectionDate")]
    if when > today:
        return [finding("BAD_COLLECTION_DATE", "ERROR", f"collection_date '{value}' is in the future.",
                        "Use the date the sample was collected.", "VAL", "SEQ_DESCR_BadCollectionDate")]
    return []


def check_geo_loc_name(value: str | None) -> list[dict]:
    if not value or value.strip().lower() in MISSING_VALUES:
        return []
    country = value.split(":", 1)[0].strip()
    if country not in APPROVED_GEO_NAMES:
        close = sorted(n for n in APPROVED_GEO_NAMES if n.lower() == country.lower())
        hint = f" Did you mean '{close[0]}'?" if close else ""
        return [finding("BAD_GEO_LOC_NAME", "ERROR", f"'{country}' is not on the INSDC geo_loc_name list.{hint}",
                        "Use an approved country, ocean or sea name, then a colon and the region, e.g. 'Nigeria: Plateau State, Jos'.",
                        "COUNTRY", "SEQ_DESCR_BadGeoLocNameCode")]
    if ":" not in value and "," in value:
        return [finding("BAD_GEO_LOC_NAME", "ERROR", f"'{value}' separates the region with a comma.",
                        "Write the country, a colon, then the region: 'Country: region, locality'.", "SOURCE", "SEQ_DESCR_BadGeoLocNameCode")]
    return []


# ------------------------------------------------------------- metadata ----

def check_metadata(metadata) -> list[dict]:
    found = []
    found += check_collection_date(metadata.collection_date)
    found += check_geo_loc_name(metadata.geographic_location)
    organism = (metadata.organism or "").strip()
    if organism and not re.fullmatch(r"[A-Z][a-z]+ phage \S.*", organism):
        found.append(finding("PHAGE_ORGANISM_NAME", "WARNING",
                             f"Organism '{organism}' does not follow the '<host genus> phage <name>' form used in NCBI Taxonomy "
                             "(for example 'Escherichia phage T7').",
                             "Name the phage after its host genus, e.g. 'Xanthomonas phage vB_XcaM_JOS1'. NCBI may adjust names during taxonomy review.",
                             "ANNOT", "TAX_LOOKUP_MISMATCH"))
    for key in ("organism", "host", "isolation_source", "collection_date", "geographic_location", "isolate", "strain"):
        value = getattr(metadata, key, None)
        if value and re.search(r"[\[\]]", str(value)):
            found.append(finding("BRACKET_IN_SOURCE_VALUE", "ERROR", f"{key} contains square brackets, which break FASTA source modifiers.",
                                 "Remove [ and ] from the value.", "SOURCE", "SEQ_DESCR"))
    if metadata.submitter_email and not re.fullmatch(r"[^@\s]+@[^@\s]+\.[^@\s]+", metadata.submitter_email):
        found.append(finding("BAD_SUBMITTER_EMAIL", "ERROR", f"Submitter email '{metadata.submitter_email}' is not a valid address.",
                             "Correct the email address.", "VAL", "GENERIC.MissingPubInfo"))
    prefix = getattr(metadata, "locus_tag_prefix", None)
    if prefix:
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9]{2,11}", prefix):
            found.append(finding("BAD_LOCUS_TAG_PREFIX", "ERROR", f"Locus tag prefix '{prefix}' must be 3-12 letters or digits and start with a letter.",
                                 "Use the prefix assigned when you registered your BioProject.", "ANNOT", "BAD_LOCUS_TAG_FORMAT"))
    else:
        found.append(finding("NO_LOCUS_TAG_PREFIX", "WARNING", "No registered locus tag prefix: genes are written without locus_tags.",
                             "Register a BioProject and add its locus tag prefix to the metadata (locus_tag_prefix). "
                             "Genome submissions are expected to have a locus_tag on every gene.", "ANNOT", "MISSING_LOCUS_TAGS"))
    template = getattr(metadata, "submission_template", None)
    if template and not Path(template).expanduser().is_file():
        found.append(finding("SUBMISSION_TEMPLATE_NOT_FOUND", "ERROR", f"Submission template '{template}' was not found.",
                             "Give the path to the .sbt file downloaded from NCBI's Submission Template form.", "VAL", "GENERIC.MissingPubInfo"))
    elif not template:
        found.append(finding("LOCAL_SUBMISSION_TEMPLATE", "WARNING",
                         "submission.sbt was generated by PhageMine with names only; it lacks affiliation and contact details.",
                         "Create the official template with NCBI's Submission Template form (submit.ncbi.nlm.nih.gov/genbank/template/submission/) "
                         "and set submission_template in the metadata to its path.", "VAL", "GENERIC.MissingPubInfo"))
    return found


# ------------------------------------------------------------- sequence ----

def check_sequence(genome: str) -> list[dict]:
    found = []
    sequence = genome.upper()
    invalid = sorted(set(sequence) - IUPAC_DNA)
    if invalid:
        found.append(finding("INVALID_CHARACTERS", "ERROR", f"Sequence contains characters that are not DNA codes: {''.join(invalid)}.",
                             "Remove them; only A, C, G, T and IUPAC ambiguity codes are allowed.", "SEQ", "SEQ_INST.InvalidResidue"))
    if len(sequence) < 200:
        found.append(finding("SHORT_SEQUENCE", "ERROR", f"Sequence is {len(sequence)} bp; GenBank removes sequences under 200 bp from genome submissions.",
                             "Submit the assembled genome, not short fragments.", "DISC", "SHORT_SEQUENCES"))
    if sequence[:1] == "N" or sequence[-1:] == "N":
        found.append(finding("TERMINAL_NS", "ERROR", "The sequence begins or ends with N.", "Trim terminal Ns.", "VAL", "SEQ_INST.TerminalNs"))
    for match in re.finditer(r"N{100,}", sequence):
        found.append(finding("INTERNAL_N_RUN", "ERROR", f"A run of {len(match.group())} Ns at {match.start() + 1}-{match.end()}.",
                             "Close the gap, or annotate it as an assembly_gap with linkage evidence; do not delete Ns.",
                             "VAL", "SEQ_INST.InternalNsInSeqRaw", start=match.start() + 1, end=match.end()))
    ambiguous = sum(1 for base in sequence if base not in "ACGT")
    if sequence and ambiguous / len(sequence) > 0.5:
        found.append(finding("MOSTLY_AMBIGUOUS", "ERROR", "More than half of the sequence is ambiguous.", "Improve the assembly.", "SEQ", "DISC_PERCENT_N"))
    elif sequence and ambiguous / len(sequence) > 0.05:
        found.append(finding("HIGH_PERCENT_N", "WARNING", f"{ambiguous / len(sequence):.1%} of bases are ambiguous.",
                             "Check assembly quality before submitting.", "DISC", "DISC_PERCENT_N"))
    missing = [base for base in "ACGT" if base not in sequence]
    if sequence and missing:
        found.append(finding("ZERO_BASECOUNT", "WARNING", f"The sequence contains no {', '.join(missing)}.",
                             "This usually indicates a problem with the sequence.", "DISC", "ZERO_BASECOUNT"))
    return found


# --------------------------------------------------------- product names ----

_DB_ID = re.compile(r"\b(?:PF\d{5}|COG\d{4}|IPR\d{6}|TIGR\d{5}|VOG\d+|phrog_?\d+|GO:\d{7}|PM_\d{6}|[A-Z][A-Z0-9]{2,11}_\d{3,})\b", re.I)
_DUF = re.compile(r"\bDUF\d+\b")
_EC = re.compile(r"\bEC[ :]?\d|\b\d+\.\d+\.\d+\.(?:\d+|-)\b")


def check_product(name: str) -> list[tuple[str, str, str]]:
    """(code, severity, message) for one product name, per NOMEN and DISC SUSPECT_PRODUCT_NAMES."""
    text = (name or "").strip()
    lower = text.lower()
    issues = []
    if not text:
        return [("EMPTY_PRODUCT", "ERROR", "CDS has no product name; every CDS needs one.")]
    if lower in {"hypothetical protein", "uncharacterized protein"}:
        return []
    if re.search(r"hypothetical|uncharacteri[sz]ed|unknown (?:function|protein)|protein of unknown function", lower):
        issues.append(("HYPOTHETICAL_VARIANT", "ERROR", f"'{text}': use exactly 'hypothetical protein' when the function is unknown."))
    if re.search(r"\b(?:predicted|possible|probable|potential|candidate)\b", lower):
        issues.append(("QUALIFIER_WORDING", "ERROR", f"'{text}': use 'putative' instead of predicted/possible/probable/potential/candidate."))
    if re.search(r"\b(?:homolog(?:ue)?|ortholog|paralog)\b", lower):
        issues.append(("EVOLUTIONARY_IMPLICATION", "ERROR", f"'{text}': replace 'homolog' with '-like protein', 'putative X' or 'X'."))
    if _DB_ID.search(text) or (_DUF.search(text) and not re.search(r"domain-containing protein|family protein", lower)):
        issues.append(("DATABASE_ID_IN_PRODUCT", "ERROR", f"'{text}': remove database identifiers, accessions or locus tags from the product name."))
    if _EC.search(text):
        issues.append(("EC_NUMBER_IN_PRODUCT", "ERROR", f"'{text}': move the EC number to an EC_number qualifier."))
    if re.search(r"\b(?:partial|fragment|truncated|frameshift(?:ed)?|pseudogene|pseudo)\b", lower):
        issues.append(("SUSPICIOUS_PHRASE", "ERROR", f"'{text}': partial/fragment/truncated/frameshift do not belong in a product name; consider a pseudo gene with a note."))
    if re.search(r"[%@=\\;]|\.$", text) or '"' in text:
        issues.append(("PUNCTUATION", "ERROR", f"'{text}': remove %, @, =, \\, ;, quotes and trailing periods."))
    if "," in text or re.search(r"(?<![A-Za-z]):(?![A-Za-z])|\w: ", text):
        issues.append(("PUNCTUATION", "WARNING", f"'{text}': commas and colons are allowed only inside chemical or enzyme names."))
    if re.search(r"\bgene\b", lower):
        issues.append(("GENE_IN_PRODUCT", "WARNING", f"'{text}': name the protein, not the gene."))
    if re.search(r"\b\d+(?:\.\d+)?\s*k?da\b", lower):
        issues.append(("MOLECULAR_WEIGHT", "WARNING", f"'{text}': avoid molecular weights in product names."))
    if re.search(r"mitochondri|chloroplast|golgi", lower):
        issues.append(("ORGANELLE_NAME", "WARNING", f"'{text}': organelle names do not belong in phage protein names."))
    if re.search(r"sulph|haem|colour|characteris|\bfibre\b", lower):
        issues.append(("SPELLING", "WARNING", f"'{text}': use American spelling."))
    if re.search(r"\bprotein protein\b|ase protein\b", lower):
        issues.append(("REDUNDANT_PROTEIN", "WARNING", f"'{text}': do not repeat 'protein' or add it after '-ase'."))
    if re.match(r"[A-Z][a-z]", text) and not re.match(r"[A-Z][a-z]{2}[A-Z]", text):
        issues.append(("CAPITALIZATION", "WARNING", f"'{text}': start with lowercase unless the first word is an acronym or gene symbol."))
    if len(text) > 100:
        issues.append(("DESCRIPTIVE_NAME", "WARNING", f"Product name is {len(text)} characters; use a concise protein name and put description in a note."))
    return issues


# -------------------------------------------------------------- features ----

def check_features(genome: str, proteins: list, products: dict[str, str]) -> list[dict]:
    found = []
    length = len(genome)
    for protein in proteins:
        where = {"protein_id": protein.protein_id, "start": protein.start, "end": protein.end}
        for code, severity, message in check_product(products.get(protein.protein_id, "")):
            found.append(finding(code, severity, message, "Revise the product name (see NCBI protein naming guidelines).",
                                 "NOMEN", "SUSPECT_PRODUCT_NAMES", **where))
        five_at_end = protein.start == 1 if protein.strand == "+" else protein.end == length
        three_at_end = protein.end == length if protein.strand == "+" else protein.start == 1
        if (protein.partial_5prime and not five_at_end) or (protein.partial_3prime and not three_at_end):
            found.append(finding("INTERNAL_PARTIAL_CDS", "ERROR", "A CDS is marked partial but does not reach the end of the sequence.",
                                 "In prokaryote-like genomes, partial CDSs are allowed only at a sequence end or gap. Make it complete, extend it, or mark the gene pseudo.",
                                 "DISC", "PARTIAL_CDS_COMPLETE_SEQUENCE", **where))
        if not protein.partial_5prime and len(protein.cds) >= 3 and protein.cds[:3].upper() not in TABLE11_STARTS:
            found.append(finding("NO_START_CODON", "ERROR", f"Complete CDS begins with {protein.cds[:3]}, which is not a start codon in translation table 11.",
                                 "Correct the start, mark the 5' end partial if it is at a sequence end, or mark the gene pseudo.",
                                 "VAL", "SEQ_FEAT.StartCodon", **where))
        if not (protein.partial_5prime or protein.partial_3prime) and len(protein.sequence or "") < 50:
            found.append(finding("SHORT_PROTEIN", "WARNING", f"Protein is {len(protein.sequence or '')} aa (under 50).",
                                 "Confirm it is a real gene; NCBI flags proteins under 50 aa.", "DISC", "SHORT_PROT_SEQUENCES", **where))
    ordered = sorted(proteins, key=lambda p: (p.start, -p.end))
    for i, outer in enumerate(ordered):
        for inner in ordered[i + 1:]:
            if inner.start > outer.end:
                break
            if outer.start <= inner.start and inner.end <= outer.end and (inner.start, inner.end) != (outer.start, outer.end):
                found.append(finding("CONTAINED_CDS", "WARNING", f"{inner.protein_id} lies entirely within {outer.protein_id}.",
                                     "Remove CDSs that do not encode real proteins.", "DISC", "CONTAINED_CDS",
                                     protein_id=inner.protein_id, start=inner.start, end=inner.end))
            elif inner.start <= outer.end:
                a, b = products.get(outer.protein_id, ""), products.get(inner.protein_id, "")
                if a and a == b and a not in {"hypothetical protein", "uncharacterized protein"}:
                    found.append(finding("OVERLAPPING_SIMILAR_CDS", "WARNING", f"{outer.protein_id} and {inner.protein_id} overlap and share the product '{a}'.",
                                         "This often indicates a frameshifted gene; if so, annotate one pseudo gene across the span with a note.",
                                         "DISC", "OVERLAPPING_CDS", protein_id=inner.protein_id, start=inner.start, end=inner.end))
    return found


def review(genome: str, proteins: list, products: dict[str, str], metadata) -> dict:
    """Run every rule. Returns {"errors": [...], "warnings": [...]}."""
    findings = check_metadata(metadata) + check_sequence(genome) + check_features(genome, proteins, products)
    return {"errors": [f for f in findings if f["severity"] == "ERROR"],
            "warnings": [f for f in findings if f["severity"] != "ERROR"]}
