from __future__ import annotations

from .models import Protein

CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L","ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V","TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P","ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A","TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E","TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R","AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G"}


def translate(cds: str) -> str:
    return "".join(CODON_TABLE.get(cds[i:i+3], "X") for i in range(0, len(cds) - 2, 3)).rstrip("*")


def predict_orfs(genome_id: str, sequence: str, min_orf_nt: int = 90) -> list[Protein]:
    """Deliberately simple forward-strand demonstration caller, not production gene prediction."""
    proteins: list[Protein] = []
    stops = {"TAA", "TAG", "TGA"}
    for frame in range(3):
        start: int | None = None
        for i in range(frame, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if start is None and codon == "ATG":
                start = i
            elif start is not None and codon in stops:
                cds = sequence[start:i+3]
                if len(cds) >= min_orf_nt:
                    proteins.append(Protein(genome_id, "", start + 1, i + 3, "+", cds, translate(cds), "simple_orf_demo"))
                start = None
    proteins.sort(key=lambda protein: protein.start)
    for index, protein in enumerate(proteins, 1):
        protein.protein_id = f"PM_{index:06d}"
    return proteins
