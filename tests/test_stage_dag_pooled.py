import unittest

from phagemine.models import EvidenceLevel, Protein
from phagemine.pooled import deduplicate_proteins, remap_evidence
from phagemine.stage_dag import downstream, stale_stages, stage_record


def protein(genome, pid, sequence):
    return Protein(protein_id=pid, genome_id=genome, start=1, end=len(sequence) * 3,
                   strand="+", cds="ATG" * len(sequence), sequence=sequence,
                   gene_call_source="TEST",
                   annotation_level=EvidenceLevel.COMPUTATIONAL)


class StageAndPooledTests(unittest.TestCase):
    def test_swissprot_invalidation_is_downstream_only(self):
        self.assertEqual(downstream("SWISSPROT"), {"SWISSPROT", "EVIDENCE_FUSION", "RANKING_MINING", "QC_REPORTING", "GENBANK"})

    def test_changed_resource_marks_only_dependents(self):
        previous = {"inputs": {"GENE_PREDICTION": {"input_checksums": {"g": "1"}}, "SWISSPROT": {"resource_checksums": {"db": "old"}}}}
        current = {"GENE_PREDICTION": {"input_checksums": {"g": "1"}}, "SWISSPROT": {"resource_checksums": {"db": "new"}}}
        stale = stale_stages(previous, current)
        self.assertIn("SWISSPROT", stale)
        self.assertNotIn("PFAM", stale)

    def test_exact_deduplication_and_lossless_remap(self):
        proteins = [protein("G1", "P1", "MKK"), protein("G2", "P2", "MKK"), protein("G3", "P3", "MNN")]
        reps, occurrences = deduplicate_proteins(proteins)
        self.assertEqual([p.protein_id for p in reps], ["P1", "P3"])
        mapped = remap_evidence({"P1": [{"source": "Pfam", "supports": True, "provenance": {}}]}, reps, occurrences)
        self.assertEqual(set(mapped), {"P1", "P2"})
        self.assertTrue(all(e["provenance"]["pooled"] for e in mapped["P2"]))

    def test_stage_record_is_explicit(self):
        row = stage_record("PFAM", status="SUCCESS_WITH_HITS", input_checksums={"proteins": "x"})
        self.assertEqual(row["dependencies"], ["GENE_PREDICTION"])
        self.assertIn("output_checksums", row)


if __name__ == "__main__":
    unittest.main()
