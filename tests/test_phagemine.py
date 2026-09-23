import json
import csv
import hashlib
import os
import stat
import shutil
import tempfile
import unittest
import pytest
from io import StringIO
from unittest.mock import patch
from pathlib import Path

from phagemine.genome import predict_orfs, translate
from phagemine.gene_prediction import DemoORFPredictor, PHANOTATEPredictor
from phagemine.io import read_fasta
from phagemine.models import Evidence, EvidenceLevel, Protein
from phagemine.pipeline import _evidence_progress_summary, run
from phagemine.genbank import feature_table, product_name, table2asn_status, validate, write_package
from phagemine.annotation import MockEvidenceBackend
from phagemine.cli import main
from phagemine.models import SubmissionMetadata
from phagemine.genome_representation import GenomeRepresentation, Orientation, Rotation, Topology
from phagemine.sequencing_provenance import SequencingPlatform, SequencingProvenance
from phagemine.pfam import PfamHMMAdapter
from phagemine.vog import VOGHMMAdapter
from phagemine.phrogs import MMSEQS_FORMAT, PHROGSMMseqsAdapter
from phagemine.swissprot import SwissProtEvidenceAdapter
from phagemine.progress import ProgressReporter
from phagemine.resume import _load_source, RESUME_STAGES, checkpoint_reusable, reclassify
from phagemine.fusion import attach_gene_call_assessments, classify_protein, classify_proteins, normalize_function, write_classification
from phagemine.context import build_context, write_context
from phagemine.compare import compare
from phagemine.batch import discover_inputs, batch
from phagemine.reconciliation import GeneModel, gene_call_review, reconcile_models
from phagemine.adjudication import adjudicate
from phagemine.alternative_evidence import alternative_models, extract_translation, cache_key
from phagemine.pfam import PfamHMMAdapter
from phagemine.vog import VOGHMMAdapter
from phagemine.swissprot import SwissProtEvidenceAdapter
from phagemine.phrogs import PHROGSMMseqsAdapter
from phagemine.benchmark import benchmark, import_phold, import_prokka, import_pharokka, import_phagemine, compare_models, metrics, truth_metrics, classify_product_relation
from phagemine.reporting import write_checkpoint_snapshot, write_stage_checkpoint, prefix_checkpoint_artifacts, update_comparative_report
from phagemine.evidence import EvidenceAdapterResult
from phagemine.mining import mine
from phagemine.resources import EvidenceResourceManager, ResourceStatus, ResourceType, default_registry_path
from phagemine.hallmarks import assess_hallmarks
from phagemine.review import build_annotation_review


ROOT = Path(__file__).resolve().parents[1]
RESUME_SOURCE_FIXTURE = ROOT / "tests" / "fixtures" / "resume_source"


class PhageMineTests(unittest.TestCase):
    def test_hoc_like_product_is_preserved_with_separate_display_note(self):
        result = classify_protein(self._fusion_protein([
            self._fusion_e("PHROGs", "Hoc-like head decoration", identifier="phrog_2973")
        ]))
        self.assertEqual(result["display_product"], "hoc-like head decoration protein")
        self.assertIn("capsid-display candidate", result["biotechnology_relevance"])

    def test_major_head_synonym_normalizes_to_major_capsid(self):
        result = classify_protein(self._fusion_protein([
            self._fusion_e("PHROGs", "major head protein", identifier="phrog_247")
        ]))
        self.assertEqual(result["display_product"], "major capsid protein")

    def test_primary_html_receives_completed_inphared_numerical_taxonomy(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); (root / "report.html").write_text("<html><body><h1>Report</h1></body></html>")
            update_comparative_report(root, {"inphared": {"unique_reference_summaries": [{
                "reference_description": "Reference phage", "reference_accessions": "AB123",
                "host_genus": "Salmonella", "phage_family": "Sarkviridae", "phage_genus": "Jerseyvirus",
                "intergenomic_similarity_percent": 96.2, "query_aligned_percent": 99.0,
                "reference_aligned_percent": 98.0, "taxonomic_interpretation": "CONSISTENT_WITH_SAME_SPECIES_THRESHOLD",
            }]}})
            report = (root / "report.html").read_text()
            self.assertIn("INPHARED nearest-reference nucleotide comparison", report)
            self.assertIn("96.20%", report)
            self.assertIn("Mash distance is not converted to similarity", report)

            # Re-applying the same comparative result must replace rather than
            # duplicate the INPHARED section.
            update_comparative_report(root, {"inphared": {"unique_reference_summaries": [{
                "reference_description": "Reference phage", "reference_accessions": "AB123",
                "host_genus": "Salmonella", "phage_family": "Sarkviridae", "phage_genus": "Jerseyvirus",
                "intergenomic_similarity_percent": 96.2, "query_aligned_percent": 99.0,
                "reference_aligned_percent": 98.0,
                "taxonomic_interpretation": "CONSISTENT_WITH_SAME_SPECIES_THRESHOLD",
            }]}})
            report = (root / "report.html").read_text()
            self.assertEqual(
                report.count("id='inphared-numerical-taxonomy'"),
                1,
            )

    def test_orf_reconciliation_boundary_and_strand_categories(self):
        p=GeneModel("PHANOTATE","P1",100,400,"+")
        self.assertEqual(reconcile_models([p],[GeneModel("Prodigal","D1",100,400,"+")])[0]["conflict_type"],"EXACT_CONCORDANCE")
        self.assertEqual(reconcile_models([p],[GeneModel("Prodigal","D1",130,400,"+")])[0]["conflict_type"],"START_DISCORDANCE")
        self.assertEqual(reconcile_models([p],[GeneModel("Prodigal","D1",100,370,"+")])[0]["conflict_type"],"STOP_DISCORDANCE")
        self.assertEqual(reconcile_models([p],[GeneModel("Prodigal","D1",130,370,"+")])[0]["conflict_type"],"START_AND_STOP_DISCORDANCE")
        self.assertEqual(reconcile_models([p],[GeneModel("Prodigal","D1",100,400,"-")])[0]["conflict_type"],"STRAND_DISCORDANCE")

    def test_orf_reconciliation_caller_specific_and_deterministic(self):
        p=GeneModel("PHANOTATE","P1",100,400,"+"); d=GeneModel("Prodigal","D1",800,1000,"+")
        rows=reconcile_models([p],[d]); self.assertEqual({r["conflict_type"] for r in rows},{"PHANOTATE_ONLY","PRODIGAL_ONLY"})
        self.assertEqual(rows,reconcile_models([p],[d]))

    def test_orf_reconciliation_handles_unsorted_prodigal_models(self):
        phanotate=[GeneModel("PHANOTATE","P1",100,400,"+")]
        prodigal=[GeneModel("Prodigal","D2",800,1000,"+"),GeneModel("Prodigal","D1",100,400,"+")]
        rows=reconcile_models(phanotate,prodigal)
        self.assertEqual(rows[0]["prodigal_id"],"D1")
        self.assertEqual([row["prodigal_id"] for row in rows if row["conflict_type"]=="PRODIGAL_ONLY"],["D2"])

    def test_gene_call_review_flags_short_unsupported_caller_specific_orf(self):
        rows=reconcile_models([GeneModel("PHANOTATE","P1",100,180,"+")],[])
        review=gene_call_review(rows, {"P1": []}, {"P1": 27})[0]
        self.assertEqual(review["gene_call_confidence"], "LOW")
        self.assertEqual(review["review_flag"], "POSSIBLE_FALSE_CALL")

    def test_gene_call_review_preserves_exact_concordant_orf(self):
        rows=reconcile_models([GeneModel("PHANOTATE","P1",100,400,"+")],[GeneModel("Prodigal","D1",100,400,"+")])
        review=gene_call_review(rows, {}, {"P1": 100})[0]
        self.assertEqual(review["gene_call_confidence"], "HIGH")
        self.assertEqual(review["review_flag"], "NONE")

    def test_gene_call_fields_remain_separate_from_functional_review(self):
        classification = classify_protein(self._fusion_protein())
        attached = attach_gene_call_assessments([classification], [{
            "protein_id": "P", "gene_call_confidence": "HIGH", "review_flag": "NONE",
        }])[0]
        self.assertEqual(attached["gene_call_confidence"], "HIGH")
        self.assertEqual(attached["gene_call_review_flag"], "NONE")
        self.assertEqual(attached["functional_review_flag"], "REVIEW_REQUIRED")

    def test_exact_concordance_adjudication_confirms_both_boundaries_and_strand(self):
        row = reconcile_models(
            [GeneModel("PHANOTATE", "P1", 100, 400, "+")],
            [GeneModel("Prodigal", "D1", 100, 400, "+")],
        )[0]
        result = adjudicate([row], {
            "P1": [{"supports": True, "evidence_strength": "STRONG"}],
            "D1": [],
        })[0]
        self.assertEqual(result["start_decision"], "START_CONFIRMED")
        self.assertEqual(result["stop_decision"], "STOP_CONFIRMED")
        self.assertEqual(result["strand_decision"], "STRAND_CONFIRMED")
        self.assertEqual(result["evidence_agreement"], "COORDINATE_CONCORDANCE")
        self.assertFalse(result["manual_review"])

    def test_hallmark_check_never_converts_non_detection_to_absence(self):
        rows=assess_hallmarks([{"protein_id":"P1","functional_state":"PROBABLE_FUNCTION","display_product":"major capsid protein","confidence":"MODERATE"}])
        by_name={row["hallmark"]:row for row in rows}
        self.assertEqual(by_name["major_capsid"]["status"],"DETECTED")
        self.assertEqual(by_name["portal"]["status"],"NOT_ESTABLISHED")
        self.assertIn("not evidence of biological absence",by_name["portal"]["interpretation"])

    def test_annotation_review_combines_function_gene_call_and_hallmark_flags(self):
        protein=self._fusion_protein(); classification=classify_protein(protein)
        review=build_annotation_review([protein],[classification],[{"review_flag":"POSSIBLE_FALSE_CALL","protein_id":"P","phanotate_id":"P","gene_call_confidence":"LOW","rationale":"short unsupported"}],[{"hallmark":"major_capsid","status":"NOT_ESTABLISHED","interpretation":"not established"}])
        self.assertEqual({row["review_type"] for row in review},{"FUNCTION_ASSIGNMENT","GENE_CALL","HALLMARK_NOT_ESTABLISHED"})

    def test_orf_adjudication_is_observational_and_model_specific(self):
        row=reconcile_models([GeneModel("PHANOTATE","P1",100,400,"+")],[GeneModel("Prodigal","D1",130,400,"+")])[0]
        result=adjudicate([row], {"P1":[{"supports":True,"evidence_strength":"STRONG","source":"Swiss-Prot"}],"D1":[]})[0]
        self.assertEqual(result["start_decision"],"START_AMBIGUOUS")
        self.assertEqual(result["adjudication_status"],"INSUFFICIENT_COMPARATIVE_EVIDENCE")
        self.assertTrue(result["manual_review"])

    def test_alternative_translation_is_strand_aware_and_model_specific(self):
        cds, protein=extract_translation("ATGAAATAG",1,9,"+")
        self.assertEqual(cds,"ATGAAATAG"); self.assertEqual(protein,"MK")
        rows=[{"locus_id":"L1","conflict_type":"START_DISCORDANCE","prodigal_id":"D1","prodigal_start":1,"prodigal_end":9,"prodigal_strand":"+"}]
        models=alternative_models(rows,"ATGAAATAG")
        self.assertEqual(models[0]["alternative_model_id"],"ALT_PRODIGAL_00001")
        self.assertNotEqual(models[0]["translation_sha256"],"")

    def test_alternative_cache_identity_supports_all_adapters(self):
        model={"translation_sha256":"abc"}
        adapters=[PfamHMMAdapter(evalue_threshold=1e-5),VOGHMMAdapter(evalue_threshold=1e-5),SwissProtEvidenceAdapter(evalue_threshold=1e-5),PHROGSMMseqsAdapter(evalue_threshold=1e-5)]
        keys=[cache_key(model,a) for a in adapters]
        self.assertEqual(keys,[cache_key(model,a) for a in adapters]); self.assertEqual(len(set(keys)),4)
        self.assertNotEqual(cache_key(model,PfamHMMAdapter(evalue_threshold=1e-3)),keys[0])

    def test_batched_alternative_evidence_maps_protein_ids_without_model_attribute(self):
        from phagemine.alternative_evidence import acquire_alternative_evidence
        class Adapter:
            name='TEST'
            def provenance(self): return {'adapter':'TEST','thresholds':{'x':1}}
            def analyze(self, proteins):
                from phagemine.evidence import EvidenceAdapterResult
                from phagemine.models import Evidence, EvidenceLevel
                hits=[]
                for p in proteins:
                    if p.protein_id != 'ALT_PRODIGAL_00002':
                        e=Evidence('TEST','statement',EvidenceLevel.COMPUTATIONAL,'TEST','1',supports=True,identifier=p.protein_id); e.provenance['protein_id']=p.protein_id; hits.append(e)
                return EvidenceAdapterResult('TEST','REAL',evidence=hits)
        models=[{'alternative_model_id':f'ALT_PRODIGAL_{i:05d}','locus_id':f'L{i}','caller_id':str(i),'start':1,'end':9,'strand':'+','cds':'ATGAAATAG','protein_sequence':'MK','translation_sha256':str(i),'evidence':[]} for i in range(1,4)]
        acquire_alternative_evidence(models,[Adapter()])
        self.assertEqual(len(models[0]['evidence']),1); self.assertEqual(len(models[1]['evidence']),0); self.assertEqual(len(models[2]['evidence']),1)

    def test_alternative_search_status_counts_match_model_evidence(self):
        from phagemine.alternative_evidence import acquire_alternative_evidence
        from phagemine.evidence import EvidenceAdapterResult
        from phagemine.models import Evidence, EvidenceLevel
        class Adapter:
            name='TEST'
            def provenance(self): return {'adapter':'TEST','thresholds':{}}
            def analyze(self, proteins):
                hits=[]
                for p in proteins:
                    n=2 if p.protein_id == 'ALT_PRODIGAL_00001' else (1 if p.protein_id == 'ALT_PRODIGAL_00003' else 0)
                    for i in range(n):
                        e=Evidence('TEST','statement',EvidenceLevel.COMPUTATIONAL,'TEST',str(i),supports=(p.protein_id != 'ALT_PRODIGAL_00003'))
                        e.provenance['protein_id']=p.protein_id; hits.append(e)
                return EvidenceAdapterResult('TEST','REAL',evidence=hits)
        models=[{'alternative_model_id':f'ALT_PRODIGAL_{i:05d}','locus_id':f'L{i}','caller_id':str(i),'start':1,'end':9,'strand':'+','cds':'ATGAAATAG','protein_sequence':'MK','translation_sha256':str(i),'evidence':[]} for i in range(1,4)]
        with tempfile.TemporaryDirectory() as td:
            status=Path(td)/'status.tsv'; acquire_alternative_evidence(models,[Adapter()],status_output=status)
            rows=list(csv.DictReader(status.open(),delimiter='\t'))
        self.assertEqual([int(r['records_returned']) for r in rows],[2,0,1])
        self.assertEqual([int(r['accepted_records']) for r in rows],[2,0,0])
        self.assertEqual(rows[1]['search_status'],'SEARCH_EXECUTED_ZERO_HITS')
        self.assertEqual(sum(int(r['records_returned']) for r in rows),3)

    def test_fresh_alternative_status_uses_canonical_source_partition(self):
        """Fresh batched results with legacy source spelling still balance."""
        from phagemine.alternative_evidence import acquire_alternative_evidence
        from phagemine.evidence import EvidenceAdapterResult
        from phagemine.models import Evidence, EvidenceLevel
        class Adapter:
            name='SwissProt'
            def provenance(self): return {'thresholds': {'x': 1}}
            def analyze(self, proteins):
                hits=[]
                for p in proteins:
                    e=Evidence('sequence_similarity','hit',EvidenceLevel.CURATED,
                               'SwissProt','1',supports=True,identifier=p.protein_id)
                    e.provenance['protein_id']=p.protein_id; hits.append(e)
                return EvidenceAdapterResult('SwissProt','REAL',evidence=hits)
        models=[{'alternative_model_id':f'ALT_PRODIGAL_{i:05d}','locus_id':f'L{i}',
                 'caller_id':str(i),'start':1,'end':9,'strand':'+','cds':'ATGAAATAG',
                 'protein_sequence':'MK','translation_sha256':str(i),'evidence':[]} for i in range(1,3)]
        with tempfile.TemporaryDirectory() as td:
            status=Path(td)/'status.tsv'; acquire_alternative_evidence(models,[Adapter()],status_output=status)
            rows=list(csv.DictReader(status.open(),delimiter='\t'))
        self.assertEqual(sum(int(r['records_returned']) for r in rows),2)
        self.assertEqual(sum(int(r['accepted_records']) for r in rows),2)

    def test_pmf_comparison_statistics_pairwise_core_and_presence(self):
        from phagemine.family_compare import compare_database
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)/'db'; root.mkdir()
            families=[{'family_id':'PMF-000001','family_database_version':'v1'},{'family_id':'PMF-000002','family_database_version':'v1'}]
            (root/'families.json').write_text(json.dumps(families))
            with (root/'family_members.tsv').open('w') as h:
                h.write('family_id\tmember_id\tsequence_sha256\tsource_genome_id\tsource_protein_id\thost_genus\toriginal_annotation\n')
                h.write('PMF-000001\tA|P1\tx\tA\tA|P1\tPseudomonas\t\nPMF-000001\tB|P1\ty\tB\tB|P1\tSalmonella\t\nPMF-000002\tA|P2\tz\tA\tA|P2\tPseudomonas\t\n')
            summary=compare_database(root,Path(td)/'out',core_genomes=['A','B'])
            self.assertEqual(summary['statistics']['total_proteins'],3)
            self.assertEqual(summary['statistics']['singleton_families'],1)
            self.assertEqual(summary['core_family_count'],1)
            pair=list(csv.DictReader((Path(td)/'out'/'pmf_pairwise_genomes.tsv').open(),delimiter='\t'))[0]
            self.assertEqual(pair['shared_pmf_count'],'1')
            self.assertTrue((Path(td)/'out'/'pmf_presence_absence.tsv').exists())

    def test_pmf_sensitivity_rows_preserve_effective_configuration(self):
        from phagemine.family_compare import sensitivity_rows
        summary={'statistics':{'total_families':3,'singleton_families':1,'multi_member_families':2,'family_size_distribution':{'1':1,'2':2}},'minimum_coverage':0.5,'coverage_mode':0,'clustering_mode':1,'backend':'MMSEQS2','core_family_count':2,'host_group_counts':{'CROSS_HOST_COHORT':1}}
        row=sensitivity_rows({0.2:summary})[0]
        self.assertEqual(row['minimum_identity'],0.2); self.assertEqual(row['minimum_coverage'],0.5)
        self.assertEqual(row['coverage_mode'],0); self.assertEqual(row['clustering_mode'],1); self.assertEqual(row['backend'],'MMSEQS2')

    def test_pmf_enrichment_conservative_states(self):
        from phagemine.family_enrichment import enrich_database
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); db=root/'db'; db.mkdir(); result=root/'A'; result.mkdir(); out=root/'out'
            (db/'families.json').write_text(json.dumps([{'family_id':'PMF-000001','member_count':1,'genome_count':1,'host_genera':['Pseudomonas'],'source_genome_ids':['A']}]))
            (db/'family_members.tsv').write_text('family_id\tmember_id\tsequence_sha256\tsource_genome_id\tsource_protein_id\thost_genus\toriginal_annotation\nPMF-000001\tA|P1\tx\tA\tA|P1\tPseudomonas\t\n')
            (result/'evidence.json').write_text(json.dumps([{'protein_id':'P1','annotation':'hypothetical protein','evidence':[]}]))
            (result/'functional_classification.json').write_text(json.dumps([{'protein_id':'P1','functional_state':'CONSERVED_UNKNOWN','proposed_function':None}]))
            rows=enrich_database(db,[result],out)
            self.assertEqual(rows[0]['family_functional_status'],'UNKNOWN_SINGLETON')
            self.assertTrue((out/'pmf_member_evidence.tsv').exists())

    def test_pmf_priority_only_ranks_unknown_multi_genome(self):
        from phagemine.family_priority import prioritize
        with tempfile.TemporaryDirectory() as td:
            db=Path(td); (db/'pmf_functional_enrichment.tsv').write_text('family_id\tfamily_functional_status\tmember_count\tgenome_count\tsource_genomes\thost_genera\tevidence_sources_present\tmember_annotations\nA\tUNKNOWN_MULTI_GENOME\t2\t2\tA,B\tPseudomonas\tPfam,VOGDB\t\nB\tUNKNOWN_SINGLETON\t1\t1\tA\tPseudomonas\t\t\nC\tCONFLICTING_MEMBER_ANNOTATIONS\t2\t2\tA,B\tPseudomonas\tPfam\tx;y\n')
            out=Path(td)/'out'; rows=prioritize(db,out)
            self.assertEqual([r['family_id'] for r in rows],['A'])
            self.assertNotIn('novel', (out/'pmf_unresolved_priority.tsv').read_text().lower())

    def test_benchmark_import_and_coordinate_metrics(self):
        with tempfile.TemporaryDirectory() as temp:
            g=Path(temp)/"x.gff"; g.write_text("##gff-version 3\ng\tProkka\tCDS\t100\t400\t.\t+\t0\tID=x1;product=hypothetical protein\n")
            a=import_prokka(g); b=import_pharokka(g); self.assertEqual(compare_models(a,b)[0]["relationship"],"EXACT_MATCH"); self.assertEqual(metrics(a,b)["exact_f1"],1.0)

    def test_benchmark_imports_phold_genbank_and_compares_products(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); gbk=root/"phold.gbk"
            gbk.write_text('LOCUS       test\nFEATURES             Location/Qualifiers\n     CDS             complement(100..400)\n                     /locus_tag="P1"\n                     /product="major capsid protein"\n                     /translation="MPEP\n                     TIDE"\nORIGIN\n//\n')
            phold=import_phold(gbk)
            self.assertEqual((phold[0]["start"],phold[0]["end"],phold[0]["strand"]),(100,400,"-"))
            self.assertEqual(phold[0]["product"],"major capsid protein")
            self.assertEqual(phold[0]["protein_length"],8)
            other=[{**phold[0],"tool":"PHAGEMINE","locus_id":"PM_1","product":"putative major capsid protein"}]
            benchmark({"PHAGEMINE":other,"Phold":phold},root/"comparison")
            rows=list(csv.DictReader((root/"comparison"/"functional_comparison.tsv").open(),delimiter="\t"))
            self.assertEqual(rows[0]["status"],"EXACT_PRODUCT_AGREEMENT")

    def test_benchmark_exact_match_precedes_boundary_contact(self):
        a=[{'start':2,'end':313,'strand':'-','locus_id':'A'},{'start':313,'end':894,'strand':'-','locus_id':'B'}]
        b=[{'start':313,'end':894,'strand':'-','locus_id':'X'}]
        rows=compare_models(a,b)
        exact=[r for r in rows if r['relationship']=='EXACT_MATCH']
        self.assertEqual(len(exact),1); self.assertEqual(exact[0]['method_a']['locus_id'],'B')
        self.assertFalse(any(r['relationship']=='START_AND_STOP_DIFFERENCE' and r.get('method_b') for r in rows))

    def test_benchmark_imports_phagemine_and_fails_on_missing_required_files(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); (root/"genes.gff3").write_text("g\tPhageMine\tCDS\t1\t9\t.\t+\t0\tID=PM_000001\n"); (root/"annotation.tsv").write_text("protein_id\tannotation\nPM_000001\tunknown\n"); (root/"functional_classification.tsv").write_text("protein_id\tfunctional_state\tproposed_function\tconfidence\nPM_000001\tUNRESOLVED\t\tNONE\n")
            self.assertEqual(import_phagemine(root)[0]["tool"],"PHAGEMINE")
            (root/"functional_classification.tsv").unlink()
            legacy=import_phagemine(root)
            self.assertEqual(legacy[0]["functional_classification_status"],"LEGACY_OUTPUT_NOT_AVAILABLE")
            (root/"genes.gff3").unlink()
            with self.assertRaises(ValueError): import_phagemine(root)

    def test_swissprot_unavailable_is_explicit(self):
        adapter=SwissProtEvidenceAdapter(database_path=None, diamond="/nonexistent/diamond")
        result=adapter.analyze([])
        self.assertEqual(result.status,"UNAVAILABLE")
        self.assertIn("unavailable", result.message.lower())

    def test_unavailable_checkpoint_is_invalidated_when_resource_is_available(self):
        checkpoint={"status":"REAL","provenance":{"status":"UNAVAILABLE","adapter":"SwissProtEvidenceAdapter"}}
        self.assertFalse(checkpoint_reusable(checkpoint,{"status":"REAL","adapter":"SwissProtEvidenceAdapter"}))
        valid={"status":"REAL","provenance":{"status":"REAL","adapter":"SwissProtEvidenceAdapter","swissprot_version":"v1","thresholds":{"evalue":1e-5}}}
        self.assertTrue(checkpoint_reusable(valid,{"status":"REAL","adapter":"SwissProtEvidenceAdapter","swissprot_version":"v1","thresholds":{"evalue":1e-5}}))
        self.assertFalse(checkpoint_reusable(valid,{"status":"REAL","adapter":"SwissProtEvidenceAdapter","swissprot_version":"v2","thresholds":{"evalue":1e-5}}))
    def test_batch_discovers_deterministically_and_runs_isolated_samples(self):
        with tempfile.TemporaryDirectory() as temp:
            inp=Path(temp)/"in"; out=Path(temp)/"out"; inp.mkdir()
            genome=(ROOT/"examples/demo_phage.fasta").read_text(); (inp/"z.fa").write_text(genome.replace("demo_phage","z")); (inp/"a.fasta").write_text(genome.replace("demo_phage","a")); (inp/"skip.txt").write_text("x")
            self.assertEqual([p.name for p in discover_inputs(inp)], ["a.fasta","z.fa"])
            def fake_run(path, destination, progress=None, **kwargs):
                self.assertTrue((Path(destination) / "sample_status.json").is_file())
                progress.start("input/genome validation"); progress.finish()
                (Path(destination)/"analysis_genome.fasta").write_text(genome)
                (Path(destination)/"functional_classification.json").write_text(json.dumps([{"functional_state":"UNRESOLVED"}]))
                for name, value in (("evidence.json", []), ("genomic_context.json", []), ("modules.json", []), ("run_manifest.json", {})):
                    (Path(destination)/name).write_text(json.dumps(value))
            with patch("phagemine.batch.run", fake_run):
                rows=batch(inp,out,gene_predictor="demo")
            self.assertEqual([r["sample_id"] for r in rows], ["a","z"])
            self.assertTrue((out/"a"/"run_manifest.json").is_file()); self.assertTrue((out/"z"/"functional_classification.json").is_file())
            self.assertTrue((out/"batch_manifest.json").is_file())

    def test_batch_sample_uses_standard_single_run_scientific_outputs(self):
        with tempfile.TemporaryDirectory() as temp:
            inp, out = Path(temp)/"in", Path(temp)/"out"; inp.mkdir()
            genome=(ROOT/"examples/demo_phage.fasta").read_text(); (inp/"a.fasta").write_text(genome)
            observed_commands = []
            def fake_run(path, destination, progress=None, **kwargs):
                observed_commands.append(kwargs.get("command"))
                d=Path(destination); (d/"analysis_genome.fasta").write_text(genome)
                (d/"functional_classification.json").write_text("[]"); (d/"genomic_context.json").write_text("[]"); (d/"modules.json").write_text("[]")
                for name in ("run_manifest.json","evidence.json","proteins.faa","genes.gff3","cds.fna","annotation.tsv","functional_classification.tsv"):
                    (d/name).write_text("{}" if name.endswith("json") else "")
            with patch("phagemine.batch.run", fake_run):
                batch(inp, out, gene_predictor="demo")
            sample=out/"a"
            self.assertEqual(observed_commands, ["annotate"])
            self.assertTrue((sample/"annotation.tsv").is_file())
            self.assertFalse((sample/"a_proteins.faa").exists())

    def test_batch_propagates_threads_to_each_single_genome_run(self):
        with tempfile.TemporaryDirectory() as temp:
            inp, out = Path(temp) / "in", Path(temp) / "out"
            inp.mkdir()
            genome = (ROOT / "examples/demo_phage.fasta").read_text()
            (inp / "a.fasta").write_text(genome)
            observed = []

            def fake_run(path, destination, progress=None, **kwargs):
                observed.append(kwargs.get("threads"))
                destination = Path(destination)
                (destination / "analysis_genome.fasta").write_text(genome)
                for name, value in (
                    ("functional_classification.json", []), ("genomic_context.json", []),
                    ("modules.json", []), ("evidence.json", []), ("run_manifest.json", {}),
                ):
                    (destination / name).write_text(json.dumps(value))

            with patch("phagemine.batch.run", fake_run):
                batch(inp, out, gene_predictor="demo", threads=4)
            self.assertEqual(observed, [4])

    def test_batch_continues_after_failed_sample(self):
        with tempfile.TemporaryDirectory() as temp:
            inp=Path(temp)/"in"; out=Path(temp)/"out"; inp.mkdir()
            genome=(ROOT/"examples/demo_phage.fasta").read_text(); (inp/"a.fasta").write_text(genome.replace("demo_phage","a")); (inp/"bad.fasta").write_text(">a\nNNN\n")
            def fake_run(path, destination, progress=None, **kwargs):
                if Path(path).name == "bad.fasta":
                    raise RuntimeError("synthetic failure")
                (Path(destination)/"analysis_genome.fasta").write_text(genome)
                (Path(destination)/"functional_classification.json").write_text(json.dumps([{"functional_state":"UNRESOLVED"}]))
                for name, value in (("evidence.json", []), ("genomic_context.json", []), ("modules.json", []), ("run_manifest.json", {})):
                    (Path(destination)/name).write_text(json.dumps(value))
            with patch("phagemine.batch.run", fake_run):
                rows=batch(inp,out,gene_predictor="demo")
            self.assertEqual([r["status"] for r in rows], ["SUCCESS","FAILED"])
            self.assertTrue((out/"a"/"run_manifest.json").is_file())
            self.assertTrue((out/"bad"/"sample_status.json").is_file())

    def test_batch_persists_checkpoints_and_reuses_valid_sample(self):
        with tempfile.TemporaryDirectory() as temp:
            inp, out = Path(temp)/"in", Path(temp)/"out"; inp.mkdir()
            genome=(ROOT/"examples/demo_phage.fasta").read_text(); (inp/"a.fasta").write_text(genome)
            calls=[]
            def fake_run(path, destination, progress=None, **kwargs):
                calls.append(Path(path).name)
                progress.start("gene prediction"); progress.finish()
                d=Path(destination); (d/"analysis_genome.fasta").write_text(genome)
                (d/"functional_classification.json").write_text(json.dumps([{"functional_state":"UNRESOLVED"}]))
                for name, value in (("evidence.json", []), ("genomic_context.json", []), ("modules.json", []), ("run_manifest.json", {})):
                    (d/name).write_text(json.dumps(value))
            with patch("phagemine.batch.run", fake_run):
                batch(inp, out, gene_predictor="demo")
                rows=batch(inp, out, resume_existing=True, gene_predictor="demo")
            self.assertEqual(calls, ["a.fasta"])
            self.assertEqual(rows[0]["status"], "REUSED")
            state=json.loads((out/"a"/"sample_status.json").read_text())
            self.assertEqual(state["status"], "REUSED")
            self.assertTrue((out/"a"/"checkpoints").is_dir())

    def test_batch_recovers_failed_evidence_complete_sample(self):
        with tempfile.TemporaryDirectory() as temp:
            inp, out = Path(temp)/"in", Path(temp)/"out"; inp.mkdir()
            genome=(ROOT/"examples/demo_phage.fasta").read_text(); fasta=inp/"a.fasta"; fasta.write_text(genome)
            sample=out/"a"; (sample/"checkpoints").mkdir(parents=True); (sample/"logs").mkdir()
            digest=hashlib.sha256(fasta.read_bytes()).hexdigest()
            (sample/"sample_status.json").write_text(json.dumps({"sample_id":"a","input_sha256":digest,"status":"FAILED","failed_stage":"functional classification / evidence fusion","error_message":"synthetic","checkpoints":{"evidence_integration":{"status":"COMPLETE"}}}))
            with patch("phagemine.batch.run", side_effect=AssertionError("pipeline must not restart")), patch("phagemine.batch.recover_evidence_complete") as recover:
                def recovered(destination, progress):
                    (sample/"analysis_genome.fasta").write_text(genome)
                    (sample/"functional_classification.json").write_text("[]")
                    (sample/"modules.json").write_text("[]")
                recover.side_effect=recovered
                rows=batch(inp, out, resume_existing=True, gene_predictor="demo")
            self.assertEqual(rows[0]["status"], "SUCCESS")
            recover.assert_called_once()

    def test_evidence_complete_checkpoint_contains_recovery_artifacts(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); fasta=root/"input.fasta"; fasta.write_text(">g\nATGATGATG\n")
            protein=self._fusion_protein(); representation=GenomeRepresentation.original("g", "ATGATGATG")
            checkpoint=root/"checkpoints"/"evidence_complete"
            write_checkpoint_snapshot(root, checkpoint, representation, SequencingProvenance(), [protein], {"gene_caller":{"name":"PHANOTATE"},"evidence_adapters":[{"adapter":"x","status":"REAL"}],"input_sha256":hashlib.sha256(fasta.read_bytes()).hexdigest()}, fasta)
            self.assertTrue(all((checkpoint/name).is_file() for name in ("original_input.fasta","analysis_genome.fasta","genome_representation.json","proteins.faa","evidence.json","checkpoint_manifest.json")))

    def test_stage_checkpoint_persists_evidence_and_provenance_atomically(self):
        with tempfile.TemporaryDirectory() as temp:
            write_stage_checkpoint(Path(temp)/"checkpoints", "pfam", [{"identifier":"PF1"}], {"database":"pfam.hmm","threshold":1e-5})
            payload=json.loads((Path(temp)/"checkpoints"/"pfam"/"evidence.json").read_text())
            manifest=json.loads((Path(temp)/"checkpoints"/"pfam"/"checkpoint_manifest.json").read_text())
            self.assertEqual(payload[0]["identifier"], "PF1")
            self.assertEqual(manifest["provenance"]["threshold"], 1e-5)

    def test_batch_checkpoint_aliases_are_sample_prefixed_and_legacy_remains(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); (root/"checkpoints"/"evidence_complete").mkdir(parents=True)
            (root/"checkpoints"/"evidence_complete"/"evidence.json").write_text("[{\"protein_id\":\"P1\"}]")
            (root/"checkpoints"/"evidence_complete"/"checkpoint_manifest.json").write_text("{}")
            prefix_checkpoint_artifacts(root, "phage_a")
            self.assertTrue((root/"checkpoints/evidence_complete/phage_a_evidence.json").is_file())
            self.assertTrue((root/"checkpoints/evidence_complete/phage_a_checkpoint.json").is_file())
            self.assertTrue(json.loads((root/"checkpoints/evidence_complete/phage_a_evidence.json").read_text())[0]["sample_id"] == "phage_a")
    def test_compare_synthetic_phrog_anchor_and_mmseqs_fallback(self):
        with tempfile.TemporaryDirectory() as temp:
            roots=[]
            for suffix, phrog in (("a", True), ("b", False)):
                root=Path(temp)/suffix; root.mkdir(); pid="P1"; seq="M"*10
                (root/"proteins.faa").write_text(f">{pid}\n{seq}\n")
                (root/"genes.gff3").write_text("##gff-version 3\n")
                cls={"protein_id":pid,"functional_state":"CONSERVED_UNKNOWN","proposed_function":None,"functional_category":None,"conservation_status":"STRONGLY_CONSERVED"}
                ctx={"protein_id":pid,"gene_order_index":0,"start":1,"end":30,"strand":"+","functional_state":"CONSERVED_UNKNOWN","proposed_function":None,"functional_category":None,"conservation_status":"STRONGLY_CONSERVED","functional_module":"head_and_packaging","module_id":"m1"}
                ev={"protein_id":pid,"genome_id":"g","start":1,"end":30,"strand":"+","sequence":seq,"evidence":[{"supports":True,"source":"PHROGs","identifier":"1"}] if phrog else []}
                (root/"functional_classification.json").write_text(json.dumps([cls])); (root/"genomic_context.json").write_text(json.dumps([ctx])); (root/"modules.json").write_text("[]"); (root/"evidence.json").write_text(json.dumps([ev])); (root/"run_manifest.json").write_text("{}")
                roots.append(root)
            rows, links=compare(roots, Path(temp)/"out")
            self.assertEqual(rows[0]["synteny_status"], "STRONGLY_CONSERVED_CONTEXT")
            self.assertEqual(links[0]["reference_protein_id"], "P1")
            self.assertIn(links[0]["orthology_method"], {"PHROG_ANCHOR", "EXACT_SEQUENCE_MATCH"})
            self.assertIn(links[0]["orthology_method"], {"PHROG_ANCHOR", "MMSEQS2", "EXACT_SEQUENCE_MATCH"})

    def test_compare_exact_fallback_is_not_labeled_mmseqs(self):
        with tempfile.TemporaryDirectory() as temp:
            roots=[]
            for suffix in ("a", "b"):
                root=Path(temp)/suffix; root.mkdir(); seq="M"*10
                (root/"proteins.faa").write_text(f">P1\n{seq}\n")
                (root/"genes.gff3").write_text("##gff-version 3\n")
                (root/"functional_classification.json").write_text(json.dumps([{"protein_id":"P1","functional_state":"UNRESOLVED","proposed_function":None,"functional_category":None,"conservation_status":"NOT_ESTABLISHED"}]))
                (root/"genomic_context.json").write_text(json.dumps([{"protein_id":"P1","gene_order_index":0,"start":1,"end":30,"strand":"+","functional_state":"UNRESOLVED","proposed_function":None,"functional_category":None,"conservation_status":"NOT_ESTABLISHED","functional_module":None}]))
                (root/"modules.json").write_text("[]"); (root/"evidence.json").write_text(json.dumps([{"protein_id":"P1","genome_id":"g","start":1,"end":30,"strand":"+","sequence":seq,"evidence":[]}]))
                (root/"run_manifest.json").write_text("{}")
                roots.append(root)
            _, links=compare(roots, Path(temp)/"out", mmseqs="/unused/mmseqs")
            self.assertEqual(links[0]["orthology_method"], "EXACT_SEQUENCE_MATCH")
            self.assertFalse(links[0]["provenance"]["mmseqs_executed"])

    @pytest.mark.integration
    def test_compare_real_mmseqs_provenance_for_non_phrog_match(self):
        mmseqs = shutil.which("mmseqs")
        if not mmseqs:
            self.skipTest("MMseqs2 unavailable")
        with tempfile.TemporaryDirectory() as temp:
            roots=[]
            for suffix in ("a", "b"):
                root=Path(temp)/suffix; root.mkdir(); seq="MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAN" * 3
                (root/"proteins.faa").write_text(f">P1\n{seq}\n")
                (root/"genes.gff3").write_text("##gff-version 3\n")
                common={"protein_id":"P1","functional_state":"UNRESOLVED","proposed_function":None,"functional_category":None,"conservation_status":"NOT_ESTABLISHED"}
                ctx={**common,"gene_order_index":0,"start":1,"end":90,"strand":"+","functional_module":None}
                ev={"protein_id":"P1","genome_id":"g","start":1,"end":90,"strand":"+","sequence":seq,"evidence":[]}
                (root/"functional_classification.json").write_text(json.dumps([common])); (root/"genomic_context.json").write_text(json.dumps([ctx])); (root/"modules.json").write_text("[]"); (root/"evidence.json").write_text(json.dumps([ev])); (root/"run_manifest.json").write_text("{}")
                roots.append(root)
            _, links=compare(roots, Path(temp)/"out", mmseqs=mmseqs)
            self.assertEqual(links[0]["orthology_method"], "MMSEQS2")
            self.assertTrue(links[0]["provenance"]["mmseqs_executed"])
            self.assertIn("mmseqs_version", links[0]["provenance"])
            self.assertIn("commands", links[0]["provenance"])
    def _context_protein(self, protein_id, start, sequence="M" * 30):
        return Protein("g", protein_id, start, start + len(sequence) * 3 - 1, "+", "ATG" * len(sequence), sequence, "test")

    def _context_classification(self, protein_id, function=None, category=None, state="UNRESOLVED", conservation="NOT_ESTABLISHED"):
        return {"protein_id": protein_id, "functional_state": state, "proposed_function": function, "functional_category": category, "conservation_status": conservation}

    def test_context_coherent_unknown_module_and_unknown_preserved(self):
        proteins = [self._context_protein("A", 1), self._context_protein("B", 100), self._context_protein("C", 199)]
        classifications = [self._context_classification("A", "portal protein", state="PROBABLE_FUNCTION"), self._context_classification("B", state="CONSERVED_UNKNOWN", conservation="STRONGLY_CONSERVED"), self._context_classification("C", "major capsid protein", state="PROBABLE_FUNCTION")]
        records, modules = build_context(proteins, classifications)
        target = next(r for r in records if r["protein_id"] == "B")
        self.assertEqual(target["functional_module"], "head_and_packaging")
        self.assertTrue(target["conserved_unknown"])
        self.assertIsNone(target["proposed_function"])
        self.assertEqual(len(modules), 1)

    def test_context_mixed_neighbors_and_edges(self):
        proteins = [self._context_protein("A", 1), self._context_protein("B", 100), self._context_protein("C", 199), self._context_protein("D", 298)]
        classifications = [self._context_classification("A", "portal protein", state="PROBABLE_FUNCTION"), self._context_classification("B"), self._context_classification("C", "integrase", state="PROBABLE_FUNCTION"), self._context_classification("D", "tail protein", state="PROBABLE_FUNCTION")]
        records, modules = build_context(proteins, classifications)
        self.assertIsNone(next(r for r in records if r["protein_id"] == "B")["functional_module"])
        self.assertTrue(any(r["ambiguous_module_boundary"] for r in records))
        self.assertIsNone(records[0]["upstream_protein_id"])
        self.assertIsNone(records[-1]["downstream_protein_id"])

    def test_context_output_deterministic_and_does_not_mutate_fusion(self):
        proteins = [self._context_protein("A", 1), self._context_protein("B", 100)]
        classifications = [self._context_classification("A", "portal protein", state="PROBABLE_FUNCTION"), self._context_classification("B", state="CONSERVED_UNKNOWN", conservation="CONSERVED")]
        snapshot = json.loads(json.dumps(classifications, sort_keys=True))
        records, modules = build_context(proteins, classifications)
        self.assertEqual(classifications, snapshot)
        self.assertEqual((records, modules), build_context(proteins, classifications))
        with tempfile.TemporaryDirectory() as temp:
            write_context(temp, records, modules)
            self.assertTrue((Path(temp) / "genomic_context.tsv").is_file())
            self.assertTrue((Path(temp) / "genomic_context.json").is_file())
            self.assertTrue((Path(temp) / "modules.tsv").is_file())
            self.assertTrue((Path(temp) / "modules.json").is_file())
    def _fusion_protein(self, evidence=()):
        return Protein("g", "P", 1, 30, "+", "ATG" * 10, "M" * 10, "test", evidence=list(evidence))

    def _fusion_e(self, source, desc, strength="STRONG", supports=True, modality="test", category=None, identifier="x"):
        metrics = {"functional_category": category} if category else {}
        return Evidence(modality, "statement", EvidenceLevel.COMPUTATIONAL, source, "1", status="REAL", supports=supports, identifier=identifier, evidence_strength=strength, description=desc, metrics=metrics)

    def test_fusion_classification_rules_and_determinism(self):
        unknown = classify_protein(self._fusion_protein([self._fusion_e("PHROGs", "unknown function", category="unknown")]))
        self.assertEqual(unknown["functional_state"], "CONSERVED_UNKNOWN")
        probable = classify_protein(self._fusion_protein([self._fusion_e("PHROGs", "RNA polymerase")]))
        self.assertEqual(probable["functional_state"], "PROBABLE_FUNCTION")
        orthology_pair = classify_protein(self._fusion_protein([self._fusion_e("PHROGs", "toxin", identifier="p"), self._fusion_e("VOGDB", "toxin", identifier="v")]))
        self.assertEqual(orthology_pair["functional_state"], "PROBABLE_FUNCTION")
        self.assertEqual(orthology_pair["supporting_source_count"], 2)
        self.assertEqual(orthology_pair["supporting_modality_count"], 2)
        known = classify_protein(self._fusion_protein([self._fusion_e("Swiss-Prot", "major capsid protein"), self._fusion_e("PHROGs", "Major capsid protein {ECO:0001}", identifier="y")]))
        self.assertEqual(known["functional_state"], "KNOWN_FUNCTION")
        self.assertEqual(known["display_classification"], "Specific function strongly supported")
        self.assertIn("PHROGs:y", known["best_evidence"])
        self.assertEqual((known["evidence_tier"],known["evidence_tier_label"]),(2,"Curated and independently corroborated function"))
        broad = classify_protein(self._fusion_protein([self._fusion_e("PHROGs", "unknown function", category="head and packaging")]))
        self.assertEqual(broad["functional_state"], "FUNCTIONAL_CLASS_ONLY")
        unresolved = classify_protein(self._fusion_protein())
        self.assertEqual(unresolved["functional_state"], "UNRESOLVED")
        rejected = classify_protein(self._fusion_protein([self._fusion_e("Swiss-Prot", "DNA polymerase", supports=False)]))
        self.assertEqual(rejected["functional_state"], "UNRESOLVED")
        conflict = classify_protein(self._fusion_protein([self._fusion_e("Swiss-Prot", "integrase"), self._fusion_e("PHROGs", "major capsid protein", identifier="y")]))
        self.assertEqual(conflict["functional_state"], "CONFLICTING_EVIDENCE")
        same_source = classify_protein(self._fusion_protein([self._fusion_e("PHROGs", "integrase", identifier="a"), self._fusion_e("PHROGs", "integrase", identifier="b")]))
        self.assertEqual(same_source["supporting_source_count"], 1)
        self.assertEqual(same_source["supporting_modality_count"], 1)
        rejected_alt = classify_protein(self._fusion_protein([self._fusion_e("Swiss-Prot", "integrase"), self._fusion_e("PHROGs", "major capsid protein", supports=False, identifier="rejected")]))
        self.assertNotEqual(rejected_alt["functional_state"], "CONFLICTING_EVIDENCE")
        hypothetical = classify_protein(self._fusion_protein([self._fusion_e("VOGDB", "REFSEQ hypothetical protein")]))
        self.assertEqual(hypothetical["functional_state"], "CONSERVED_UNKNOWN")
        unavailable = classify_protein(self._fusion_protein())
        self.assertEqual(unavailable["functional_state"], "UNRESOLVED")
        original = self._fusion_e("PHROGs", "RNA polymerase")
        before = original.to_dict() if hasattr(original, "to_dict") else original.__dict__.copy()
        classify_protein(self._fusion_protein([original]))
        self.assertEqual(before, original.__dict__)
        self.assertEqual(normalize_function("Major capsid protein {ECO:0001}"), "major capsid protein")
        self.assertEqual(unknown, classify_protein(self._fusion_protein([self._fusion_e("PHROGs", "unknown function", category="unknown")])))

    def test_fusion_moderate_informative_reason_uses_selected_evidence(self):
        evidence = self._fusion_e("PHROGs", "RNA polymerase", strength="MODERATE", identifier="moderate-1")
        result = classify_protein(self._fusion_protein([evidence]))
        self.assertEqual(result["functional_state"], "PROBABLE_FUNCTION")
        self.assertIn("PHROGs", result["confidence_reasons"][0])
        self.assertIn("MODERATE", result["confidence_reasons"][0])

    def test_fusion_reduces_mitochondrial_bcs1_to_family_level(self):
        result = classify_protein(self._fusion_protein([
            self._fusion_e("VOGDB", "sp|P32839|BCS1_YEAST Mitochondrial chaperone BCS1"),
            self._fusion_e("Pfam", "ATPase family associated with various cellular activities (AAA)"),
        ]))
        self.assertEqual(result["proposed_function"], "bcs1-like aaa-family atpase")
        self.assertNotIn("mitochondrial", result["proposed_function"])
        self.assertTrue(any("organelle-specific" in flag for flag in result["ambiguity_flags"]))

    def test_fusion_removes_mimivirus_locus_from_band7_product(self):
        result = classify_protein(self._fusion_protein([
            self._fusion_e("VOGDB", "sp|Q5UP73|YR614_MIMIV Putative band 7 family protein R614"),
            self._fusion_e("Pfam", "SPFH domain / Band 7 family"),
        ]))
        self.assertEqual(result["proposed_function"], "band 7/spfh family protein")
        self.assertNotIn("r614", result["proposed_function"])

    def test_fusion_reduces_opg_locus_to_independent_kelch_architecture(self):
        result = classify_protein(self._fusion_protein([
            self._fusion_e("VOGDB", "sp|A0A7H0DN20|PG047_MONPV Immune evasion protein OPG047"),
            self._fusion_e("Pfam", "Kelch-repeats beta-propeller domain"),
        ]))
        self.assertEqual(result["proposed_function"], "kelch-repeat beta-propeller protein")
        self.assertNotIn("immune evasion", result["proposed_function"])

    def test_fusion_suppresses_uncorroborated_opg_function(self):
        result = classify_protein(self._fusion_protein([
            self._fusion_e("VOGDB", "Immune evasion protein OPG047"),
        ]))
        self.assertIsNone(result["proposed_function"])
        self.assertEqual(result["functional_state"], "CONSERVED_UNKNOWN")

    def test_fusion_never_promotes_domain_only_evidence_to_probable_function(self):
        result = classify_protein(self._fusion_protein([
            self._fusion_e("Pfam", "Anthrax toxin lethal factor, middle domain"),
        ]))
        self.assertEqual(result["functional_state"], "FUNCTIONAL_CLASS_ONLY")
        self.assertIsNone(result["proposed_function"])
        self.assertEqual(result["display_product"], "hypothetical protein")
        self.assertEqual(result["domain_summary"], "Detected domain evidence: domain x")

    def test_fusion_suppresses_taxon_specific_products_but_preserves_raw_evidence(self):
        for unsafe in ("Interferon gamma", "Apolipoprotein CIII", "Centromere kinetochore component CENP-T", "Male sterility protein"):
            evidence = self._fusion_e("VOGDB", unsafe)
            result = classify_protein(self._fusion_protein([evidence]))
            self.assertNotEqual(result["functional_state"], "PROBABLE_FUNCTION")
            self.assertIsNone(result["proposed_function"])
            self.assertEqual(evidence.description, unsafe)

    def test_probable_function_always_has_a_defensible_product(self):
        result = classify_protein(self._fusion_protein([
            self._fusion_e("PHROGs", "portal protein"),
            self._fusion_e("VOGDB", "major capsid protein"),
        ]))
        self.assertNotEqual(result["functional_state"], "PROBABLE_FUNCTION")
        self.assertIsNone(result["proposed_function"])

    def test_fusion_outputs_cover_integration_fixture(self):
        proteins = [self._fusion_protein() for _ in range(99)]
        for i, protein in enumerate(proteins, 1): protein.protein_id = f"P{i:03d}"
        with tempfile.TemporaryDirectory() as temp:
            write_classification(temp, proteins)
            self.assertEqual(len(json.loads((Path(temp) / "functional_classification.json").read_text())), 99)
            self.assertEqual(len((Path(temp) / "functional_classification.tsv").read_text().splitlines()), 100)
            self.assertEqual([x["protein_id"] for x in classify_proteins(proteins)], [x["protein_id"] for x in classify_proteins(proteins)])

    def test_fusion_conserved_unknown_from_phrogs_and_vogdb(self):
        result = classify_protein(self._fusion_protein([self._fusion_e("PHROGs", "unknown function"), self._fusion_e("VOGDB", "hypothetical protein")]))
        self.assertEqual(result["functional_state"], "CONSERVED_UNKNOWN")

    def test_fusion_phrogs_and_vogdb_same_function_is_probable(self):
        result = classify_protein(self._fusion_protein([self._fusion_e("PHROGs", "portal protein"), self._fusion_e("VOGDB", "portal protein")]))
        self.assertEqual(result["functional_state"], "PROBABLE_FUNCTION")
        self.assertNotEqual(result["functional_state"], "KNOWN_FUNCTION")

    def test_fusion_category_only(self):
        result = classify_protein(self._fusion_protein([self._fusion_e("PHROGs", "unknown function", category="head and packaging")]))
        self.assertEqual(result["functional_state"], "FUNCTIONAL_CLASS_ONLY")

    def test_fusion_supporting_source_and_modality_counts(self):
        result = classify_protein(self._fusion_protein([self._fusion_e("PHROGs", "portal protein"), self._fusion_e("PHROGs", "portal protein", identifier="2"), self._fusion_e("VOGDB", "portal protein")]))
        self.assertEqual(result["supporting_source_count"], 2)
        self.assertEqual(result["supporting_modality_count"], 2)

    def test_fusion_rejected_alternative_cannot_conflict(self):
        result = classify_protein(self._fusion_protein([self._fusion_e("Swiss-Prot", "integrase"), self._fusion_e("PHROGs", "major capsid protein", supports=False)]))
        self.assertNotEqual(result["functional_state"], "CONFLICTING_EVIDENCE")

    def test_fusion_unavailable_is_not_negative_evidence(self):
        result = classify_protein(self._fusion_protein([self._fusion_e("PHROGs", "portal protein")]))
        self.assertEqual(result["functional_state"], "PROBABLE_FUNCTION")

    def test_fusion_ambiguous_labels_are_not_conflict(self):
        result = classify_protein(self._fusion_protein([self._fusion_e("PHROGs", "DNA-dependent RNA polymerase"), self._fusion_e("VOGDB", "T7 RNA polymerase")]))
        self.assertNotEqual(result["functional_state"], "CONFLICTING_EVIDENCE")
        self.assertTrue(result["ambiguity_flags"])

    def test_fusion_quantitative_phrogs_winner_drives_product_and_best_evidence(self):
        weak = self._fusion_e("PHROGs", "head closure Hc1", identifier="phrog_161")
        weak.metrics.update({"percent_identity": 0.31, "query_coverage": 0.55, "bit_score": 80, "evalue": 1e-18})
        winner = self._fusion_e("PHROGs", "tail completion or Neck1 protein", identifier="phrog_17193")
        winner.metrics.update({"percent_identity": 0.849, "query_coverage": 0.95, "bit_score": 231, "evalue": 2.5e-68})
        result = classify_protein(self._fusion_protein([weak, winner]))
        self.assertEqual(result["proposed_function"], "tail completion or neck1 protein")
        self.assertIn("phrog_17193", result["best_evidence"])
        self.assertEqual(len(result["product_alternatives"]), 2)

    def test_fusion_phrogs_major_head_normalizes_to_major_capsid(self):
        evidence = self._fusion_e("PHROGs", "major head protein", identifier="phrog_247")
        evidence.metrics.update({"percent_identity": 0.546, "query_coverage": 0.9, "bit_score": 397, "evalue": 1e-118})
        result = classify_protein(self._fusion_protein([evidence]))
        self.assertEqual(result["proposed_function"], "major capsid protein")

    def test_fusion_hoc_product_is_separate_from_immunoglobulin_domain_note(self):
        phrog = self._fusion_e("PHROGs", "Hoc-like head decoration", identifier="phrog_2973")
        phrog.metrics.update({"percent_identity": 0.772, "query_coverage": 0.9, "bit_score": 212, "evalue": 3e-63})
        result = classify_protein(self._fusion_protein([phrog, self._fusion_e("Pfam", "immunoglobulin I-set domain")]))
        self.assertEqual(result["proposed_function"], "hoc-like head decoration protein")
        self.assertIn("immunoglobulin", result["domain_summary"])
        self.assertIn("capsid-display", result["biotechnology_relevance"])

    def test_fusion_suppresses_ijeoma_taxon_inappropriate_domain_summaries(self):
        for description in (
            "fungal deubiquitinating enzyme PH domain",
            "apoptosis antagonizing transcription factor",
            "spore coat protein U domain",
            "protein involved in starch initiation 1 C-terminal domain",
            "ZNF598 C2H2 zinc finger domain",
            "queuosine salvage protein",
        ):
            result = classify_protein(self._fusion_protein([self._fusion_e("Pfam", description)]))
            self.assertNotIn(description.lower(), result["display_product"])

    def test_organism_specific_domain_text_is_note_only_and_neutralized(self):
        for description in (
            "coronavirus replicase NSP2, C-terminal",
            "influenza C hemagglutinin stalk",
            "Salmonella virulence plasmid 28.1 kDa A protein",
            "spirochaete C-terminal Tudor-like domain",
        ):
            result = classify_protein(self._fusion_protein([self._fusion_e("Pfam", description)]))
            self.assertEqual(result["display_product"], "hypothetical protein")
            self.assertNotIn(description.lower(), result["domain_summary"].lower())
            self.assertFalse(result["confidence_calibrated"])

    def test_truth_metrics_separate_precision_coverage_and_abstention(self):
        truth = [
            {"start": 1, "end": 90, "strand": "+", "product": "major capsid protein"},
            {"start": 100, "end": 180, "strand": "+", "product": "endolysin"},
        ]
        predicted = [
            {"start": 1, "end": 90, "strand": "+", "product": "major head protein"},
            {"start": 100, "end": 180, "strand": "+", "product": "hypothetical protein"},
        ]
        result = truth_metrics(predicted, truth)
        self.assertEqual(result["correct_product_assertions"], 1)
        self.assertEqual(result["product_abstentions"], 1)
        self.assertEqual(result["asserted_product_precision"], 1.0)
        self.assertEqual(result["functional_coverage"], 0.5)

    def test_ijeoma_diagnostic_domain_and_compatible_function_rules(self):
        polymerase = classify_protein(self._fusion_protein([
            self._fusion_e("Pfam", "DNA polymerase family A")
        ]))
        self.assertEqual(polymerase["proposed_function"], "family-a dna polymerase")
        self.assertEqual(polymerase["confidence"], "MODERATE")

        resolvase_phrog = self._fusion_e("PHROGs", "Holliday junction resolvase", identifier="phrog_5702")
        resolvase_phrog.metrics.update({"percent_identity": 0.648, "query_coverage": 0.8, "bit_score": 120, "evalue": 4e-23})
        endonuclease = self._fusion_e("PHROGs", "endonuclease", identifier="phrog_93")
        endonuclease.metrics.update({"percent_identity": 0.7, "query_coverage": 0.9, "bit_score": 140, "evalue": 1e-30})
        resolvase = classify_protein(self._fusion_protein([
            resolvase_phrog, endonuclease, self._fusion_e("Pfam", "VRR-Nuc domain")
        ]))
        self.assertEqual(resolvase["proposed_function"], "holliday junction resolvase")
        self.assertIn("VRR-Nuc", resolvase["diagnostic_domain_rule"])

    def test_ijeoma_regulatory_and_ssdna_rules(self):
        excisionase = self._fusion_e("PHROGs", "excisionase and transcriptional regulator", identifier="phrog_66", category="integration and excision")
        terminase = self._fusion_e("VOGDB", "terminase small subunit", identifier="VOG04175")
        result = classify_protein(self._fusion_protein([
            excisionase, terminase, self._fusion_e("Pfam", "helix-turn-helix domain")
        ]))
        self.assertEqual(result["proposed_function"], "excisionase and transcriptional regulator")
        self.assertTrue(result["conflict_resolved_by_corroboration"])
        ssdna = classify_protein(self._fusion_protein([
            self._fusion_e("Pfam", "Enterobacter phage Enc34, ssDNA-binding protein")
        ]))
        self.assertEqual(ssdna["proposed_function"], "single-stranded dna-binding protein")
        self.assertEqual(ssdna["confidence"], "LOW")

    def test_ijeoma_unsafe_mixed_domains_are_rewritten_without_inferring_false_orf(self):
        ninh = classify_protein(self._fusion_protein([
            self._fusion_e("Pfam", "phage NinH protein"),
            self._fusion_e("Pfam", "transposase"),
        ]))
        self.assertEqual(ninh["display_product"], "hypothetical protein")
        self.assertIn("ninh-like domain", ninh["domain_summary"])
        short = self._fusion_protein([self._fusion_e("Pfam", "P-type ATPase actuator domain and DISARM protein DrmE, C-terminal domain")])
        short.sequence = "M" * 63
        result = classify_protein(short)
        self.assertEqual(result["display_product"], "hypothetical protein")
        self.assertIn("atpase-related domain", result["domain_summary"])
        self.assertEqual(result["review_flag"], "NONE")

    def test_repeated_phrogs_hnh_profiles_support_hnh_endonuclease(self):
        evidence = []
        for identifier, evalue, coverage, score in (
            ("phrog_99", 6.935e-20, 0.946, 86.0),
            ("phrog_4578", 3.7e-23, 0.518, 82.8),
            ("phrog_3593", 5.6e-12, 0.633, 46.8),
        ):
            hit = self._fusion_e("PHROGs", "HNH endonuclease", identifier=identifier)
            hit.metrics.update({"evalue": evalue, "query_coverage": coverage, "bit_score": score})
            evidence.append(hit)
        evidence.append(self._fusion_e("Pfam", "HNH endonuclease domain"))
        result = classify_protein(self._fusion_protein(evidence))
        self.assertEqual(result["proposed_function"], "hnh endonuclease")
        self.assertEqual(result["confidence"], "MODERATE")

    def test_annotated_phrogs_consensus_outranks_unannotated_profile(self):
        unknown = self._fusion_e("PHROGs", None, identifier="phrog_38319")
        unknown.metrics.update({"evalue": 1.6e-42, "query_coverage": 0.549, "bit_score": 161})
        first = self._fusion_e("PHROGs", "terminase large subunit", identifier="phrog_10250")
        first.metrics.update({"evalue": 1.7e-37, "query_coverage": 0.637, "bit_score": 146})
        second = self._fusion_e("PHROGs", "terminase large subunit", identifier="phrog_2")
        second.metrics.update({"evalue": 6.3e-36, "query_coverage": 0.912, "bit_score": 125})
        result = classify_protein(self._fusion_protein([
            unknown, first, second,
            self._fusion_e("Pfam", "terminase RNaseH-like domain"),
        ]))
        self.assertEqual(result["proposed_function"], "terminase large subunit")
        self.assertEqual(result["confidence"], "MODERATE")

    def test_danladi_chimallin_nuclear_shell_synonyms_converge(self):
        phrog = self._fusion_e("PHROGs", "nuclear shell protein", identifier="phrog_4795")
        phrog.metrics.update({"query_coverage": 0.981, "sequence_identity": 0.398,
                              "bit_score": 461, "evalue": 1.526e-141})
        vog = self._fusion_e("VOGDB", "sp|F8SJT5|CHMA_BPPA3 Chimallin", identifier="VOG27585")
        vog.metrics.update({"query_coverage": 0.997, "bit_score": 796.8, "evalue": 4.3e-239})
        swiss = self._fusion_e("Swiss-Prot", "sp|F8SJT5|CHMA_BPPA3 Chimallin",
                               strength="MODERATE", identifier="F8SJT5")
        swiss.metrics.update({
            "reviewed": True, "organism": "Pseudomonas phage PA3.",
            "percent_identity": 65.9, "query_coverage": 575 / 621,
            "subject_coverage": 575 / 602, "evalue": 3.92e-274,
            "bit_score": 765,
        })
        result = classify_protein(self._fusion_protein([phrog, vog, swiss]))
        self.assertEqual(result["proposed_function"], "chimallin")
        self.assertEqual(result["product_supporting_source_count"], 3)
        self.assertTrue(result["selected_by_curated_phage_anchor"])
        self.assertEqual(result["review_flag"], "NONE")

    def test_danladi_gp144_curated_phage_anchor_beats_partial_tail_label_but_keeps_review(self):
        tail = self._fusion_e("PHROGs", "tail fiber protein", identifier="phrog_3201")
        tail.metrics.update({"sequence_identity": 0.494, "query_coverage": 0.718,
                             "bit_score": 179, "evalue": 1.961e-50})
        vog = self._fusion_e("VOGDB", "sp|Q8SD18|ENLYS_BPDPK Endolysin gp144",
                             identifier="VOG15058")
        vog.metrics.update({"query_coverage": 200 / 259, "bit_score": 228.1,
                            "evalue": 2.6e-67})
        swiss = self._fusion_e("Swiss-Prot",
                               "Endolysin gp144 {ECO:0000303|PubMed:28461978}",
                               strength="MODERATE", identifier="Q8SD18")
        swiss.metrics.update({
            "reviewed": True, "organism": "Pseudomonas phage phiKZ.",
            "percent_identity": 56.8, "query_coverage": 257 / 259,
            "subject_coverage": 257 / 260, "evalue": 2.54e-88,
            "bit_score": 265,
        })
        result = classify_protein(self._fusion_protein([tail, vog, swiss]))
        self.assertEqual(result["proposed_function"], "endolysin")
        self.assertTrue(result["selected_by_curated_phage_anchor"])
        self.assertIn("Swiss-Prot:Q8SD18", result["best_evidence"])
        self.assertEqual(result["product_supporting_source_count"], 2)
        self.assertEqual(result["review_flag"], "REVIEW_REQUIRED")

    def test_danladi_gp181_curated_structural_hydrolase_not_silently_overruled_by_tail_fiber(self):
        tail = self._fusion_e("PHROGs", "tail fiber protein", identifier="phrog_3201")
        tail.metrics.update({"sequence_identity": 0.60, "query_coverage": 0.95,
                             "bit_score": 600, "evalue": 1e-150})
        vog = self._fusion_e("VOGDB",
                             "sp|Q8SCY1|EXLYS_BPDPK Peptidoglycan hydrolase gp181",
                             identifier="VOG42004")
        vog.metrics.update({"query_coverage": 1.0, "bit_score": 2827.4, "evalue": 0.0})
        swiss = self._fusion_e("Swiss-Prot",
                               "Peptidoglycan hydrolase gp181 {ECO:0000305}",
                               strength="MODERATE", identifier="Q8SCY1")
        swiss.metrics.update({
            "reviewed": True, "organism": "Pseudomonas phage phiKZ.",
            "percent_identity": 38.5,
            "alignment_length": 2171,
            "query_length": 2256,
            "subject_length": 2237,
            "query_coverage": 0.9321808510638298,
            "subject_coverage": 0.9400983459991059,
            "evalue": 0.0,
            "bit_score": 1316.0,
        })
        result = classify_protein(self._fusion_protein([tail, vog, swiss]))
        self.assertEqual(result["proposed_function"], "structural peptidoglycan hydrolase")
        self.assertTrue(result["selected_by_curated_phage_anchor"])
        self.assertIn("Swiss-Prot:Q8SCY1", result["best_evidence"])
        self.assertEqual(result["product_supporting_source_count"], 2)
        self.assertEqual(result["review_flag"], "REVIEW_REQUIRED")

    def test_head_and_capsid_maturation_protease_are_semantic_synonyms(self):
        phrog = self._fusion_e("PHROGs", "head maturation protease", identifier="phrog_3014")
        swiss = self._fusion_e("Swiss-Prot", "Capsid maturation protease {ECO:0000305}",
                               strength="MODERATE", identifier="Q8SCY7")
        swiss.metrics.update({
            "reviewed": True, "organism": "Pseudomonas phage phiKZ.",
            "percent_identity": 60.0, "query_coverage": 0.95,
            "subject_coverage": 0.95, "evalue": 1e-100,
        })
        result = classify_protein(self._fusion_protein([phrog, swiss]))
        self.assertEqual(result["proposed_function"], "capsid maturation protease")
        self.assertNotIn(
            "distinct strong labels could not be confidently established as biologically incompatible",
            result["ambiguity_flags"],
        )
        self.assertEqual(result["review_flag"], "NONE")

    def test_lower_identity_curated_phage_anchor_requires_near_full_length_coverage(self):
        swiss = self._fusion_e("Swiss-Prot", "Peptidoglycan hydrolase gp181",
                               strength="MODERATE", identifier="Q8SCY1")
        swiss.metrics.update({
            "reviewed": True, "organism": "Pseudomonas phage phiKZ.",
            "percent_identity": 38.5,
            "query_coverage": 0.70,
            "subject_coverage": 0.70,
            "evalue": 0.0,
            "bit_score": 600,
        })
        tail = self._fusion_e("PHROGs", "tail fiber protein", identifier="phrog_3201")
        tail.metrics.update({
            "percent_identity": 0.60, "query_coverage": 0.95,
            "bit_score": 700, "evalue": 1e-180,
        })
        result = classify_protein(self._fusion_protein([tail, swiss]))
        self.assertFalse(result["selected_by_curated_phage_anchor"])
        self.assertEqual(result["proposed_function"], "tail fiber protein")
        self.assertEqual(result["review_flag"], "NONE")

    def test_product_evidence_tier_counts_only_sources_supporting_selected_product(self):
        phrog = self._fusion_e("PHROGs", "RNA polymerase", identifier="p")
        pfam = self._fusion_e("Pfam", "RNA polymerase beta subunit domain", identifier="PFTEST")
        result = classify_protein(self._fusion_protein([phrog, pfam]))
        self.assertEqual(result["proposed_function"], "rna polymerase")
        self.assertEqual(result["supporting_source_count"], 2)
        self.assertEqual(result["product_supporting_source_count"], 1)
        self.assertEqual(result["evidence_tier"], 4)

    def test_context_refinement_never_invents_minor_tail_without_candidate(self):
        proteins=[]
        for index,(description,start) in enumerate((("tail protein",1),("unknown function",301),("minor tail protein",601)),1):
            evidence=self._fusion_e("PHROGs",description,identifier=f"p{index}")
            protein=self._fusion_protein([evidence]); protein.protein_id=f"P{index}"; protein.start=start; protein.end=start+299
            proteins.append(protein)
        results=classify_proteins(proteins)
        self.assertEqual(results[1]["proposed_function"],"putative tail-associated protein")
        self.assertEqual(results[1]["confidence"],"LOW")
        self.assertEqual(results[1]["context_support"]["module"],"tail")

    def test_context_refinement_uses_minor_tail_only_when_candidate_exists(self):
        left=self._fusion_protein([self._fusion_e("PHROGs","tail protein")]); left.protein_id="P1"; left.start=1; left.end=300
        target=self._fusion_protein([
            self._fusion_e("PHROGs","minor tail protein",identifier="minor"),
            self._fusion_e("VOGDB","portal protein",identifier="portal"),
        ]); target.protein_id="P2"; target.start=301; target.end=600
        right=self._fusion_protein([self._fusion_e("PHROGs","tail terminator")]); right.protein_id="P3"; right.start=601; right.end=900
        result=classify_proteins([left,target,right])[1]
        self.assertEqual(result["proposed_function"],"putative minor tail protein")

    def test_specific_head_morphogenesis_term_beats_broad_head_term(self):
        broad=self._fusion_e("PHROGs","head protein",identifier="broad")
        broad.metrics.update({"percent_identity":0.6,"query_coverage":0.9,"bit_score":220,"evalue":1e-60})
        specific=self._fusion_e("PHROGs","head morphogenesis",identifier="specific")
        specific.metrics.update({"percent_identity":0.42,"query_coverage":0.8,"bit_score":180,"evalue":1e-50})
        result=classify_protein(self._fusion_protein([broad,specific]))
        self.assertEqual(result["proposed_function"],"head morphogenesis protein")

    def test_benchmark_product_semantics_do_not_reward_domain_verbosity(self):
        self.assertEqual(classify_product_relation("major capsid protein","major head protein"),"EXACT_PRODUCT_AGREEMENT")
        self.assertEqual(classify_product_relation("conserved phage protein of unknown function","hypothetical protein"),"BOTH_FUNCTION_UNKNOWN")
        self.assertEqual(classify_product_relation("hypothetical protein containing ASCH domain","hypothetical protein"),"PHAGEMINE_DOMAIN_INFORMATION_ONLY")
        self.assertEqual(classify_product_relation("single-stranded DNA-binding protein","Gp2.5-like ssDNA binding protein and ssDNA annealing protein"),"EQUIVALENT_FUNCTION")

    def test_fusion_resume_output_layer_regenerates_both_files(self):
        proteins = [self._fusion_protein()]
        with tempfile.TemporaryDirectory() as temp:
            write_classification(temp, proteins)
            self.assertTrue((Path(temp) / "functional_classification.json").is_file())
            self.assertTrue((Path(temp) / "functional_classification.tsv").is_file())
    def test_resume_progress_uses_resume_stage_ledger_without_duplicate_labels(self):
        stream = StringIO()
        progress = ProgressReporter(stream=stream)
        progress.STAGES = RESUME_STAGES
        progress.start("PHROGs")
        progress.finish("accepted hits")
        output = stream.getvalue()
        self.assertIn("[PhageMine 0/11] RUNNING: PHROGs", output)
        self.assertIn("[PhageMine 1/11] DONE: PHROGs", output)
        self.assertNotIn("RUNNING: RUNNING:", output)

    def test_resume_source_loader_reuses_proteins_and_evidence_verbatim(self):
        source = RESUME_SOURCE_FIXTURE
        manifest, representation, proteins, sequence, _ = _load_source(source)
        persisted = json.loads((source / "evidence.json").read_text())
        self.assertEqual(len(proteins), 2)
        self.assertEqual([protein.sequence for protein in proteins], [record["sequence"] for record in persisted])
        self.assertEqual(proteins[0].protein_id, "PM_000001")
        self.assertEqual(proteins[0].evidence[0].statement, persisted[0]["evidence"][0]["statement"])
        self.assertEqual(proteins[0].evidence[0].provenance, persisted[0]["evidence"][0]["provenance"])
        self.assertEqual(proteins[0].evidence[0].provenance["adapter"], "VOGHMMAdapter")
        self.assertEqual(manifest["evidence_adapters"][0]["provenance"]["threshold_mode"], "GA")
        self.assertEqual(representation.analysis_sequence, sequence)

    def test_resume_source_loader_rejects_protein_sequence_mismatch(self):
        with tempfile.TemporaryDirectory() as temp:
            source = Path(temp) / "source"
            shutil.copytree(RESUME_SOURCE_FIXTURE, source)
            fasta = source / "proteins.faa"
            text = fasta.read_text()
            fasta.write_text(text.replace("MISQDKFEYEISAMK", "MSSQDKFEYEISAMK", 1))
            with self.assertRaisesRegex(ValueError, "protein sequence/length mismatch"):
                _load_source(source)

    def test_resume_source_loader_rejects_missing_artifact(self):
        with tempfile.TemporaryDirectory() as temp:
            source = Path(temp) / "source"
            shutil.copytree(RESUME_SOURCE_FIXTURE, source)
            (source / "evidence.json").unlink()
            with self.assertRaisesRegex(ValueError, "missing source artifacts"):
                _load_source(source)

    def test_reclassify_reuses_evidence_without_database_adapters(self):
        with tempfile.TemporaryDirectory() as temp, \
             patch("phagemine.phrogs.PHROGSMMseqsAdapter.analyze") as phrogs, \
             patch("phagemine.swissprot.SwissProtEvidenceAdapter.analyze") as swissprot, \
             patch("phagemine.pfam.PfamHMMAdapter.analyze") as pfam, \
             patch("phagemine.vog.VOGHMMAdapter.analyze") as vog:
            output = Path(temp) / "reclassified"
            count = reclassify(RESUME_SOURCE_FIXTURE, output)
            self.assertEqual(count, 2)
            self.assertTrue((output / "annotation.tsv").is_file())
            self.assertTrue((output / "functional_classification.json").is_file())
            self.assertTrue((output / "annotated_proteins.faa").is_file())
            self.assertTrue((output / "genbank_submission").is_dir())
            self.assertFalse(phrogs.called or swissprot.called or pfam.called or vog.called)
            manifest = json.loads((output / "run_manifest.json").read_text())
            self.assertEqual(manifest["command"], "reclassify")
            self.assertEqual(manifest["reclassification"]["database_searches_run"], [])

    def test_reclassify_cli_is_visible_and_requires_new_output(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "reclassified"
            self.assertEqual(main(["reclassify", str(RESUME_SOURCE_FIXTURE), "--output", str(output)]), 0)
            self.assertTrue((output / "annotation.tsv").is_file())
    def test_progress_stage_transitions_and_non_tty_output(self):
        stream = StringIO()
        progress = ProgressReporter(stream=stream)
        progress.start("Pfam")
        progress.finish("2 accepted hits / 1 proteins")
        text = stream.getvalue()
        self.assertIn("RUNNING: Pfam", text)
        self.assertIn("DONE: Pfam", text)
        self.assertIn("2 accepted hits", text)

    def test_unavailable_evidence_progress_is_not_reported_as_zero_hits(self):
        result = EvidenceAdapterResult(
            "PfamHMMAdapter",
            "UNAVAILABLE",
            message="Pfam/HMMER unavailable; no domain evidence was fabricated.",
        )
        summary = _evidence_progress_summary(result, "Pfam/HMMER unavailable")
        self.assertIn("UNAVAILABLE", summary)
        self.assertNotIn("0 accepted hits", summary)

    def test_completed_zero_hit_search_remains_distinct_from_unavailable(self):
        result = EvidenceAdapterResult("PfamHMMAdapter", "SUCCESS_NO_HIT")
        summary = _evidence_progress_summary(result, "Pfam/HMMER unavailable")
        self.assertEqual(summary, "0 accepted hits / 0 proteins")

    def test_single_run_consumes_registered_pmfdb_and_inphared_resources(self):
        from phagemine.resources import ResourceType
        def registered(_manager, kind):
            if kind == ResourceType.PMFDB:
                return {"path": "/db/pmfdb", "resource_type": "PMFDB"}
            if kind == ResourceType.INPHARED_GENOMES:
                return {"path": "/db/inphared", "resource_type": "INPHARED_GENOMES", "provenance": {}}
            return None
        comparative = {"pmfdb": {"status": "COMPLETE"}, "inphared": {"status": "SKIPPED", "matches": []}}
        resolved = {"path": "/db/reference_phage_genomes.fna", "version": "2026-04-07",
                    "runtime_paths": {"reference_fasta": "/db/reference_phage_genomes.fna",
                                      "mash_index": "/db/inphared.msh", "metadata": "/db/genome_metadata.tsv",
                                      "manifest": "/db/genome_manifest.json"}}
        comparison = {"status": "COMPLETE", "matches": [], "unique_reference_summaries": []}
        with tempfile.TemporaryDirectory() as temp, \
             patch("phagemine.resources.EvidenceResourceManager.find", new=registered), \
             patch("phagemine.pipeline.resolve_validated_inphared", return_value=(resolved, "ready")), \
             patch("phagemine.inphared.compare_genomes", return_value=comparison) as compare, \
             patch("phagemine.discovery.build_discovery_outputs", return_value=comparative) as build:
            output = Path(temp) / "run"
            run(ROOT / "examples/demo_phage.fasta", output, command="run", predictor=DemoORFPredictor())
            self.assertTrue(build.called)
            kwargs = build.call_args.kwargs
            self.assertEqual(kwargs["pmfdb"], "/db/pmfdb")
            self.assertIsNone(kwargs["inphared"])
            self.assertEqual(compare.call_count, 1)
            manifest = json.loads((output / "run_manifest.json").read_text())
            self.assertEqual(manifest["comparative_analysis"]["pmfdb"]["status"], "COMPLETE")
            self.assertEqual(manifest["comparative_analysis"]["inphared"]["status"], "COMPLETE")

    def test_annotate_consumes_registered_inphared_without_discovery_workflow(self):
        from phagemine.resources import ResourceType

        def registered(_manager, kind):
            if kind == ResourceType.INPHARED_GENOMES:
                return {
                    "path": "/db/reference_phage_genomes.fna",
                    "version": "2026-04-07",
                    "resource_type": "INPHARED_GENOMES",
                    "provenance": {
                        "mash_index_path": "/db/inphared.msh",
                        "metadata_path": "/db/genome_metadata.tsv",
                        "reference_manifest_path": "/db/genome_manifest.json",
                    },
                }
            return None

        comparison = {
            "status": "COMPLETE",
            "matches": [],
            "unique_reference_summaries": [],
        }

        resolved = {"path": "/db/reference_phage_genomes.fna", "version": "2026-04-07",
                    "runtime_paths": {"reference_fasta": "/db/reference_phage_genomes.fna",
                                      "mash_index": "/db/inphared.msh", "metadata": "/db/genome_metadata.tsv",
                                      "manifest": "/db/genome_manifest.json"}}
        with tempfile.TemporaryDirectory() as temp, \
             patch("phagemine.pipeline.resolve_validated_inphared", return_value=(resolved, "ready")), \
             patch("phagemine.inphared.compare_genomes", return_value=comparison) as compare, \
             patch("phagemine.discovery.build_discovery_outputs") as discovery:
            output = Path(temp) / "annotate"

            run(
                ROOT / "examples/demo_phage.fasta",
                output,
                command="annotate",
                predictor=DemoORFPredictor(),
            )

            self.assertTrue(compare.called)
            self.assertFalse(discovery.called)

            kwargs = compare.call_args.kwargs
            self.assertEqual(kwargs["mash_index"], "/db/inphared.msh")
            self.assertEqual(kwargs["metadata"], "/db/genome_metadata.tsv")
            self.assertEqual(
                kwargs["reference_fasta"],
                "/db/reference_phage_genomes.fna",
            )
            self.assertEqual(
                kwargs["output"],
                output / "comparative",
            )

            manifest = json.loads((output / "run_manifest.json").read_text())
            self.assertEqual(
                manifest["comparative_analysis"]["inphared"]["status"],
                "COMPLETE",
            )

    def test_annotate_skips_incomplete_registered_inphared_without_crashing(self):
        from phagemine.resources import EvidenceResourceManager, ResourceType
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            reference = root / "reference_phage_genomes.fna"
            reference.write_text(">REF1\nACGT\n")
            registry = root / "resources.json"
            manager = EvidenceResourceManager(registry)
            manager.register("incomplete", ResourceType.INPHARED_GENOMES, reference)
            output = root / "annotate"
            progress_stream = StringIO()
            with patch.dict(os.environ, {"PHAGEMINE_REGISTRY_PATH": str(registry)}), \
                 patch("phagemine.resources._tool_available", return_value=True):
                run(ROOT / "examples/demo_phage.fasta", output, command="annotate",
                    predictor=DemoORFPredictor(), progress=ProgressReporter(stream=progress_stream))
            manifest = json.loads((output / "run_manifest.json").read_text())
            self.assertEqual(manifest["comparative_analysis"]["inphared"]["status"], "SKIPPED")
            self.assertIn("inphared.msh", manifest["comparative_analysis"]["inphared"]["reason"])
            self.assertIn("SKIPPED: INPHARED genome comparison", progress_stream.getvalue())

    def test_annotate_skips_inphared_when_mash_is_missing(self):
        reason = "INPHARED resource is not operationally READY: required executable unavailable: mash"
        with tempfile.TemporaryDirectory() as temp, \
             patch("phagemine.pipeline.resolve_validated_inphared", return_value=(None, reason)), \
             patch("phagemine.inphared.compare_genomes") as compare:
            output = Path(temp) / "annotate"
            run(ROOT / "examples/demo_phage.fasta", output, command="annotate", predictor=DemoORFPredictor())
            manifest = json.loads((output / "run_manifest.json").read_text())
            self.assertEqual(manifest["comparative_analysis"]["inphared"]["status"], "SKIPPED")
            self.assertIn("required executable unavailable: mash",
                          manifest["comparative_analysis"]["inphared"]["reason"])
            self.assertFalse(compare.called)

    def test_run_with_only_inphared_does_not_invoke_pmf_discovery(self):
        resolved = {"path": "/db/reference_phage_genomes.fna", "version": "2026-04-07",
                    "runtime_paths": {"reference_fasta": "/db/reference_phage_genomes.fna",
                                      "mash_index": "/db/inphared.msh", "metadata": "/db/genome_metadata.tsv",
                                      "manifest": "/db/genome_manifest.json"}}
        comparison = {"status": "COMPLETE", "matches": [], "unique_reference_summaries": []}
        with tempfile.TemporaryDirectory() as temp, \
             patch("phagemine.pipeline.resolve_validated_inphared", return_value=(resolved, "ready")), \
             patch("phagemine.inphared.compare_genomes", return_value=comparison) as compare, \
             patch("phagemine.discovery.build_discovery_outputs") as discovery:
            output = Path(temp) / "run"
            run(ROOT / "examples/demo_phage.fasta", output, command="run", predictor=DemoORFPredictor())
            self.assertEqual(compare.call_count, 1)
            self.assertFalse(discovery.called)

    def test_batch_validates_inphared_once_and_compares_each_genome_once(self):
        resolved = {"path": "/db/reference_phage_genomes.fna", "version": "2026-04-07",
                    "runtime_paths": {"reference_fasta": "/db/reference_phage_genomes.fna",
                                      "mash_index": "/db/inphared.msh", "metadata": "/db/genome_metadata.tsv",
                                      "manifest": "/db/genome_manifest.json"}}
        comparison = {"status": "COMPLETE", "matches": [], "unique_reference_summaries": []}
        with tempfile.TemporaryDirectory() as temp, \
             patch("phagemine.batch.resolve_validated_inphared", return_value=(resolved, "ready")) as resolve, \
             patch("phagemine.pipeline.resolve_validated_inphared") as pipeline_resolve, \
             patch("phagemine.inphared.compare_genomes", return_value=comparison) as compare:
            input_dir = Path(temp) / "inputs"; input_dir.mkdir()
            for index in range(3):
                shutil.copy2(ROOT / "examples/demo_phage.fasta", input_dir / f"sample-{index}.fasta")
            rows = batch(input_dir, Path(temp) / "batch", gene_predictor="demo")
            self.assertEqual([row["status"] for row in rows], ["SUCCESS"] * 3)
            self.assertEqual(resolve.call_count, 1)
            self.assertFalse(pipeline_resolve.called)
            self.assertEqual(compare.call_count, 3)

    def test_progress_skip_is_explicit(self):
        stream = StringIO()
        progress = ProgressReporter(stream=stream)
        progress.start("INPHARED genome comparison")
        progress.skip("INPHARED genomes not configured")

        output = stream.getvalue()
        self.assertIn("SKIPPED: INPHARED genome comparison", output)
        self.assertIn("INPHARED genomes not configured", output)
        self.assertEqual(progress.completed, 1)

    def test_progress_quiet_mode(self):
        stream = StringIO()
        progress = ProgressReporter(stream=stream, quiet=True)
        progress.start("Pfam")
        progress.finish("summary")
        self.assertEqual(stream.getvalue(), "")

    @pytest.mark.integration
    def test_progress_does_not_enter_scientific_outputs(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "progress-output"
            run(ROOT / "examples/demo_phage.fasta", output, predictor=DemoORFPredictor(), progress=ProgressReporter(stream=StringIO()))
            self.assertNotIn("RUNNING", (output / "run_manifest.json").read_text())
            self.assertNotIn("DONE", (output / "annotation.tsv").read_text())

    def test_run_defaults_to_current_directory_and_orf_reconciliation(self):
        with tempfile.TemporaryDirectory() as temp:
            fasta=Path(temp)/"ijeoma.fasta"; fasta.write_text(">Ijeoma\nATGAAATAG\n")
            captured={}
            def fake_run(*args, **kwargs):
                captured["output"]=args[1]; captured["reconcile_orfs"]=kwargs["reconcile_orfs"]; return 1
            with patch("phagemine.cli.run", side_effect=fake_run), patch("phagemine.preflight.preflight_resources"), patch("phagemine.cli.Path.cwd", return_value=Path(temp)):
                self.assertEqual(main(["run",str(fasta),"--gene-predictor","demo","--quiet"]),0)
            self.assertEqual(captured["output"],str(Path(temp)/"ijeoma_phagemine_results"))
            self.assertTrue(captured["reconcile_orfs"])
    def test_translation(self):
        self.assertEqual(translate("ATGGCTTAA"), "MA")

    def test_fasta_validation(self):
        genome_id, sequence = read_fasta(ROOT / "examples/demo_phage.fasta")
        self.assertEqual(genome_id, "demo_phage")
        self.assertTrue(sequence.startswith("ATG"))

    def test_orfs_receive_stable_ids(self):
        _, sequence = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs("demo", sequence)
        self.assertEqual(len(proteins), 5)
        self.assertEqual(proteins[0].protein_id, "PM_000001")

    def test_evidence_has_provenance(self):
        evidence = Evidence("domain", "signal", EvidenceLevel.WEAK, "source", "v1")
        self.assertEqual(evidence.source_version, "v1")
        self.assertTrue(evidence.supports)

    @pytest.mark.integration
    def test_full_pipeline_generates_explainable_outputs(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "run"
            count = run(ROOT / "examples/demo_phage.fasta", output, predictor=DemoORFPredictor())
            self.assertEqual(count, 5)
            for filename in ("original_input.fasta", "analysis_genome.fasta", "genome_representation.json", "sequencing_provenance.json", "scientific_validation_status.json", "genes.gff3", "proteins.faa", "annotated_proteins.faa", "cds.fna", "annotation.tsv", "evidence.json", "candidate_ranking.tsv", "hallmark_completeness.tsv", "annotation_review.tsv", "report.md", "report.html", "run_manifest.json", "quality_control.json"):
                self.assertTrue((output / filename).exists(), filename)
            scientific_status = json.loads((output / "scientific_validation_status.json").read_text())
            self.assertFalse(scientific_status["functional_strength_empirically_calibrated"])
            self.assertEqual(scientific_status["feature_scope"]["tRNA"], "NOT_CALLED")
            with (output / "annotation.tsv").open() as handle:
                header = next(csv.reader(handle, delimiter="\t"))
                self.assertEqual(header[:8], ["protein_id", "start", "end", "strand", "length_aa", "gene", "product", "proposed_function"])
                self.assertIn("domain_note", header)
                self.assertIn("confidence", header)
                self.assertIn("confidence_calibrated", header)
            named_headers=[line for line in (output / "annotated_proteins.faa").read_text().splitlines() if line.startswith(">")]
            self.assertEqual(len(named_headers), count)
            self.assertTrue(all(" product=\"" in line and " coordinates=" in line and " strand=" in line for line in named_headers))
            report = (output / "report.md").read_text()
            self.assertIn("Computational hypothesis only", report)
            evidence = json.loads((output / "evidence.json").read_text())
            self.assertTrue(any(item["missing_evidence"] for item in evidence))
            for filename in ("genome.fsa", "features.tbl", "proteins.faa", "sequencing_provenance.json", "validation.json", "provenance.json", "README.txt"):
                self.assertTrue((output / "genbank_submission" / filename).exists(), filename)

    def test_genbank_feature_table_and_validation(self):
        genome_id, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs(genome_id, genome)
        validation = validate(genome_id, genome, proteins, {"input_sha256": "fixture"})
        self.assertTrue(validation["valid"])
        table = feature_table(genome_id, proteins)
        self.assertIn(f"{proteins[0].start}\t{proteins[0].end}\tCDS", table)
        self.assertIn("gnl|PhageMine|PM_000001", table)
        self.assertIn("product\thypothetical protein", table)
        self.assertEqual(validation["provenance"]["input_sha256"], "fixture")

    def test_genbank_reverse_strand_coordinates(self):
        genome_id, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        protein = predict_orfs(genome_id, genome)[0]
        protein.strand = "-"
        protein.cds = protein.cds.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        protein.sequence = translate(protein.cds)
        table = feature_table(genome_id, [protein])
        self.assertIn(f"{protein.end}\t{protein.start}\tCDS", table)
        validation = validate(genome_id, genome, [protein], {"input_sha256": "fixture"})
        # Reversing an otherwise valid CDS changes its biological 3-prime end.
        # Coordinate rendering remains correct, but a complete CDS without a
        # terminal stop must now fail NCBI pre-submission validation.
        self.assertFalse(validation["valid"])
        self.assertTrue(any(item["code"] == "NCBI_CDS_NOSTOP" for item in validation["errors"]))

    def test_genbank_reports_fasta_and_claim_problems(self):
        genome_id, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs(genome_id, genome)
        proteins[0].cds = "ATGAAA"
        invalid = validate(genome_id, genome, proteins, {"input_sha256": "fixture"})
        self.assertFalse(invalid["valid"])
        self.assertIn("cds_fasta_mismatch", {item["code"] for item in invalid["errors"]})

    def test_genbank_warns_on_unsupported_mock_claims(self):
        genome_id, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs(genome_id, genome)
        MockEvidenceBackend().annotate(proteins)
        warnings = validate(genome_id, genome, proteins, {"input_sha256": "fixture"})["warnings"]
        self.assertIn("mock_or_unsupported_function", {item["code"] for item in warnings})

    def test_genbank_cli_is_independently_callable(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "genbank"
            self.assertEqual(main(["genbank", str(ROOT / "examples/demo_phage.fasta"), "--gene-predictor", "demo", "--output", str(output)]), 0)
            self.assertTrue((output / "genbank_submission" / "validation.json").exists())

    def test_phanotate_adapter_parses_coordinates_strands_and_translation(self):
        genome_id, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = PHANOTATEPredictor.parse_output(genome_id, genome, "# START STOP FRAME CONTIG\n1 150 + demo_phage\n311 475 + demo_phage\n")
        self.assertEqual([protein.protein_id for protein in proteins], ["PM_000001", "PM_000002"])
        self.assertEqual(proteins[0].gene_call_source, "PHANOTATE")
        self.assertEqual(proteins[0].locus_tag, "PM_000001")
        self.assertEqual(proteins[0].start_codon, "ATG")
        self.assertEqual(proteins[0].stop_codon, "TAA")
        self.assertEqual(proteins[0].sequence, translate(genome[:150]))
        self.assertEqual(proteins[1].gene_call_parameters["reported_strand"], "+")

    def test_phanotate_adapter_handles_reverse_strand(self):
        sequence = "TTATTTCAT"
        proteins = PHANOTATEPredictor.parse_output("reverse", sequence, "9 1 - reverse\n")
        self.assertEqual(proteins[0].strand, "-")
        self.assertEqual(proteins[0].cds, "ATGAAATAA")
        self.assertEqual(proteins[0].sequence, "MK")

    def test_phanotate_unavailable_gives_actionable_error(self):
        predictor = PHANOTATEPredictor("/definitely/not/a/phanotate")
        with self.assertRaisesRegex(RuntimeError, "PHANOTATE is required"):
            predictor.predict("demo", "ATGAAATAA", ROOT / "examples/demo_phage.fasta")

    def test_original_representation_preserves_authoritative_sequence(self):
        representation = GenomeRepresentation.original("assembly", "ATGCCCTAA")
        self.assertEqual(representation.original_sequence, "ATGCCCTAA")
        self.assertEqual(representation.analysis_sequence, "ATGCCCTAA")
        self.assertEqual(representation.topology, Topology.UNKNOWN)
        self.assertEqual(representation.orientation, Orientation.ORIGINAL)
        self.assertEqual(representation.rotation, Rotation.NONE)
        self.assertEqual(representation.transform_history, ())

    def test_reverse_complement_representation_has_explicit_history(self):
        representation = GenomeRepresentation.original("assembly", "ATGCCCTAA").with_reverse_complement("User requested comparison to reference orientation", [{"type": "reference_alignment", "status": "supporting"}], {"accession": "REF_1"})
        self.assertEqual(representation.original_sequence, "ATGCCCTAA")
        self.assertEqual(representation.analysis_sequence, "TTAGGGCAT")
        self.assertEqual(representation.orientation, Orientation.REVERSE_COMPLEMENT)
        self.assertEqual(representation.transform_history[0].operation, "reverse_complement")
        self.assertEqual(representation.manifest()["reference"]["accession"], "REF_1")

    def test_circular_rotation_representation_records_coordinate_scope(self):
        representation = GenomeRepresentation.original("assembly", "AAACCCGGG", topology=Topology.CIRCULAR).with_rotation(4, "User selected coordinate 4 as analysis origin")
        self.assertEqual(representation.analysis_sequence, "CCCGGGAAA")
        self.assertEqual(representation.rotation, Rotation.ROTATED)
        self.assertEqual(representation.transform_history[0].parameters["analysis_origin"], 4)
        with self.assertRaises(ValueError):
            GenomeRepresentation.original("linear", "AAACCC").with_rotation(2, "not allowed")

    def test_combined_transform_history_and_analysis_coordinate_scope(self):
        representation = GenomeRepresentation.original("assembly", "ATGAAATAA", topology=Topology.CIRCULAR).with_reverse_complement("explicit test transform").with_rotation(2, "explicit test origin")
        self.assertEqual([event.operation for event in representation.transform_history], ["reverse_complement", "rotate"])
        self.assertEqual(representation.analysis_sequence, "TATTTCATT")
        self.assertIn("analysis_sequence_id", representation.manifest()["coordinate_scope"])

    @pytest.mark.integration
    def test_pipeline_preserves_original_fasta_and_coordinates_use_analysis_representation(self):
        genome_id, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        representation = GenomeRepresentation.original(genome_id, genome, topology=Topology.CIRCULAR).with_rotation(151, "Fixture-only explicit rotation")
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "represented"
            run(ROOT / "examples/demo_phage.fasta", output, predictor=DemoORFPredictor(), representation=representation)
            self.assertEqual((output / "original_input.fasta").read_bytes(), (ROOT / "examples/demo_phage.fasta").read_bytes())
            manifest = json.loads((output / "run_manifest.json").read_text())
            self.assertEqual(manifest["genome_representation"]["analysis_sequence_id"], representation.analysis_sequence_id)
            gff = (output / "genes.gff3").read_text()
            self.assertIn(representation.analysis_sequence_id + "\tPhageMine\tCDS", gff)
            self.assertIn("Analysis sequence", (output / "report.md").read_text())
            annotations = list(csv.DictReader((output / "annotation.tsv").open(), delimiter="\t"))
            self.assertTrue(all(row["start"] and row["end"] and row["strand"] for row in annotations))

    def test_illumina_provenance(self):
        provenance = SequencingProvenance.from_dict({"sequencing_platform": "ILLUMINA", "assembler": "SPAdes"})
        self.assertEqual(provenance.sequencing_platform, SequencingPlatform.ILLUMINA)
        self.assertEqual(provenance.assembler, "SPAdes")

    def test_oxford_nanopore_provenance(self):
        provenance = SequencingProvenance.from_dict({"sequencing_platform": "OXFORD_NANOPORE", "polishing_method": "Medaka"})
        self.assertEqual(provenance.sequencing_platform, SequencingPlatform.OXFORD_NANOPORE)
        self.assertEqual(provenance.polishing_method, "Medaka")

    def test_pacbio_provenance(self):
        provenance = SequencingProvenance.from_dict({"sequencing_platform": "PACBIO", "read_type": "HiFi"})
        self.assertEqual(provenance.sequencing_platform, SequencingPlatform.PACBIO)
        self.assertEqual(provenance.read_type, "HiFi")

    def test_hybrid_provenance(self):
        provenance = SequencingProvenance.from_dict({"sequencing_platform": "HYBRID", "assembly_method": "hybrid de novo"})
        self.assertEqual(provenance.sequencing_platform, SequencingPlatform.HYBRID)
        self.assertEqual(provenance.assembly_method, "hybrid de novo")

    @pytest.mark.integration
    def test_unknown_provenance_does_not_invent_metadata(self):
        provenance = SequencingProvenance()
        self.assertEqual(provenance.sequencing_platform, SequencingPlatform.UNKNOWN)
        self.assertIsNone(provenance.assembler)
        self.assertIsNone(provenance.raw_reads_available)

    def test_pfam_domtblout_parser_preserves_domain_evidence(self):
        _, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs("demo", genome)
        proteins[0].protein_id = "PM_000001"
        proteins[0].sequence = "M" * 104
        adapter = PfamHMMAdapter(database_version="Pfam-test-1")
        evidence = adapter.parse_domtblout((ROOT / "tests/fixtures/pfam_domtblout.txt").read_text(), proteins, {"status": "REAL", "fixture": True})
        self.assertEqual(len(evidence), 3)
        record = next(item for item in evidence if item.provenance["protein_id"] == "PM_000001")
        self.assertEqual(record.identifier, "PF01813.21")
        self.assertEqual(record.family_name, "ATP-synt_D")
        self.assertEqual(record.source, "Pfam")
        self.assertEqual(record.status, "REAL")
        self.assertEqual(record.evidence_strength, "WEAK")
        self.assertEqual(record.description, "ATP synthase subunit D")
        self.assertEqual(record.coordinates, {"start": 12, "end": 62})
        self.assertAlmostEqual(record.metrics["independent_domain_e_value"], 6.7e-7)
        self.assertAlmostEqual(record.metrics["conditional_domain_e_value"], 0.00082)
        self.assertEqual(record.metrics["query_protein_id"], "PM_000001")
        self.assertEqual(record.metrics["envelope_coordinates"], {"start": 10, "end": 102})
        self.assertEqual(record.provenance["protein_id"], "PM_000001")
        weak = next(item for item in evidence if item.provenance["protein_id"] == "PM_000002")
        self.assertEqual(weak.evidence_strength, "WEAK")
        self.assertFalse("Putative" in weak.description or "protein" == weak.description)

    def test_pfam_thresholds_retain_rejected_raw_hits(self):
        _, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs("demo", genome)
        adapter = PfamHMMAdapter(evalue_threshold=1e-5, coverage_threshold=0.4)
        evidence = adapter.parse_domtblout((ROOT / "tests/fixtures/pfam_domtblout.txt").read_text(), proteins, {"status": "REAL"})
        rejected = next(item for item in evidence if item.provenance["protein_id"] == "PM_000002")
        self.assertEqual(rejected.evidence_strength, "REJECTED")
        self.assertFalse(rejected.supports)
        self.assertIn("raw_hmmscan_row", rejected.metrics)

    def test_pfam_ga_mode_classifies_hits_strong_and_records_mode(self):
        adapter = PfamHMMAdapter(threshold_mode="GA", hmmscan="hmmscan")
        self.assertEqual(adapter.threshold_mode, "GA")
        self.assertTrue(adapter.trusted_cutoff)
        self.assertIn("threshold_mode", adapter.provenance())
        self.assertEqual(adapter.provenance()["threshold_mode"], "GA")
        self.assertTrue(adapter.provenance()["trusted_cutoff"])
        command = adapter._search_command(Path("domtblout"), Path("proteins.faa"))
        self.assertIn("--cut_ga", command)

    def test_pfam_ga_parser_marks_reported_hits_strong_and_mining_can_use_them(self):
        _, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs("demo", genome)
        proteins[0].sequence = "M" * 104
        adapter = PfamHMMAdapter(threshold_mode="GA", database_version="Pfam-test-1")
        evidence = adapter.parse_domtblout((ROOT / "tests/fixtures/pfam_domtblout.txt").read_text(), proteins)
        ga = next(item for item in evidence if item.provenance["protein_id"] == "PM_000001")
        self.assertEqual(ga.evidence_strength, "STRONG")
        self.assertTrue(ga.supports)
        proteins[0].evidence.append(ga)
        mine(proteins)
        self.assertGreater(proteins[0].biological_interest, 0)

    def test_short_orf_qc_evidence_is_idempotent(self):
        protein = Protein(
            "demo", "PM_SHORT", 1, 90, "+",
            "ATG" + ("AAA" * 28) + "TAA",
            "M" * 29,
            "PHANOTATE",
        )

        mine([protein])
        mine([protein])

        qc = [
            evidence for evidence in protein.evidence
            if evidence.modality == "quality"
            and evidence.source == "phagemine-qc"
        ]

        self.assertEqual(len(qc), 1)
        self.assertFalse(qc[0].supports)
        self.assertEqual(qc[0].metrics["protein_length"], 29)
        self.assertEqual(protein.score_components["quality_penalty"], 15)

    def test_weak_only_pfam_does_not_enable_ranking(self):
        _, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs("demo", genome)
        adapter = PfamHMMAdapter(threshold_mode="NONE")
        evidence = adapter.parse_domtblout((ROOT / "tests/fixtures/pfam_domtblout.txt").read_text(), proteins)
        proteins[0].evidence.append(next(item for item in evidence if item.provenance["protein_id"] == "PM_000001"))
        mine(proteins)
        self.assertEqual(proteins[0].biological_interest, 0)
        self.assertEqual(proteins[0].evidence_diversity, "Low")

    @pytest.mark.integration
    def test_pipeline_persists_accepted_ga_pfam_evidence(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "ga-persistence"
            _, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
            proteins = predict_orfs("demo", genome)
            accepted = Evidence(
                "domain", "Pfam domain evidence", EvidenceLevel.COMPUTATIONAL,
                "Pfam", "Pfam-test", status="REAL", supports=True,
                identifier="PF00001", evidence_strength="STRONG",
                family_name="Test_family", description="Test domain",
                metrics={"query_protein_id": proteins[0].protein_id,
                         "independent_domain_e_value": 1e-8},
                provenance={"threshold_mode": "GA", "trusted_cutoff": True},
            )
            result = EvidenceAdapterResult(
                "PfamHMMAdapter", "REAL", evidence=[accepted],
                provenance={"threshold_mode": "GA", "trusted_cutoff": True},
            )
            with patch("phagemine.pipeline.PfamHMMAdapter") as adapter:
                adapter.return_value.name = PfamHMMAdapter.name
                adapter.return_value.analyze.return_value = result
                run(ROOT / "examples/demo_phage.fasta", output,
                    predictor=DemoORFPredictor(), pfam_threshold_mode="GA")
            records = json.loads((output / "evidence.json").read_text())
            persisted = [e for p in records for e in p["evidence"] if e["source"] == "Pfam"]
            self.assertEqual(len(persisted), 1)
            self.assertEqual(persisted[0]["status"], "REAL")
            self.assertEqual(persisted[0]["provenance"]["threshold_mode"], "GA")
            self.assertTrue(persisted[0]["provenance"]["trusted_cutoff"])
            manifest = json.loads((output / "checkpoints" / "evidence_complete" / "checkpoint_manifest.json").read_text())
            self.assertEqual(manifest["evidence_adapters"][0]["adapter"], PfamHMMAdapter.name)

    def test_pfam_unavailable_has_no_fabricated_evidence(self):
        result = PfamHMMAdapter("/definitely/missing/Pfam-A.hmm", "/definitely/missing/hmmscan").analyze([])
        self.assertEqual(result.status, "UNAVAILABLE")
        self.assertEqual(result.evidence, [])
        self.assertIn("no domain evidence was fabricated", result.message)

    def test_vog_parser_maps_orthology_and_annotations(self):
        _, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs("demo", genome)
        with tempfile.TemporaryDirectory() as temp:
            annotations = Path(temp) / "vog.annotations.tsv"
            annotations.write_text("#GroupName\tProteinCount\tSpeciesCount\tFunctionalCategory\tConsensusFunctionalDescription\nVOG00001\t1\t1\tS\tstructural viral protein\n")
            adapter = VOGHMMAdapter(annotations_path=annotations, evalue_threshold=1e-5, coverage_threshold=0.5)
            evidence = adapter.parse_domtblout((ROOT / "tests/fixtures/vog_domtblout.txt").read_text(), proteins, {"status": "REAL"})
        strong = next(item for item in evidence if item.identifier == "VOG00001")
        self.assertEqual(strong.source, "VOGDB")
        self.assertEqual(strong.evidence_strength, "STRONG")
        self.assertEqual(strong.metrics["query_coordinates"], {"start": 8, "end": 96})
        self.assertEqual(strong.metrics["vog_id"], "VOG00001")
        self.assertEqual(strong.metrics["functional_category"], "S")
        self.assertEqual(strong.description, "structural viral protein")
        self.assertEqual(strong.metrics["consensus_functional_description"], "structural viral protein")
        self.assertEqual(strong.provenance["threshold_mode"], "MANUAL")

    def test_vog_unavailable_does_not_fabricate_evidence(self):
        result = VOGHMMAdapter("/definitely/missing/vog.hmm").analyze([])
        self.assertEqual(result.status, "UNAVAILABLE")
        self.assertEqual(result.evidence, [])
        self.assertIn("no VOG evidence was fabricated", result.message)

    def test_vog_rejected_rows_are_not_supportive(self):
        _, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs("demo", genome)
        evidence = VOGHMMAdapter(evalue_threshold=1e-5, coverage_threshold=0.5).parse_domtblout(
            (ROOT / "tests/fixtures/vog_domtblout.txt").read_text(), proteins)
        accepted = next(item for item in evidence if item.identifier == "VOG00001")
        rejected = next(item for item in evidence if item.identifier == "VOG00002")
        self.assertTrue(accepted.supports)
        self.assertEqual(accepted.evidence_strength, "STRONG")
        self.assertFalse(rejected.supports)
        self.assertEqual(rejected.evidence_strength, "REJECTED")

    def test_registered_vog_version_propagates(self):
        with tempfile.TemporaryDirectory() as temp:
            registry = Path(temp) / "resources.json"
            hmm = Path(temp) / "VOGDB.hmm"
            hmm.write_text("synthetic")
            manager = EvidenceResourceManager(registry)
            manager.register("VOGDB-236", ResourceType.VOGDB, hmm, version="236")
            registered = manager.find(ResourceType.VOGDB)
            adapter = VOGHMMAdapter(registered["path"], database_version=registered["version"], hmmscan="/missing/hmmscan")
            self.assertEqual(adapter.provenance()["vogdb_version"], "236")
            _, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
            proteins = predict_orfs("demo", genome)
            row = adapter.parse_domtblout((ROOT / "tests/fixtures/vog_domtblout.txt").read_text(), proteins)[0]
            self.assertEqual(row.source_version, "236")
            self.assertEqual(row.provenance["vogdb_version"], "236")

    def test_vog_is_distinct_evidence_source(self):
        evidence = Evidence("viral_orthology", "VOG evidence", EvidenceLevel.COMPUTATIONAL, "VOGDB", "test", status="REAL", evidence_strength="STRONG", provenance={"protein_id": "PM_000001"})
        self.assertNotEqual(evidence.source, "Pfam")

    def test_phrogs_parser_preserves_phrogs_semantics(self):
        _, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs("demo", genome)
        with tempfile.TemporaryDirectory() as temp:
            annotations = Path(temp) / "phrogs.tsv"
            annotations.write_text("phrog\tcolor\tannot\tcategory\n1\t#fff\tunknown function\tunknown\n2\t#000\tDNA-associated protein\tDNA, RNA and nucleotide metabolism\n")
            adapter = PHROGSMMseqsAdapter(annotations_path=annotations, evalue_threshold=1e-5, coverage_threshold=0.5)
            evidence = adapter.parse_tabular((ROOT / "tests/fixtures/phrogs_mmseqs.tsv").read_text(), proteins)
        strong = next(item for item in evidence if item.identifier == "1")
        rejected = next(item for item in evidence if item.identifier == "2")
        self.assertEqual(strong.modality, "phage_orthology")
        self.assertEqual(strong.metrics["functional_category"], "unknown")
        self.assertIsNone(strong.description)
        self.assertTrue(strong.supports)
        self.assertEqual(strong.metrics["mmseqs_score"], 75.0)
        self.assertEqual(strong.metrics["sequence_identity"], 32.5)
        self.assertEqual(strong.metrics["bit_score"], 75.0)
        self.assertEqual(strong.metrics["percent_identity"], 32.5)
        self.assertEqual(strong.metrics["raw_mmseqs_row"].split("\t")[0], "PM_000001")
        self.assertEqual(adapter.provenance()["output_format"], MMSEQS_FORMAT)
        self.assertEqual(strong.provenance["threshold_mode"], "MANUAL")
        self.assertEqual(rejected.evidence_strength, "REJECTED")
        self.assertFalse(rejected.supports)

    def test_phrogs_unavailable_does_not_fabricate_evidence(self):
        result = PHROGSMMseqsAdapter("/missing/phrogs", mmseqs="/missing/mmseqs").analyze([])
        self.assertEqual(result.status, "UNAVAILABLE")
        self.assertEqual(result.evidence, [])

    def test_phrogs_manual_identity_and_alignment_thresholds_reject(self):
        _, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs("demo", genome)
        adapter = PHROGSMMseqsAdapter(
            evalue_threshold=None, coverage_threshold=None,
            identity_threshold=40.0, alignment_length_threshold=80)
        evidence = adapter.parse_tabular(
            (ROOT / "tests/fixtures/phrogs_mmseqs.tsv").read_text(), proteins)
        self.assertTrue(evidence)
        self.assertTrue(all(not item.supports for item in evidence))
        self.assertTrue(all(item.evidence_strength == "REJECTED" for item in evidence))

    def test_phrogs_conflicting_strong_annotations_are_preserved(self):
        _, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs("demo", genome)
        with tempfile.TemporaryDirectory() as temp:
            annotations = Path(temp) / "phrogs.tsv"
            annotations.write_text(
                "phrog\tcolor\tannot\tcategory\n"
                "1\t#fff\tintegrase\tintegration and excision\n"
                "2\t#000\tportal protein\thead and packaging\n")
            adapter = PHROGSMMseqsAdapter(annotations_path=annotations)
            text = (
                "PM_000001\t1\t50\t90\t1\t90\t104\t1\t90\t100\t1e-20\t100\t0.86\t0.9\n"
                "PM_000001\t2\t45\t85\t2\t86\t104\t3\t87\t120\t1e-15\t90\t0.82\t0.71\n")
            evidence = adapter.parse_tabular(text, proteins)
        adapter._mark_conflicts(evidence)
        self.assertEqual({item.description for item in evidence}, {"integrase", "portal protein"})
        self.assertTrue(all(item.metrics["conflict"] for item in evidence))

    def test_swissprot_parser_metadata_and_strengths(self):
        _, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs("demo", genome)
        with tempfile.TemporaryDirectory() as temp:
            dat = Path(temp) / "uniprot.dat"
            dat.write_text("ID   001R_FRG3G              Reviewed;         100 AA.\nAC   Q6GZX4;\nDE   RecName: Full=Curated viral protein;\nOS   Test virus.\nOX   NCBI_TaxID=12345;\nGN   Name=geneA;\nPE   1: Evidence at protein level;\n//\nID   002L_FRG3G              Reviewed;         100 AA.\nAC   Q6GZX3;\nDE   RecName: Full=Second curated protein;\n//\n")
            adapter = SwissProtEvidenceAdapter(metadata_path=dat)
            evidence = adapter.parse_tabular((ROOT / "tests/fixtures/swissprot_diamond.tsv").read_text(), proteins)
        strong = next(item for item in evidence if item.identifier == "Q6GZX4")
        moderate = next(item for item in evidence if item.identifier == "Q6GZX3")
        rejected = next(item for item in evidence if item.provenance["protein_id"] == "PM_000003")
        self.assertEqual(strong.evidence_strength, "STRONG")
        self.assertTrue(strong.supports)
        self.assertEqual(strong.metrics["entry_name"], "001R_FRG3G")
        self.assertEqual(strong.metrics["protein_name"], "Curated viral protein")
        self.assertEqual(strong.metrics["taxonomy_id"], "12345")
        self.assertAlmostEqual(strong.metrics["query_coverage"], 0.8)
        self.assertAlmostEqual(strong.metrics["subject_coverage"], 0.8)
        self.assertAlmostEqual(strong.metrics["evalue"], 1e-30)
        self.assertLessEqual(strong.metrics["query_coverage"], 1.0)
        self.assertLessEqual(strong.metrics["subject_coverage"], 1.0)
        self.assertEqual(moderate.evidence_strength, "MODERATE")
        self.assertTrue(moderate.supports)
        self.assertEqual(rejected.evidence_strength, "REJECTED")
        self.assertFalse(rejected.supports)

    def test_swissprot_full_length_coverage_is_one(self):
        _, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs("demo", genome)
        row = "PM_000001\tsp|Q6GZX4|001R_FRG3G\t80.0\t120\t120\t120\t1\t120\t1\t120\t1e-30\t120.0\n"
        record = SwissProtEvidenceAdapter().parse_tabular(row, proteins)[0]
        self.assertEqual(record.metrics["query_coverage"], 1.0)
        self.assertEqual(record.metrics["subject_coverage"], 1.0)

    def test_swissprot_subject_id_and_unavailable(self):
        self.assertEqual(SwissProtEvidenceAdapter.parse_subject_id("sp|Q6GZX4|001R_FRG3G"), ("Q6GZX4", "001R_FRG3G"))
        result = SwissProtEvidenceAdapter("/missing/db.dmnd", diamond="/missing/diamond").analyze([])
        self.assertEqual(result.status, "UNAVAILABLE")
        self.assertEqual(result.evidence, [])

    def test_swissprot_diamond_fields_precede_threads_option(self):
        adapter = SwissProtEvidenceAdapter("swissprot.dmnd", diamond="diamond", threads=4)
        command = adapter._search_command(Path("proteins.faa"), Path("diamond.tsv"))
        self.assertLess(command.index("bitscore"), command.index("--threads"))
        self.assertEqual(command[command.index("--threads") + 1], "4")

    @pytest.mark.integration
    def test_pipeline_manifest_marks_pfam_unavailable_and_no_mock_by_default(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "real-no-db"
            run(ROOT / "examples/demo_phage.fasta", output, predictor=DemoORFPredictor(), pfam_path="/definitely/missing/Pfam-A.hmm")
            manifest = json.loads((output / "run_manifest.json").read_text())
            self.assertEqual(manifest["evidence_adapters"][0]["status"], "UNAVAILABLE")
            self.assertNotIn("MOCK", {adapter["status"] for adapter in manifest["evidence_adapters"]})
            evidence = json.loads((output / "evidence.json").read_text())
            self.assertFalse(any(item["source"] == "mock-phage-evidence" for protein in evidence for item in protein["evidence"]))

    @pytest.mark.integration
    def test_no_evidence_does_not_assign_misleading_ordinal_ranks(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "insufficient"
            run(ROOT / "examples/demo_phage.fasta", output, predictor=DemoORFPredictor())
            manifest = json.loads((output / "run_manifest.json").read_text())
            self.assertEqual(manifest["discovery_ranking"]["status"], "INSUFFICIENT_EVIDENCE")
            ranking = (output / "candidate_ranking.tsv").read_text().splitlines()
            self.assertTrue(all(line.startswith("NA\t") for line in ranking[1:]))
            self.assertIn("INSUFFICIENT_EVIDENCE", (output / "report.md").read_text())

    def test_resource_manager_registers_pfam_and_persists_metadata(self):
        with tempfile.TemporaryDirectory() as temp:
            registry = Path(temp) / "resources.json"
            database = Path(temp) / "Pfam-A.hmm"
            database.write_text("synthetic HMM")
            manager = EvidenceResourceManager(registry)
            manager.register("PFAM", ResourceType.PFAM, database, version="Pfam-test-1", required_tools=["hmmscan"])
            loaded = EvidenceResourceManager(registry).get("PFAM")
            self.assertEqual(loaded["resource_type"], "PFAM")
            self.assertEqual(loaded["status"], ResourceStatus.READY.value)
            self.assertEqual(loaded["version"], "Pfam-test-1")
            self.assertIsNone(loaded["checksum"])

    def test_resource_manager_supports_types_missing_paths_and_removal(self):
        with tempfile.TemporaryDirectory() as temp:
            registry = Path(temp) / "resources.json"
            manager = EvidenceResourceManager(registry)
            manager.register("SWISSPROT", ResourceType.SWISSPROT, Path(temp) / "missing.fasta")
            manager.register("REFSEQ", "REFSEQ", Path(temp) / "missing-refseq")
            self.assertEqual({item["resource_type"] for item in manager.list()}, {"SWISSPROT", "REFSEQ"})
            self.assertTrue(all(item["status"] == "UNAVAILABLE" for item in manager.list()))
            self.assertTrue(manager.unregister("SWISSPROT"))
            self.assertFalse(manager.unregister("SWISSPROT"))

    def test_resource_manager_does_not_download_or_mutate_database(self):
        with tempfile.TemporaryDirectory() as temp:
            registry = Path(temp) / "resources.json"
            database = Path(temp) / "dummy.hmm"
            database.write_text("unchanged")
            manager = EvidenceResourceManager(registry)
            manager.register("PFAM", "PFAM", database)
            self.assertEqual(database.read_text(), "unchanged")
            self.assertFalse((Path(temp) / "download").exists())

    def test_pfam_registered_fallback_and_explicit_precedence(self):
        with tempfile.TemporaryDirectory() as temp:
            registry = Path(temp) / "resources.json"
            registered = Path(temp) / "registered.hmm"
            explicit = Path(temp) / "explicit.hmm"
            registered.write_text("registered")
            explicit.write_text("explicit")
            manager = EvidenceResourceManager(registry)
            manager.register("PFAM", "PFAM", registered)
            found = manager.find(ResourceType.PFAM)
            self.assertEqual(found["path"], str(registered))
            # The pipeline's precedence decision is explicit path first; this mirrors the adapter input contract.
            explicit_adapter = PfamHMMAdapter(explicit, "/missing/hmmscan")
            self.assertEqual(explicit_adapter.pfam_path, explicit)

    def test_default_registry_is_user_level_not_repository(self):
        registry = default_registry_path()
        self.assertNotIn(str(ROOT), str(registry))

    def test_registry_environment_override_is_atomic_and_keeps_backup(self):
        with tempfile.TemporaryDirectory() as temp:
            registry = Path(temp) / "isolated" / "resources.json"
            with patch.dict(os.environ, {"PHAGEMINE_REGISTRY_PATH": str(registry)}):
                self.assertEqual(default_registry_path(), registry)
                resource = Path(temp) / "resource.hmm"
                resource.write_text("data")
                manager = EvidenceResourceManager()
                manager.register("first", "PFAM", resource)
                manager.register("second", "VOGDB", resource)
            self.assertTrue(registry.is_file())
            self.assertTrue(registry.with_suffix(".json.bak").is_file())

    @pytest.mark.integration
    def test_fasta_analysis_without_sequencing_metadata_uses_unknown(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "unknown-provenance"
            run(ROOT / "examples/demo_phage.fasta", output, predictor=DemoORFPredictor())
            manifest = json.loads((output / "run_manifest.json").read_text())
            self.assertEqual(manifest["sequencing_provenance"]["sequencing_platform"], "UNKNOWN")

    @pytest.mark.integration
    def test_sequencing_provenance_is_preserved_and_separate_from_representation(self):
        genome_id, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        representation = GenomeRepresentation.original(genome_id, genome)
        provenance = SequencingProvenance.from_dict({"sequencing_platform": "OXFORD_NANOPORE", "assembler": "Flye", "raw_reads_available": True, "metadata_source": "lab notebook"})
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "ont"
            run(ROOT / "examples/demo_phage.fasta", output, predictor=DemoORFPredictor(), representation=representation, sequencing_provenance=provenance)
            manifest = json.loads((output / "run_manifest.json").read_text())
            self.assertEqual(manifest["sequencing_provenance"]["sequencing_platform"], "OXFORD_NANOPORE")
            self.assertEqual(manifest["genome_representation"]["topology"], "UNKNOWN")
            self.assertEqual(manifest["genome_representation"]["orientation"], "ORIGINAL")
            package = json.loads((output / "genbank_submission" / "sequencing_provenance.json").read_text())
            self.assertEqual(package["assembler"], "Flye")

    def _complete_metadata(self):
        return SubmissionMetadata.from_dict(json.loads((ROOT / "examples/submission_metadata.json").read_text()))

    def _fake_table2asn(self, directory, exit_code=0):
        script = Path(directory) / f"fake_table2asn_{exit_code}.sh"
        script.write_text(f"#!/bin/sh\nmkdir -p \"$4\"\necho validation > \"$4/fake.val\"\necho stats > \"$4/fake.stats\"\nexit {exit_code}\n")
        script.chmod(script.stat().st_mode | stat.S_IXUSR)
        return str(script)

    def _submission_acceptable_proteins(self, genome_id, genome):
        proteins = predict_orfs(genome_id, genome)
        for protein in proteins:
            protein.annotation = "capsid protein"
            protein.annotation_level = EvidenceLevel.CURATED
        return proteins

    def test_missing_metadata_is_not_submission_ready(self):
        genome_id, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        with tempfile.TemporaryDirectory() as temp:
            result = write_package(temp, genome_id, genome, predict_orfs(genome_id, genome), {"input_sha256": "fixture"}, table2asn_executable="not-installed-table2asn")
            self.assertEqual(result["submission_readiness"], "INCOMPLETE_METADATA")
            self.assertFalse(result["ncbi_table2asn_validation"]["official_ncbi_validation"])

    def test_complete_metadata_and_table2asn_success_are_ready(self):
        genome_id, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        with tempfile.TemporaryDirectory() as temp:
            result = write_package(temp, genome_id, genome, self._submission_acceptable_proteins(genome_id, genome), {"input_sha256": "fixture"}, self._complete_metadata(), self._fake_table2asn(temp))
            self.assertEqual(result["submission_readiness"], "READY")
            self.assertEqual(result["ncbi_table2asn_validation"]["state"], "passed")
            self.assertIn("table2asn_output/fake.val", result["ncbi_table2asn_validation"]["output_files"])

    def test_table2asn_failure_fails_readiness(self):
        genome_id, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        with tempfile.TemporaryDirectory() as temp:
            result = write_package(temp, genome_id, genome, predict_orfs(genome_id, genome), {"input_sha256": "fixture"}, self._complete_metadata(), self._fake_table2asn(temp, 1))
            self.assertEqual(result["submission_readiness"], "VALIDATION_FAILED")
            self.assertEqual(result["ncbi_table2asn_validation"]["state"], "failed")

    def test_table2asn_unavailable_is_explicit(self):
        result = table2asn_status(Path(tempfile.mkdtemp()), "not-installed-table2asn")
        self.assertEqual(result["state"], "unavailable")
        self.assertFalse(result["official_ncbi_validation"])

    def test_conservative_product_and_cds_provenance(self):
        genome_id, genome = read_fasta(ROOT / "examples/demo_phage.fasta")
        proteins = predict_orfs(genome_id, genome)
        MockEvidenceBackend().annotate(proteins)
        with tempfile.TemporaryDirectory() as temp:
            write_package(temp, genome_id, genome, proteins, {"input_sha256": "fixture"})
            package = Path(temp) / "genbank_submission"
            self.assertEqual(product_name(proteins[0]), "hypothetical protein")
            self.assertNotIn("Putative structural", (package / "features.tbl").read_text())
            provenance = json.loads((package / "cds_provenance.json").read_text())
            self.assertTrue(provenance["PM_000001"]["evidence"])
