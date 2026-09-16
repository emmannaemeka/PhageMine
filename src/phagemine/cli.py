from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .pipeline import run
from .resume import reclassify, resume
from .models import SubmissionMetadata
from .resources import EvidenceResourceManager, ResourceType
from .compare import compare
from .batch import batch
from .progress import ProgressReporter
import json


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="PhageMine: annotation followed by cautious discovery mining",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""After installing PhageMine, install and verify the evidence databases:
  phagemine databases install --all
  phagemine doctor

Individual installers are also available:
  phagemine databases install pfam
  phagemine databases install vogdb
  phagemine databases install swissprot
  phagemine databases install phrogs
  phagemine databases install pmfdb
  phagemine databases install inphared""",
    )
    from . import __version__
    parser.add_argument("--version", action="version", version=__version__)
    subcommands = parser.add_subparsers(dest="command")
    gui_command = subcommands.add_parser("gui", help="Launch the optional local graphical interface")
    doctor_command = subcommands.add_parser("doctor", help="Validate executables and registered evidence resources")
    doctor_command.add_argument("--output", help="Optional JSON report path")
    doctor_command.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    doctor_command.add_argument("--deep", action="store_true", help="Open PHROGs databases and verify their operational formats")
    databases = subcommands.add_parser("databases", help="Manage registered local evidence resources")
    database_commands = databases.add_subparsers(dest="database_command", required=True)
    register = database_commands.add_parser("register", help="Register a local evidence resource")
    register.add_argument("resource_type", choices=[item.value.lower() for item in ResourceType])
    register.add_argument("path")
    register.add_argument("--name")
    register.add_argument("--version")
    register.add_argument("--checksum")
    register.add_argument("--notes")
    register.add_argument("--annotations", help="Optional annotation table path (used by VOGDB resources)")
    register.add_argument("--metadata", help="Optional metadata path (used by Swiss-Prot resources)")
    register.add_argument("--hmm-profiles", help="Optional PHROGs all_phrogs.h3m path")
    install = database_commands.add_parser("install", help="Download, prepare, register, and validate evidence databases")
    install.add_argument("resource", nargs="?", choices=("pfam", "vogdb", "swissprot", "phrogs", "pmfdb", "inphared"))
    install.add_argument("--all", action="store_true", help="Install every annotation and comparative database")
    install.add_argument("--directory", help="Database installation root (default: platform user-data directory)")
    install.add_argument("--force", action="store_true", help="Replace an existing installation of the selected release")
    install.add_argument("--keep-downloads", action="store_true", help="Keep downloaded archives after successful installation")
    install.add_argument("--dry-run", action="store_true", help="Show download estimates without installing")
    install.add_argument("--json", action="store_true", help="Emit installation results as JSON")
    install.add_argument("--threads", type=int, default=1, help="Preparation threads for MMseqs2 indexes")
    database_commands.add_parser("status", help="List registered evidence resources")
    attach_hmm = database_commands.add_parser(
        "attach-phrogs-hmm", help="Attach an existing all_phrogs.h3m to the registered PHROGs resource")
    attach_hmm.add_argument("path")
    remove = database_commands.add_parser("remove", help="Remove a resource registration")
    remove.add_argument("name")
    resume_command = subcommands.add_parser("resume", help="Resume from a completed run and run missing evidence")
    resume_command.add_argument("source")
    resume_command.add_argument("--output", required=True)
    resume_command.add_argument("--run-missing-evidence", action="store_true")
    resume_command.add_argument("--refresh-evidence", choices=("PHROGS",))
    resume_command.add_argument("--mmseqs")
    resume_command.add_argument("--phrogs")
    resume_command.add_argument("--phrogs-annotations")
    resume_command.add_argument("--phrogs-hmm", help="Existing all_phrogs.h3m profile database")
    resume_command.add_argument("--phrogs-evalue", type=float, default=1e-5)
    resume_command.add_argument("--phrogs-coverage", type=float, default=0.5)
    resume_command.add_argument("--phrogs-score", type=float)
    resume_command.add_argument("--threads", type=int, default=1)
    reclassify_command = subcommands.add_parser(
        "reclassify",
        help="Regenerate annotations from persisted evidence without database searches",
    )
    reclassify_command.add_argument("source")
    reclassify_command.add_argument("--output", required=True)
    stage_command = subcommands.add_parser("resume-stage", help="Recover from validated mixed stage checkpoints")
    stage_command.add_argument("source"); stage_command.add_argument("--swissprot"); stage_command.add_argument("--diamond")
    compare_command = subcommands.add_parser("compare", help="Compare completed PhageMine result directories offline")
    compare_command.add_argument("results", nargs="+", help="Two to ten completed PhageMine result directories")
    compare_command.add_argument("--output", required=True)
    compare_command.add_argument("--mmseqs")
    family_command = subcommands.add_parser("families", help="Offline protein-family workflows")
    family_sub = family_command.add_subparsers(dest='families_command', required=True)
    family_build = family_sub.add_parser('build'); family_build.add_argument('proteins'); family_build.add_argument('--output',required=True); family_build.add_argument('--metadata'); family_build.add_argument('--registry'); family_build.add_argument('--version',default='1.0'); family_build.add_argument('--backend',choices=('exact','mmseqs'),default='exact'); family_build.add_argument('--mmseqs',default='mmseqs'); family_build.add_argument('--minimum-identity',type=float,default=0.3); family_build.add_argument('--minimum-coverage',type=float,default=0.5); family_build.add_argument('--coverage-mode',type=int,default=0); family_build.add_argument('--clustering-mode',type=int,default=0)
    family_assign = family_sub.add_parser('assign'); family_assign.add_argument('proteins'); family_assign.add_argument('--database',required=True); family_assign.add_argument('--output',required=True); family_assign.add_argument('--mmseqs',default='mmseqs')
    family_compare = family_sub.add_parser('compare'); family_compare.add_argument('database'); family_compare.add_argument('--output',required=True); family_compare.add_argument('--core-genomes')
    family_enrich = family_sub.add_parser('enrich'); family_enrich.add_argument('database'); family_enrich.add_argument('--results',nargs='+',required=True); family_enrich.add_argument('--output',required=True)
    family_priority = family_sub.add_parser('prioritize'); family_priority.add_argument('database'); family_priority.add_argument('--output',required=True)
    family_external = family_sub.add_parser('validate-external'); family_external.add_argument('database'); family_external.add_argument('--reference-db'); family_external.add_argument('--reference-proteins'); family_external.add_argument('--reference-metadata'); family_external.add_argument('--output',required=True); family_external.add_argument('--mmseqs',default='mmseqs'); family_external.add_argument('--top-n',type=int)
    family_sensitivity = family_sub.add_parser('sensitivity'); family_sensitivity.add_argument('proteins'); family_sensitivity.add_argument('--output',required=True); family_sensitivity.add_argument('--identities',required=True); family_sensitivity.add_argument('--minimum-coverage',type=float,default=0.5); family_sensitivity.add_argument('--coverage-mode',type=int,default=0); family_sensitivity.add_argument('--clustering-mode',type=int,default=0); family_sensitivity.add_argument('--metadata'); family_sensitivity.add_argument('--backend',choices=('exact','mmseqs'),default='mmseqs'); family_sensitivity.add_argument('--mmseqs',default='mmseqs')
    benchmark_command = subcommands.add_parser("benchmark", help="Import and compare external annotation outputs")
    benchmark_command.add_argument("--output", required=True)
    benchmark_command.add_argument("--prokka-gff", help="Prokka GFF3 output")
    benchmark_command.add_argument("--pharokka-gff", help="Pharokka GFF3 output")
    benchmark_command.add_argument("--multiphate-gff", help="multiPhATE2 GFF output")
    benchmark_command.add_argument("--phold-genbank", help="Phold GenBank output")
    benchmark_command.add_argument("--phagemine-results", help="Completed PhageMine result directory")
    truth = benchmark_command.add_mutually_exclusive_group()
    truth.add_argument("--truth-gff", help="Expert-reviewed truth-set GFF3; required for accuracy claims")
    truth.add_argument("--truth-genbank", help="Expert-reviewed truth-set GenBank file; required for accuracy claims")
    batch_command = subcommands.add_parser("batch", help="Process a directory of phage FASTA files")
    batch_command.add_argument("input_dir"); batch_command.add_argument("--output", required=True); batch_command.add_argument("--mode", choices=("annotate", "discover", "both"), default="annotate"); batch_command.add_argument("--recursive", action="store_true"); batch_command.add_argument("--resume-existing", action="store_true"); batch_command.add_argument("--fail-fast", action="store_true"); batch_command.add_argument("--gene-predictor", choices=("phanotate","demo"), default="phanotate"); batch_command.add_argument("--gene-model-policy", choices=("phanotate-only", "consensus"), default="phanotate-only"); batch_command.add_argument("--gene-model-profile", choices=("standard", "extended"), default="standard"); batch_command.add_argument("--molecule-type", choices=("dna", "rna"), default="dna"); batch_command.add_argument("--segmented", action="store_true"); batch_command.add_argument("--phanotate"); batch_command.add_argument("--reconcile-orfs", action="store_true"); batch_command.add_argument("--prodigal"); batch_command.add_argument("--threads", type=int, default=1); batch_command.add_argument("--evidence", dest="evidence_profile", choices=("core", "standard", "full"), default="core")
    extract_command = subcommands.add_parser("extract", help="Retrieve stable protein records and FASTA from a completed run")
    extract_sub = extract_command.add_subparsers(dest="extract_command", required=True)
    extract_protein = extract_sub.add_parser("protein", help="Extract one protein by stable protein ID")
    extract_protein.add_argument("protein_id")
    extract_protein.add_argument("--run", required=True, help="Completed genome or batch result directory")
    extract_protein.add_argument("--protein-fasta", action="store_true")
    extract_protein.add_argument("--cds-fasta", action="store_true")
    extract_protein.add_argument("--evidence", action="store_true")
    for name in ("annotate", "mine", "run", "genbank"):
        command = subcommands.add_parser(name, help=f"Run the MVP {name} workflow")
        command.add_argument("fasta", help="Single-record phage genome FASTA")
        command.add_argument("--output", default=None, help="Output directory (default: ./<genome stem>_phagemine_results)")
        command.add_argument("--metadata", help="Optional JSON submission metadata; missing fields are not inferred")
        command.add_argument("--sequencing-provenance", help="Optional JSON sequencing/assembly provenance; absent values remain UNKNOWN")
        command.add_argument("--table2asn", help="Optional path to official NCBI table2asn executable")
        command.add_argument("--gene-predictor", choices=("phanotate", "demo"), default="phanotate", help="Gene caller; demo is for fixtures/tests only")
        command.add_argument("--gene-model-policy", choices=("phanotate-only", "consensus"), default="phanotate-only", help="Final gene-model policy (consensus is experimental; default preserves legacy behavior)")
        command.add_argument("--gene-model-profile", choices=("standard", "extended"), default="standard", help="Consensus provider profile")
        command.add_argument("--molecule-type", choices=("dna", "rna"), default="dna", help="Declared molecule type; RNA uses Pyrodigal-rv")
        command.add_argument("--segmented", action="store_true", help="Treat all FASTA records as segments of one RNA genome")
        command.add_argument("--phanotate", help="Path to PHANOTATE executable")
        command.add_argument("--reconcile-orfs", action="store_true", help="Compare PHANOTATE models with optional Prodigal models without changing annotation")
        command.add_argument("--prodigal", help="Path to Prodigal executable for ORF reconciliation")
        command.add_argument("--pfam", help="Optional local Pfam HMM database path")
        command.add_argument("--pfam-hmmscan", help="Optional hmmscan executable path")
        command.add_argument("--pfam-evalue", type=float, help="Optional Pfam domain i-Evalue threshold")
        command.add_argument("--pfam-coverage", type=float, help="Optional Pfam query coverage threshold (0-1)")
        command.add_argument("--pfam-trusted-cutoff", action="store_true", help="Use configured Pfam trusted-cutoff policy when available")
        command.add_argument("--pfam-threshold-mode", choices=("ga", "manual", "none"), help="Pfam threshold mode")
        command.add_argument("--vogdb", help="Optional prepared combined VOGDB HMM database path")
        command.add_argument("--vog-annotations", help="Optional VOGDB annotation TSV/TSV.GZ path")
        command.add_argument("--vog-hmmscan", help="Optional hmmscan executable path for VOGDB")
        command.add_argument("--vog-evalue", type=float, default=1e-5, help="VOGDB independent E-value threshold")
        command.add_argument("--vog-coverage", type=float, default=0.5, help="VOGDB query coverage threshold")
        command.add_argument("--swissprot", help="Optional Swiss-Prot DIAMOND database path")
        command.add_argument("--swissprot-metadata", help="Optional Swiss-Prot DAT metadata path")
        command.add_argument("--diamond", help="Optional DIAMOND executable path")
        command.add_argument("--swissprot-evalue", type=float, default=1e-5, help="Swiss-Prot E-value threshold")
        command.add_argument("--phrogs", help="Optional prepared PHROGs MMseqs2 profile database prefix")
        command.add_argument("--phrogs-annotations", help="Optional PHROGs annotation table path")
        command.add_argument("--phrogs-hmm", help="Optional existing all_phrogs.h3m profile database for sensitive PyHMMER search")
        command.add_argument("--mmseqs", help="Optional MMseqs2 executable path for PHROGs")
        command.add_argument("--phrogs-evalue", type=float, default=1e-5, help="Manual PHROGs E-value threshold")
        command.add_argument("--phrogs-coverage", type=float, default=0.5, help="Manual PHROGs query coverage threshold (0-1)")
        command.add_argument("--phrogs-score", type=float, help="Optional manual PHROGs MMseqs2 bit-score threshold")
        command.add_argument("--phrogs-identity", type=float, help="Optional manual PHROGs identity threshold (MMseqs2 fident percent)")
        command.add_argument("--phrogs-alignment-length", type=int, help="Optional manual PHROGs alignment-length threshold")
        command.add_argument("--quiet", action="store_true", help="Suppress progress display")
        command.add_argument("--no-progress", action="store_true", help="Disable dynamic progress rendering")
        command.add_argument("--threads", type=int, default=1, help="Threads for evidence tools")
        command.add_argument("--mock-evidence", action="store_true", help="Use demonstration evidence; fixture/testing only")
    revise = subcommands.add_parser("revise", help="Revise a completed GenBank submission without rerunning analysis")
    revise.add_argument("results_dir")
    revise.add_argument("--corrections", required=True)
    revise.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        return 0
    if args.command == "gui":
        from .gui.launcher import main as gui_main
        return gui_main([])
    if args.command == "doctor":
        from .preflight import doctor, doctor_text, write_doctor_report
        payload = write_doctor_report(args.output, deep=args.deep) if args.output else doctor(deep=args.deep)
        print(json.dumps(payload, indent=2, sort_keys=True) if args.json or args.output else doctor_text(payload))
        required = {"phanotate.py", "hmmscan", "mmseqs", "diamond", "mash"}
        executables_ready = all(x.get("status") == "READY" for x in payload["executables"]
                                if x.get("name") in required)
        deep_ready = all(item.get("status") == "READY" for item in payload.get("deep_checks", []))
        return 0 if executables_ready and deep_ready else 1
    if args.command == "families":
        from .family import build_database, assign_protein_family, write_assignments, _fasta
        if args.families_command == 'build':
            metadata=json.loads(Path(args.metadata).read_text()) if args.metadata else None
            build_database(args.proteins,args.output,metadata=metadata,registry=args.registry,database_version=args.version,backend='MMSEQS2' if args.backend=='mmseqs' else 'EXACT_SEQUENCE',mmseqs=args.mmseqs,config={'family_build':{'minimum_identity':args.minimum_identity,'minimum_coverage':args.minimum_coverage,'coverage_mode':args.coverage_mode,'clustering_mode':args.clustering_mode}})
        elif args.families_command == 'compare':
            from .family_compare import compare_database
            core=args.core_genomes.split(',') if args.core_genomes else None
            compare_database(args.database,args.output,core_genomes=core)
        elif args.families_command == 'enrich':
            from .family_enrichment import enrich_database
            enrich_database(args.database,args.results,args.output)
        elif args.families_command == 'prioritize':
            from .family_priority import prioritize
            prioritize(args.database,args.output)
        elif args.families_command == 'validate-external':
            from .family_external import validate_external
            if bool(args.reference_db)==bool(args.reference_proteins): parser.error('provide exactly one of --reference-db or --reference-proteins')
            validate_external(args.database,args.reference_proteins,args.output,args.reference_metadata,args.mmseqs,args.top_n,args.reference_db)
        elif args.families_command == 'sensitivity':
            from .family import build_database
            from .family_compare import compare_database
            import tempfile
            metadata=json.loads(Path(args.metadata).read_text()) if args.metadata else None
            out=Path(args.output); out.mkdir(parents=True,exist_ok=True); summaries={}
            for value in sorted({float(x) for x in args.identities.split(',')}):
                db=out/f'identity_{value:.3f}'
                build_database(args.proteins,db,metadata=metadata,database_version=f'identity-{value}',backend='MMSEQS2' if args.backend=='mmseqs' else 'EXACT_SEQUENCE',mmseqs=args.mmseqs,config={'family_build':{'minimum_identity':value,'minimum_coverage':args.minimum_coverage,'coverage_mode':args.coverage_mode,'clustering_mode':args.clustering_mode}})
                summaries[value]=compare_database(db,db,core_genomes=None)
                summaries[value].update({'minimum_coverage':args.minimum_coverage,'coverage_mode':args.coverage_mode,'clustering_mode':args.clustering_mode,'backend':'MMSEQS2' if args.backend=='mmseqs' else 'EXACT_SEQUENCE'})
            from .family_compare import sensitivity_rows
            rows=sensitivity_rows(summaries)
            with (out/'pmf_threshold_sensitivity.tsv').open('w') as h:
                if rows:
                    h.write('\t'.join(rows[0])+'\n'); [h.write('\t'.join(str(r[k]) for k in rows[0])+'\n') for r in rows]
            (out/'pmf_threshold_sensitivity.json').write_text(json.dumps(rows,indent=2,sort_keys=True))
        else:
            from .family import assign_protein_families
            write_assignments(assign_protein_families([{'protein_id':pid,'sequence':seq} for pid,seq in _fasta(args.proteins).items()],args.database,mmseqs=args.mmseqs),args.output)
        print(f"Protein-family workflow complete. Outputs: {args.output}")
        return 0
    if args.command == "benchmark":
        from .benchmark import import_multiphate, import_phold, import_prokka, import_pharokka, import_phagemine, benchmark
        methods={}
        if args.phagemine_results: methods["PHAGEMINE"]=import_phagemine(args.phagemine_results)
        if args.prokka_gff: methods["Prokka"]=import_prokka(args.prokka_gff)
        if args.pharokka_gff: methods["Pharokka"]=import_pharokka(args.pharokka_gff)
        if args.multiphate_gff: methods["multiPhATE2"]=import_multiphate(args.multiphate_gff)
        if args.phold_genbank: methods["Phold"]=import_phold(args.phold_genbank)
        truth_records = None
        if args.truth_gff: truth_records = import_pharokka(args.truth_gff, genome_id="reviewed_truth")
        if args.truth_genbank: truth_records = import_phold(args.truth_genbank, genome_id="reviewed_truth")
        references = {name: truth_records for name in methods} if truth_records else None
        benchmark(methods,args.output,references=references)
        print(f"PhageMine benchmark complete. Outputs: {args.output}"); return 0
    if args.command == "revise":
        from .submission_revision import revise as revise_submission
        try:
            result = revise_submission(args.results_dir, args.corrections, args.output)
        except (OSError, ValueError, RuntimeError) as exc:
            parser.error(str(exc))
        print(json.dumps(result, indent=2, sort_keys=True))
        print(f"Revised GenBank package: {Path(args.output) / 'submission_v2' / 'genbank_submission'}")
        return 0
    if args.command == "resume-stage":
        from .stage_resume import resume_stage
        try: print(json.dumps(resume_stage(args.source,args.swissprot,args.diamond),indent=2,sort_keys=True))
        except (OSError,ValueError,RuntimeError) as exc: parser.error(str(exc))
        return 0
    if args.command == "compare":
        try: compare(args.results, args.output, args.mmseqs)
        except (OSError, ValueError, RuntimeError) as exc: parser.error(str(exc))
        print(f"PhageMine comparison complete. Outputs: {args.output}")
        return 0
    if args.command == "batch":
        try: batch(args.input_dir,args.output,args.recursive,args.resume_existing,args.fail_fast,args.gene_predictor,args.phanotate,ProgressReporter(quiet=False),args.reconcile_orfs,args.prodigal,args.threads,args.evidence_profile,args.mode, gene_model_policy=args.gene_model_policy, gene_model_profile=args.gene_model_profile, molecule_type=args.molecule_type, segmented=args.segmented)
        except (OSError, ValueError, RuntimeError) as exc: parser.error(str(exc))
        print(f"PhageMine batch complete. Outputs: {args.output}"); return 0
    if args.command == "extract":
        if args.extract_command != "protein":
            parser.error("unsupported extraction target")
        from .reporting import extract_protein_record
        try:
            payload = extract_protein_record(args.run, args.protein_id)
        except (OSError, ValueError, KeyError) as exc:
            parser.error(str(exc))
        requested = args.protein_fasta or args.cds_fasta or args.evidence
        if not requested:
            print(json.dumps(payload, indent=2, sort_keys=True))
        else:
            if args.protein_fasta: print(payload["protein_fasta"], end="")
            if args.cds_fasta: print(payload["cds_fasta"], end="")
            if args.evidence: print(json.dumps(payload["evidence"], indent=2, sort_keys=True))
        return 0
    if args.command == "resume":
        try:
            count = resume(args.source, args.output, args.run_missing_evidence, args.refresh_evidence,
                           args.mmseqs, args.phrogs, args.phrogs_annotations, args.phrogs_hmm,
                           args.phrogs_evalue, args.phrogs_coverage, args.phrogs_score,
                           ProgressReporter(quiet=False), threads=args.threads)
        except (OSError, ValueError, RuntimeError) as exc:
            parser.error(str(exc))
        print(f"PhageMine resume complete: {count} predicted proteins. Outputs: {args.output}")
        return 0
    if args.command == "reclassify":
        try:
            count = reclassify(args.source, args.output, ProgressReporter(quiet=False))
        except (OSError, ValueError, RuntimeError) as exc:
            parser.error(str(exc))
        print(f"PhageMine reclassification complete: {count} predicted proteins. Outputs: {args.output}")
        print("Database searches run: none; persisted evidence was reused.")
        return 0
    if args.command == "databases":
        if args.database_command == "install":
            from .database_installer import (
                DatabaseInstallError, DatabaseInstaller, RESOURCE_ORDER,
                installation_plan, results_json,
            )
            if bool(args.resource) == bool(args.all):
                parser.error("choose one database name or --all")
            selected = RESOURCE_ORDER if args.all else (args.resource,)
            print(installation_plan(selected), file=sys.stderr)
            if args.dry_run:
                return 0
            try:
                results = DatabaseInstaller(
                    args.directory,
                    force=args.force,
                    keep_downloads=args.keep_downloads,
                    threads=args.threads,
                ).install_many(selected)
            except (OSError, DatabaseInstallError) as exc:
                parser.error(str(exc))
            if args.json:
                print(results_json(results))
            else:
                for result in results:
                    action = "installed" if result.installed else "already ready"
                    print(f"{result.resource}: {result.status} ({action}) - {result.path}")
                print("\nRun 'phagemine doctor' to verify the complete installation.")
            return 0
        manager = EvidenceResourceManager()
        if args.database_command == "register":
            name = args.name or args.resource_type.upper()
            provenance = {key: value for key, value in (("annotations_path", args.annotations),
                          ("metadata_path", args.metadata), ("hmm_profiles_path", args.hmm_profiles)) if value}
            resource = manager.register(name, args.resource_type, args.path, version=args.version, checksum=args.checksum, notes=args.notes, provenance=provenance)
            print(json.dumps(resource.metadata(), indent=2, sort_keys=True))
            return 0
        if args.database_command == "status":
            print(json.dumps(manager.list(), indent=2, sort_keys=True))
            return 0
        if args.database_command == "attach-phrogs-hmm":
            hmm_path = Path(args.path).expanduser().resolve()
            if not hmm_path.is_file():
                parser.error(f"PHROGs HMM profile database does not exist: {hmm_path}")
            current = manager.find(ResourceType.PHROGS)
            if current is None:
                parser.error("No PHROGs resource is registered; install or register PHROGs first")
            provenance = {**(current.get("provenance") or {}), "hmm_profiles_path": str(hmm_path)}
            manager.register(
                current["name"], ResourceType.PHROGS, current["path"],
                version=current.get("version"), checksum=current.get("checksum"),
                required_tools=current.get("required_tools"),
                preparation_status=current.get("preparation_status", "prepared"),
                notes=current.get("notes"), provenance=provenance)
            checked = manager.validate(current["name"])
            if checked is None or checked.get("status") != "READY":
                parser.error("PHROGs resource failed validation after attaching HMM profiles: " +
                             "; ".join((checked or {}).get("validation_errors", [])))
            print(f"Attached PHROGs profile HMM database: {hmm_path}")
            return 0
        removed = manager.unregister(args.name)
        print(f"Removed {args.name}" if removed else f"No registration found for {args.name}")
        return 0 if removed else 1
    output = args.output or str(Path.cwd() / f"{Path(args.fasta).stem}_phagemine_results")
    try:
        from .preflight import preflight_resources
        preflight_resources({
            "PFAM": {"path": args.pfam, "required_tools": [args.pfam_hmmscan or "hmmscan"]} if args.pfam else None,
            "VOGDB": {"path": args.vogdb, "required_tools": [args.vog_hmmscan or "hmmscan"], "provenance": {"annotations_path": args.vog_annotations}} if args.vogdb else None,
            "SWISSPROT": {"path": args.swissprot, "required_tools": [args.diamond or "diamond"], "provenance": {"metadata_path": args.swissprot_metadata}} if args.swissprot else None,
            "PHROGS": {"path": args.phrogs, "required_tools": [args.mmseqs or "mmseqs"], "provenance": {"annotations_path": args.phrogs_annotations, "hmm_profiles_path": args.phrogs_hmm}} if args.phrogs else None,
        })
        metadata = SubmissionMetadata.from_dict(json.loads(Path(args.metadata).read_text())) if args.metadata else None
        from .sequencing_provenance import SequencingProvenance
        sequencing_provenance = SequencingProvenance.from_dict(json.loads(Path(args.sequencing_provenance).read_text())) if args.sequencing_provenance else SequencingProvenance()
        from .gene_prediction import create_predictor
        # A complete single-genome ``run`` performs the secondary Prodigal
        # comparison automatically.  Other workflows retain the explicit
        # switch because they may intentionally be annotation-only.
        reconcile_orfs = args.reconcile_orfs or args.command == "run"
        count = run(args.fasta, output, args.command, metadata, args.table2asn, create_predictor(args.gene_predictor, args.phanotate), sequencing_provenance=sequencing_provenance, pfam_path=args.pfam, pfam_hmmscan=args.pfam_hmmscan, pfam_evalue=args.pfam_evalue, pfam_coverage=args.pfam_coverage, pfam_trusted_cutoff=args.pfam_trusted_cutoff, use_mock_evidence=args.mock_evidence, pfam_threshold_mode=args.pfam_threshold_mode, vog_path=args.vogdb, vog_annotations=args.vog_annotations, vog_hmmscan=args.vog_hmmscan, vog_evalue=args.vog_evalue, vog_coverage=args.vog_coverage, swissprot_path=args.swissprot, swissprot_metadata=args.swissprot_metadata, diamond=args.diamond, swissprot_evalue=args.swissprot_evalue, phrogs_path=args.phrogs, phrogs_annotations=args.phrogs_annotations, phrogs_hmm_path=args.phrogs_hmm, mmseqs=args.mmseqs, phrogs_evalue=args.phrogs_evalue, phrogs_coverage=args.phrogs_coverage, phrogs_score=args.phrogs_score, phrogs_identity=args.phrogs_identity, phrogs_alignment_length=args.phrogs_alignment_length, reconcile_orfs=reconcile_orfs, prodigal=args.prodigal, threads=args.threads, progress=ProgressReporter(quiet=args.quiet, no_progress=args.no_progress), gene_model_policy=args.gene_model_policy, gene_model_profile=args.gene_model_profile, molecule_type=args.molecule_type, segmented=args.segmented)
    except (OSError, ValueError, RuntimeError) as exc:
        parser.error(str(exc))
    print(f"PhageMine complete: {count} predicted proteins. Outputs: {output}")
    print(f"GenBank pre-submission package: {Path(output) / 'genbank_submission'}")
    if args.mock_evidence:
        print("WARNING: mock evidence was enabled; it is demonstration-only and not biological evidence.")
    else:
        print("Functional evidence is limited to explicitly available adapters; unavailable evidence was not fabricated.")
    return 0
