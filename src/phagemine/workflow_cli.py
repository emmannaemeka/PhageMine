"""CLI entry points for annotation review, feature imports and reproducibility."""
from __future__ import annotations
import argparse
import json
from pathlib import Path

COMMANDS = {"curate", "rna", "evidence-import", "bundle", "validate-run", "profile", "database-health", "benchmark-curated"}


def add_commands(subcommands) -> None:
    health = subcommands.add_parser("database-health", help="Offline database readiness, version and freshness checks")
    health.add_argument("--expected-versions", help="JSON mapping resource names/types to required versions")
    health.add_argument("--max-age-days", type=int, default=365)
    health.add_argument("--checksums", action="store_true")
    health.add_argument("--output")
    benchmark = subcommands.add_parser("benchmark-curated", help="Score TSV predictions against provenance-bearing curated references")
    for name in ("predictions", "references", "output"):
        benchmark.add_argument("--" + name, required=True)
    benchmark.add_argument("--synonyms"); benchmark.add_argument("--adjudications")
    curate = subcommands.add_parser("curate", help="Apply audited product/gene/note edits")
    actions = curate.add_subparsers(dest="action", required=True)
    for action in ("template", "apply"):
        command = actions.add_parser(action)
        command.add_argument("--run", required=True); command.add_argument("--output", required=True)
        if action == "apply": command.add_argument("--changes", required=True)
    rna = subcommands.add_parser("rna", help="Run tRNAscan-SE or import noncoding-RNA GFF3")
    actions = rna.add_subparsers(dest="action", required=True)
    for action in ("scan", "import"):
        command = actions.add_parser(action)
        command.add_argument("--run", required=True); command.add_argument("--output", required=True)
        if action == "scan":
            command.add_argument("--executable", default="tRNAscan-SE"); command.add_argument("--threads", type=int, default=1)
        else:
            command.add_argument("--gff", required=True); command.add_argument("--source", required=True); command.add_argument("--version", required=True)
    evidence = subcommands.add_parser("evidence-import", help="Import sequence-validated external annotation proposals")
    for name in ("run", "genbank", "version", "database-version", "output"):
        evidence.add_argument("--" + name, required=True)
    evidence.add_argument("--source", default="Phold")
    bundle = subcommands.add_parser("bundle", help="Export or verify a streaming reproducibility ZIP")
    actions = bundle.add_subparsers(dest="action", required=True)
    create = actions.add_parser("create"); create.add_argument("--run", required=True); create.add_argument("--output", required=True)
    verify = actions.add_parser("verify"); verify.add_argument("archive")
    validate = subcommands.add_parser("validate-run", help="Check result integrity without rerunning searches")
    validate.add_argument("--run", required=True); validate.add_argument("--output")
    profile = subcommands.add_parser("profile", help="Measure elapsed time and sampled process-tree memory")
    profile.add_argument("--output", required=True); profile.add_argument("argv", nargs=argparse.REMAINDER)


def dispatch(args) -> int:
    if args.command == "database-health":
        from .database_health import health
        expected = json.loads(Path(args.expected_versions).read_text()) if args.expected_versions else None
        result = health(expected_versions=expected, max_age_days=args.max_age_days, check_checksums=args.checksums)
        if args.output: Path(args.output).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        print(json.dumps(result, indent=2, sort_keys=True)); return 1 if result["status"] == "NEEDS_ATTENTION" else 0
    elif args.command == "benchmark-curated":
        from .curated_benchmark import evaluate
        result = evaluate(args.predictions, args.references, args.output, synonyms=args.synonyms, adjudications=args.adjudications)
    elif args.command == "curate":
        from .curation import apply, template
        if args.action == "template":
            output = Path(args.output).expanduser()
            if output.exists(): raise ValueError("Curation template already exists")
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(json.dumps(template(args.run), indent=2, sort_keys=True) + "\n")
            result = {"status": "TEMPLATE_CREATED", "output": str(output)}
        else: result = apply(args.run, args.changes, args.output)
    elif args.command == "rna":
        from .rna_features import import_run, scan_run
        result = scan_run(args.run, args.output, args.executable, args.threads) if args.action == "scan" else import_run(args.run, args.gff, args.output, args.source, args.version)
    elif args.command == "evidence-import":
        from .external_annotation import import_genbank
        result = import_genbank(args.run, args.genbank, args.output, source=args.source, version=args.version, database_version=args.database_version)
    elif args.command == "bundle":
        from .reproducibility import bundle, verify
        result = bundle(args.run, args.output) if args.action == "create" else verify(args.archive)
    elif args.command == "profile":
        from .operational_validation import profile
        command = args.argv[1:] if args.argv and args.argv[0] == "--" else args.argv
        result = profile(command, args.output)
        print(json.dumps(result, indent=2, sort_keys=True)); return result["exit_code"]
    else:
        from .operational_validation import validate_run
        result = validate_run(args.run)
        if args.output: Path(args.output).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True)); return 0
