#!/usr/bin/env python3
"""
crispr_cas_offtarget_predictor.cli
==================================
Author: Lakshitha Dananjaya (devsoniclk)
Co-authored-by: OpenBot <335604065+openbotgit[bot]@users.noreply.github.com>

Rich-enhanced command-line interface for CRISPR-Cas off-target prediction.
"""

import argparse
import json
import sys

from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from .engine import CRISPRPredictorOrchestrator
from .server import run_stdio_server

console = Console()


def main():
    parser = argparse.ArgumentParser(
        description="CRISPR-Cas Guide RNA Off-Target Cleavage Predictor (SpCas9 / Cas12a)"
    )
    subparsers = parser.add_subparsers(dest="command", help="Command to execute")

    # Predict command
    pred_parser = subparsers.add_parser("predict", help="Predict off-target cleavage probability")
    pred_parser.add_argument(
        "--guide", "-g", required=True, type=str, help="Guide RNA protospacer sequence (20nt)"
    )
    pred_parser.add_argument(
        "--target", "-t", required=True, type=str, help="Target genomic protospacer sequence"
    )
    pred_parser.add_argument(
        "--pam", "-p", required=True, type=str, help="PAM sequence (e.g., TGG or TTTA)"
    )
    pred_parser.add_argument(
        "--atac",
        "-a",
        type=float,
        default=None,
        help="ATAC-seq chromatin accessibility peak signal",
    )
    pred_parser.add_argument(
        "--cas",
        "-c",
        type=str,
        default="SpCas9",
        choices=["SpCas9", "Cas12a"],
        help="Endonuclease system",
    )
    pred_parser.add_argument(
        "--json", action="store_true", help="Output results in structured JSON format"
    )

    # Audit command
    audit_parser = subparsers.add_parser("audit", help="Audit codebase for secrets and provenance")
    audit_parser.add_argument("file", type=str, help="Path to source file to audit")
    audit_parser.add_argument(
        "--json", action="store_true", help="Output results in structured JSON format"
    )

    # Serve command
    subparsers.add_parser("serve", help="Run FastMCP JSON-RPC stdio server")

    # Legacy direct flags for backward compatibility
    parser.add_argument("--guide", type=str, help=argparse.SUPPRESS)
    parser.add_argument("--target", type=str, help=argparse.SUPPRESS)
    parser.add_argument("--pam", type=str, help=argparse.SUPPRESS)
    parser.add_argument("--atac", type=float, default=None, help=argparse.SUPPRESS)
    parser.add_argument(
        "--cas", type=str, default="SpCas9", choices=["SpCas9", "Cas12a"], help=argparse.SUPPRESS
    )
    parser.add_argument("--audit-file", type=str, help=argparse.SUPPRESS)
    parser.add_argument("--json", action="store_true", help=argparse.SUPPRESS)

    args = parser.parse_args()
    orchestrator = CRISPRPredictorOrchestrator()

    if args.command == "serve":
        run_stdio_server()
        return

    # Handle legacy flags or subcommands
    audit_target = args.file if args.command == "audit" else getattr(args, "audit_file", None)
    if audit_target:
        res = orchestrator.run_full_suite(audit_target)
        if args.json:
            print(json.dumps(res, indent=2))
        else:
            table = Table(title="Section 10 Codebase Security & Provenance Audit")
            table.add_column("Property", style="bold cyan")
            table.add_column("Value")
            table.add_row("Target File", res.get("target", ""))
            table.add_row(
                "Release Veto",
                "[red]YES (BLOCKED)[/red]"
                if res.get("release_veto")
                else "[green]NO (CLEARED)[/green]",
            )
            table.add_row("Secrets Found", str(len(res.get("secret_findings", []))))
            table.add_row(
                "Sovereign Provenance",
                "[green]VERIFIED (devsoniclk)[/green]"
                if res.get("provenance", {}).get("is_compliant")
                else "[red]FAIL[/red]",
            )
            console.print(table)
            if res.get("veto_reasons"):
                for r in res["veto_reasons"]:
                    console.print(f"[bold red]VETO REASON:[/bold red] {r}")
        sys.exit(1 if res.get("release_veto") else 0)

    guide = args.guide if (args.command == "predict" or args.guide) else None
    target = args.target if (args.command == "predict" or args.target) else None
    pam = args.pam if (args.command == "predict" or args.pam) else None

    if guide and target and pam:
        res = orchestrator.predict_offtarget(guide, target, pam, args.atac, args.cas)
        if args.json:
            print(json.dumps(res, indent=2))
        else:
            panel = Panel.fit(
                f"[bold cyan]Endonuclease:[/bold cyan] {res.get('endonuclease', args.cas)}\n"
                f"[bold yellow]CFD Cleavage Score:[/bold yellow] {res['cfd_score']:.6f}\n"
                f"[bold green]Effective Cleavage Risk:[/bold green] {res['effective_cleavage_score']:.6f}\n"
                f"[bold magenta]Active Off-Target:[/bold magenta] {'[bold red]YES (HIGH RISK)[/bold red]' if res.get('is_active_offtarget', False) else '[green]NO (TOLERABLE)[/green]'}\n"
                f"[bold blue]Chromatin State:[/bold blue] {res['chromatin_state']}\n"
                f"[bold white]PAM Score:[/bold white] {res.get('pam_score', 0.0):.3f}",
                title="CRISPR Cleavage Prediction Report",
                border_style="cyan",
            )
            console.print(panel)
        sys.exit(0)

    parser.print_help()


if __name__ == "__main__":
    main()
