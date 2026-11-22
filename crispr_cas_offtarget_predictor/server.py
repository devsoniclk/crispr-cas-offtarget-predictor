#!/usr/bin/env python3
"""
crispr_cas_offtarget_predictor.server
=====================================
Author: Lakshitha Dananjaya (devsoniclk)
Co-authored-by: OpenBot <335604065+openbotgit[bot]@users.noreply.github.com>

JSON-RPC 2.0 stdio FastMCP Server for CRISPR-Cas Guide RNA Off-Target Cleavage Analysis.
Exposes modular tools for autonomous AI agent swarms and bioinformatic pipelines.
"""

import json
import sys
from typing import Any

from .engine import (
    ChromatinAccessibilityEngine,
    CRISPRBulgeEngine,
    CRISPRPredictorOrchestrator,
    GuideSpecificityAggregator,
)


class CRISPROffTargetMCPServer:
    """JSON-RPC 2.0 stdio tool execution server for CRISPR off-target analysis."""

    def __init__(self):
        self.orchestrator = CRISPRPredictorOrchestrator()

    def predict_offtarget(
        self,
        guide: str,
        target: str,
        pam: str,
        atac_signal: float | None = None,
        endonuclease: str = "SpCas9",
    ) -> dict[str, Any]:
        """
        Predicts intrinsic and epigenetically-modulated off-target cleavage probability
        for a guide RNA against a target DNA protospacer and PAM.
        Supports SpCas9 (3' NGG/NAG/NGA) and Cas12a/Cpf1 (5' TTTV/TTTN).
        """
        return self.orchestrator.predict_offtarget(
            guide=guide,
            target=target,
            pam=pam,
            atac_signal=atac_signal,
            endonuclease=endonuclease,
        )

    def score_bulge(
        self,
        guide: str,
        target: str,
        pam: str,
        endonuclease: str = "SpCas9",
    ) -> dict[str, Any]:
        """
        Evaluates indel differences (RNA or DNA bulges of 1-2 nt) between guide RNA
        and target DNA, computing position-weighted penalty scores.
        """
        return CRISPRBulgeEngine.score_bulge(
            guide=guide,
            target=target,
            pam=pam,
            endonuclease=endonuclease,
        )

    def modulate_chromatin(
        self,
        intrinsic_cfd: float,
        atac_signal: float,
        kd: float = 10.0,
    ) -> dict[str, Any]:
        """
        Modulates intrinsic CFD score with cell-type specific chromatin accessibility
        (ATAC-seq/DNase-seq signal), accounting for the heterochromatin cleavage floor.
        """
        return ChromatinAccessibilityEngine.compute_effective_cleavage_risk(
            intrinsic_cfd=intrinsic_cfd,
            atac_signal=atac_signal,
            kd=kd,
        )

    def compute_guide_specificity(
        self,
        off_target_scores: list[float],
    ) -> dict[str, Any]:
        """
        Computes global guide RNA specificity score [0-100] across a set of candidate
        off-target cleavage sites.
        """
        score = GuideSpecificityAggregator.compute_specificity_score(off_target_scores)
        return {
            "specificity_score": score,
            "candidate_count": len(off_target_scores),
            "total_offtarget_mass": sum(off_target_scores),
            "tier": "HIGH_SPECIFICITY"
            if score >= 80.0
            else ("MODERATE" if score >= 50.0 else "LOW_SPECIFICITY"),
        }

    def audit_codebase(
        self,
        file_path: str,
    ) -> dict[str, Any]:
        """
        Audits a target source file for secrets, copyleft licenses, and sovereign authorship.
        """
        return self.orchestrator.run_full_suite(file_path)

    def handle_request(self, req: dict[str, Any]) -> dict[str, Any]:
        """Processes a single JSON-RPC 2.0 request dictionary."""
        req_id = req.get("id")
        method = req.get("method")
        params = req.get("params", {})

        if method == "predict_offtarget":
            res = self.predict_offtarget(
                guide=params.get("guide", ""),
                target=params.get("target", ""),
                pam=params.get("pam", ""),
                atac_signal=params.get("atac_signal"),
                endonuclease=params.get("endonuclease", "SpCas9"),
            )
        elif method == "score_bulge":
            res = self.score_bulge(
                guide=params.get("guide", ""),
                target=params.get("target", ""),
                pam=params.get("pam", ""),
                endonuclease=params.get("endonuclease", "SpCas9"),
            )
        elif method == "modulate_chromatin":
            res = self.modulate_chromatin(
                intrinsic_cfd=float(params.get("intrinsic_cfd", 0.0)),
                atac_signal=float(params.get("atac_signal", 0.0)),
                kd=float(params.get("kd", 10.0)),
            )
        elif method == "compute_guide_specificity":
            res = self.compute_guide_specificity(
                off_target_scores=params.get("off_target_scores", []),
            )
        elif method == "audit_codebase":
            res = self.audit_codebase(
                file_path=params.get("file_path", ""),
            )
        elif method == "tools/list":
            res = {
                "tools": [
                    {
                        "name": "predict_offtarget",
                        "description": "Predict SpCas9/Cas12a off-target cleavage with CFD and epigenetics",
                    },
                    {
                        "name": "score_bulge",
                        "description": "Score 1-2bp RNA/DNA indel bulges in guide-target alignments",
                    },
                    {
                        "name": "modulate_chromatin",
                        "description": "Modulate intrinsic CFD with ATAC-seq accessibility factor",
                    },
                    {
                        "name": "compute_guide_specificity",
                        "description": "Compute whole-guide specificity score [0-100]",
                    },
                    {
                        "name": "audit_codebase",
                        "description": "Audit source file for secrets and sovereign provenance",
                    },
                ]
            }
        else:
            return {
                "jsonrpc": "2.0",
                "id": req_id,
                "error": {"code": -32601, "message": f"Method {method} not found"},
            }

        return {"jsonrpc": "2.0", "id": req_id, "result": res}


def run_stdio_server():
    """Stdio JSON-RPC loop handling tool execution."""
    server = CRISPROffTargetMCPServer()
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        try:
            req = json.loads(line)
            resp = server.handle_request(req)
            sys.stdout.write(json.dumps(resp) + "\n")
            sys.stdout.flush()
        except Exception as e:  # noqa: BLE001
            err_resp = {"jsonrpc": "2.0", "error": {"code": -32603, "message": str(e)}}
            sys.stdout.write(json.dumps(err_resp) + "\n")
            sys.stdout.flush()


if __name__ == "__main__":
    run_stdio_server()
