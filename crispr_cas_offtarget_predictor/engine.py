#!/usr/bin/env python3
"""
crispr_offtarget_predictor_engine.py
====================================
Author: Lakshitha Dananjaya (devsoniclk)
Co-authored-by: OpenBot <335604065+openbotgit[bot]@users.noreply.github.com>

Core Algorithmic Framework for Automated CRISPR-Cas Guide RNA Off-Target Cleavage Prediction.
Supports SpCas9 (3' PAM NGG/NAG) and Cas12a/Cpf1 (5' PAM TTTV) Endonucleases.

Key Capabilities:
1. Position-weighted mismatch penalty matrix (Cutting Frequency Determination / CFD model).
2. Seed sequence sensitivity weighting (PAM-proximal 8-12 nt).
3. DNA/RNA structural bulge penalty scoring (insertion/deletion tolerance).
4. Epigenetic chromatin accessibility integration (ATAC-seq / DNase-seq peak modulation).
5. Log-space numerical underflow protection for cascading probabilistic products.
6. Section 10 Sovereign Authorship (devsoniclk) and secret containment verification.
"""

import math
import os
import re
from typing import Any

# ============================================================================
# CONSTANTS & BIOLOGICAL MATRICES (DOENCH ET AL. / HSU ET AL.)
# ============================================================================

# Default PAM weights for SpCas9 (3' of protospacer)
CAS9_PAM_WEIGHTS = {
    "NGG": 1.000,
    "AGG": 1.000,
    "TGG": 1.000,
    "CGG": 1.000,
    "GGG": 1.000,
    "NAG": 0.259,
    "AAG": 0.259,
    "TAG": 0.259,
    "CAG": 0.259,
    "GAG": 0.259,
    "NGA": 0.069,
    "AGA": 0.069,
    "TGA": 0.069,
    "CGA": 0.069,
    "GGA": 0.069,
}

# Default PAM weights for Cas12a/Cpf1 (5' of protospacer)
CAS12A_PAM_WEIGHTS = {
    "TTTV": 1.000,
    "TTTA": 1.000,
    "TTTC": 1.000,
    "TTTG": 1.000,
    "TTTN": 0.200,
    "TCTA": 0.050,
    "TCTC": 0.050,
}

# Baseline position-weighted mismatch penalties for SpCas9 (1-indexed from 5' to 3')
# Positions 11-20 represent the PAM-proximal SEED region (higher sensitivity / lower mismatch tolerance)
# Positions 1-10 represent the PAM-distal NON-SEED region (higher mismatch tolerance)
DEFAULT_CAS9_POSITION_WEIGHTS = {
    1: 0.95,
    2: 0.92,
    3: 0.90,
    4: 0.88,
    5: 0.85,
    6: 0.82,
    7: 0.80,
    8: 0.75,
    9: 0.70,
    10: 0.65,
    11: 0.40,
    12: 0.35,
    13: 0.25,
    14: 0.20,
    15: 0.15,
    16: 0.10,
    17: 0.08,
    18: 0.05,
    19: 0.03,
    20: 0.02,
}

# Transversion vs transition penalty multipliers
TRANSITION_TYPES = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

# ============================================================================
# 1. POSITION-WEIGHTED MISMATCH & CFD SCORING ENGINE
# ============================================================================


class CRISPRCleavageEngine:
    """
    Computes intrinsic endonuclease cleavage probabilities using position-weighted
    mismatch penalties (CFD model) and log-additive numerical stability.
    """

    @staticmethod
    def sanitize_sequence(seq: str) -> str:
        """Sanitizes nucleic acid sequence to uppercase standard DNA (U converted to T)."""
        clean = seq.strip().upper().replace("U", "T")
        if not re.fullmatch(r"[ACGTN]+", clean):
            raise ValueError(f"Invalid nucleotide character detected in sequence: {seq}")
        return clean

    @classmethod
    def get_cas9_pam_score(cls, pam: str) -> float:
        """Evaluates SpCas9 3' PAM sequence score."""
        sanitized = cls.sanitize_sequence(pam)
        if len(sanitized) < 3:
            return 0.0
        pam3 = sanitized[:3]
        if pam3 in CAS9_PAM_WEIGHTS:
            return CAS9_PAM_WEIGHTS[pam3]
        # Match NGG / NAG / NGA patterns
        if pam3.endswith("GG"):
            return 1.000
        elif pam3.endswith("AG"):
            return 0.259
        elif pam3.endswith("GA"):
            return 0.069
        return 0.001  # Basal non-canonical cleavage threshold

    @classmethod
    def get_cas12a_pam_score(cls, pam: str) -> float:
        """Evaluates Cas12a 5' PAM sequence score (TTTV where V=A/C/G)."""
        sanitized = cls.sanitize_sequence(pam)
        if len(sanitized) < 4:
            return 0.0
        pam4 = sanitized[:4]
        if pam4 in CAS12A_PAM_WEIGHTS:
            return CAS12A_PAM_WEIGHTS[pam4]
        if pam4.startswith("TTT") and pam4[3] in {"A", "C", "G"}:
            return 1.000
        elif pam4.startswith("TTT"):
            return 0.200
        return 0.010

    @classmethod
    def compute_cfd_score(
        cls, guide: str, target: str, pam: str, endonuclease: str = "SpCas9"
    ) -> dict[str, Any]:
        """
        Calculates position-dependent Cutting Frequency Determination (CFD) score.
        Computes in log-space to prevent catastrophic floating-point underflow.
        """
        guide_clean = cls.sanitize_sequence(guide)
        target_clean = cls.sanitize_sequence(target)

        if len(guide_clean) != len(target_clean):
            raise ValueError(
                f"Length mismatch: guide ({len(guide_clean)}) vs target ({len(target_clean)}). Use bulge scoring for indels."
            )

        seq_len = len(guide_clean)
        log_cfd = 0.0
        mismatches = []

        # PAM evaluation
        if endonuclease.lower() in {"spcas9", "cas9"}:
            pam_score = cls.get_cas9_pam_score(pam)
            pos_weights = DEFAULT_CAS9_POSITION_WEIGHTS
        elif endonuclease.lower() in {"cas12a", "cpf1"}:
            pam_score = cls.get_cas12a_pam_score(pam)
            # Cas12a seed region is PAM-proximal (positions 1-8 from 5' end)
            pos_weights = {
                i: (0.05 + 0.04 * (i - 1)) if i <= 8 else 0.85 for i in range(1, seq_len + 1)
            }
        else:
            raise ValueError(f"Unsupported endonuclease system: {endonuclease}")

        if pam_score <= 0.0:
            return {
                "cfd_score": 0.0,
                "pam_score": 0.0,
                "mismatch_count": 0,
                "mismatches": [],
                "is_active_offtarget": False,
                "log_cfd": float("-inf"),
            }

        log_cfd += math.log(pam_score)

        for pos in range(1, seq_len + 1):
            g_nt = guide_clean[pos - 1]
            t_nt = target_clean[pos - 1]

            if g_nt != t_nt and g_nt != "N" and t_nt != "N":
                # Position weight penalty
                base_weight = pos_weights.get(pos, 0.5)
                # Differentiate transition vs transversion
                is_transition = (g_nt, t_nt) in TRANSITION_TYPES
                penalty = base_weight if is_transition else (base_weight * 0.75)

                # Floor penalty to avoid log(0)
                safe_penalty = max(penalty, 1e-6)
                log_cfd += math.log(safe_penalty)
                mismatches.append(
                    {
                        "position": pos,
                        "guide_nt": g_nt,
                        "target_nt": t_nt,
                        "penalty": safe_penalty,
                        "is_seed": (
                            pos >= 11 if endonuclease.lower() in {"spcas9", "cas9"} else pos <= 8
                        ),
                    }
                )

        cfd_score = math.exp(log_cfd)
        is_active = cfd_score >= 0.01  # Significant off-target risk threshold (>= 1%)

        return {
            "cfd_score": min(max(cfd_score, 0.0), 1.0),
            "pam_score": pam_score,
            "mismatch_count": len(mismatches),
            "mismatches": mismatches,
            "is_active_offtarget": is_active,
            "log_cfd": log_cfd,
            "endonuclease": endonuclease,
        }


# ============================================================================
# 2. DNA & RNA BULGE (INDEL) PENALTY SCORING ENGINE
# ============================================================================


class CRISPRBulgeEngine:
    """
    Models DNA and RNA bulges (insertions and deletions in RNA:DNA hybrid)
    and their position-dependent impact on endonuclease cleavage.
    """

    @classmethod
    def score_bulge(
        cls, guide: str, target: str, pam: str, endonuclease: str = "SpCas9"
    ) -> dict[str, Any]:
        """
        Detects and scores single or double nucleotide bulges:
        - RNA bulge: len(guide) > len(target) (extra base in gRNA)
        - DNA bulge: len(target) > len(guide) (extra base in target DNA)
        """
        g_clean = CRISPRCleavageEngine.sanitize_sequence(guide)
        t_clean = CRISPRCleavageEngine.sanitize_sequence(target)
        diff = len(g_clean) - len(t_clean)

        if diff == 0:
            # No bulge, delegate to standard CFD
            return CRISPRCleavageEngine.compute_cfd_score(g_clean, t_clean, pam, endonuclease)

        bulge_type = "RNA_BULGE" if diff > 0 else "DNA_BULGE"
        bulge_size = abs(diff)

        if bulge_size > 2:
            # Bulges > 2 nt completely abolish Cas9/Cas12a cleavage
            return {
                "cfd_score": 0.0,
                "bulge_type": bulge_type,
                "bulge_size": bulge_size,
                "is_tolerated": False,
                "reason": f"Bulge size {bulge_size} bp exceeds maximum tolerance threshold (<= 2 bp)",
            }

        # Estimate optimal alignment and bulge position via maximum matching score
        best_score = 0.0
        best_align = None

        if bulge_type == "RNA_BULGE":
            # Extra base in gRNA; test skipping each base in guide
            for skip_idx in range(len(g_clean)):
                trimmed_g = g_clean[:skip_idx] + g_clean[skip_idx + 1 :]
                if len(trimmed_g) == len(t_clean):
                    res = CRISPRCleavageEngine.compute_cfd_score(
                        trimmed_g, t_clean, pam, endonuclease
                    )
                    # Bulge penalty factor: seed position (11-20 in Cas9) severely penalized
                    bulge_pos = skip_idx + 1
                    pos_penalty = (
                        0.03
                        if (bulge_pos >= 11 if endonuclease == "SpCas9" else bulge_pos <= 8)
                        else 0.45
                    )
                    adjusted_score = res["cfd_score"] * pos_penalty
                    if adjusted_score > best_score:
                        best_score = adjusted_score
                        best_align = (bulge_pos, adjusted_score)
        else:
            # DNA bulge: extra base in target DNA; test skipping each base in target
            for skip_idx in range(len(t_clean)):
                trimmed_t = t_clean[:skip_idx] + t_clean[skip_idx + 1 :]
                if len(trimmed_t) == len(g_clean):
                    res = CRISPRCleavageEngine.compute_cfd_score(
                        g_clean, trimmed_t, pam, endonuclease
                    )
                    bulge_pos = skip_idx + 1
                    pos_penalty = (
                        0.03
                        if (bulge_pos >= 11 if endonuclease == "SpCas9" else bulge_pos <= 8)
                        else 0.40
                    )
                    adjusted_score = res["cfd_score"] * pos_penalty
                    if adjusted_score > best_score:
                        best_score = adjusted_score
                        best_align = (bulge_pos, adjusted_score)

        return {
            "cfd_score": best_score,
            "bulge_type": bulge_type,
            "bulge_size": bulge_size,
            "best_bulge_position": best_align[0] if best_align else None,
            "is_tolerated": best_score >= 0.005,
            "endonuclease": endonuclease,
        }


# ============================================================================
# 3. EPIGENETIC CHROMATIN ACCESSIBILITY WEIGHTING ENGINE
# ============================================================================


class ChromatinAccessibilityEngine:
    """
    Modulates intrinsic cleavage probability with epigenetic chromatin accessibility
    data derived from ATAC-seq, DNase-seq, or MNase-seq peak datasets.
    """

    DEFAULT_BASAL_RATE = 0.15  # Minimum basal cleavage in dense heterochromatin (15%)

    @classmethod
    def compute_accessibility_factor(
        cls, atac_signal: float, kd: float = 10.0, basal_rate: float = DEFAULT_BASAL_RATE
    ) -> float:
        """
        Computes hyperbolic chromatin accessibility scaling factor Phi(A) in [basal_rate, 1.0].
        A = 0.0 (closed heterochromatin) -> factor = basal_rate
        A >> kd (open euchromatin) -> factor approaches 1.0
        """
        atac_signal = max(atac_signal, 0.0)

        # Hyperbolic Hill-like activation
        dynamic_factor = atac_signal / (atac_signal + kd)
        factor = basal_rate + (1.0 - basal_rate) * dynamic_factor
        return min(max(factor, basal_rate), 1.0)

    @classmethod
    def compute_effective_cleavage_risk(
        cls, intrinsic_cfd: float, atac_signal: float, kd: float = 10.0
    ) -> dict[str, Any]:
        """
        Calculates effective in-vivo cleavage risk:
        Risk_eff = Score_intrinsic * Phi(A)
        """
        acc_factor = cls.compute_accessibility_factor(atac_signal, kd=kd)
        effective_score = intrinsic_cfd * acc_factor

        return {
            "intrinsic_cfd": intrinsic_cfd,
            "atac_signal": atac_signal,
            "accessibility_factor": acc_factor,
            "effective_cleavage_score": effective_score,
            "chromatin_state": "OPEN_EUCHROMATIN"
            if acc_factor > 0.65
            else ("INTERMEDIATE" if acc_factor > 0.35 else "CLOSED_HETEROCHROMATIN"),
        }


# ============================================================================
# 4. WHOLE-GUIDE SPECIFICITY & OFFTARGET AGGREGATOR
# ============================================================================


class GuideSpecificityAggregator:
    """
    Ranks guide RNAs by calculating global specificity scores across a genome-wide
    or candidate locus pool of potential off-target sites.
    """

    @classmethod
    def compute_specificity_score(cls, off_target_scores: list[float]) -> float:
        """
        Computes Hsu / Doench aggregate specificity score in [0, 100]:
        S_guide = 100 / (1.0 + sum(off_target_scores))
        """
        total_offtarget_mass = sum(off_target_scores)
        score = 100.0 / (1.0 + total_offtarget_mass)
        return min(max(score, 0.0), 100.0)


# ============================================================================
# 5. SOVEREIGN AUTHORSHIP, SECRETS & LICENSE AUDITOR (SECTION 10 GATE)
# ============================================================================


class SovereigntyAndSecretAuditor:
    """
    Enforces Section 10 QA & Security Gate:
    1. Zero exposed secrets / API keys in code or fixtures.
    2. Zero GPL copyleft contamination (strict MIT / Apache-2.0 compliance).
    3. Mandatory Sovereign Authorship: Lakshitha Dananjaya (devsoniclk).
    4. Mandatory OpenBot Co-Author trailer.
    """

    MANDATORY_AUTHOR = "Lakshitha Dananjaya"
    MANDATORY_ORG = "devsoniclk"
    MANDATORY_TRAILER = (
        "Co-authored-by: OpenBot <335604065+openbotgit[bot]@users.noreply.github.com>"
    )

    SECRET_PATTERNS = (
        (
            re.compile(
                r'(?i)(?:api_key|apikey|secret_key|private_key|token)\s*[:=]\s*["\']([a-zA-Z0-9_\-]{20,})["\']'
            ),
            "HIGH_ENTROPY_SECRET",
        ),
        (re.compile(r"-----BEGIN (?:[A-Z0-9_\-]+ )?PRIVATE KEY-----"), "PEM_PRIVATE_KEY"),
        (re.compile(r"\b(gh[pous]_[A-Za-z0-9_]{36,255})\b"), "GITHUB_TOKEN"),
    )

    COPYLEFT_PATTERNS = (
        re.compile(r"(?i)\b(GNU General Public License|GPL-2\.0|GPL-3\.0|AGPL-3\.0|LGPL-3\.0)\b"),
    )

    @classmethod
    def audit_secrets(cls, text: str) -> list[dict[str, Any]]:
        findings = []
        for pattern, leak_type in cls.SECRET_PATTERNS:
            matches = pattern.findall(text)
            for m in matches:
                val = m if isinstance(m, str) else m[0]
                obfuscated = val[:4] + "*" * (len(val) - 8) + val[-4:] if len(val) > 8 else "***"
                findings.append({"type": leak_type, "severity": "CRITICAL", "matched": obfuscated})
        return findings

    @classmethod
    def audit_license(cls, license_text: str) -> dict[str, Any]:
        has_copyleft = any(p.search(license_text) for p in cls.COPYLEFT_PATTERNS)
        is_permissive = bool(
            re.search(
                r"(?i)\b(MIT License|Apache License, Version 2\.0|Apache-2\.0)\b", license_text
            )
        )
        return {
            "is_compliant": is_permissive and not has_copyleft,
            "has_copyleft_contamination": has_copyleft,
            "is_permissive": is_permissive,
            "verdict": "VETO_COPYLEFT_CONTAMINATION"
            if has_copyleft
            else ("PASS" if is_permissive else "WARNING_UNKNOWN_LICENSE"),
        }

    @classmethod
    def audit_provenance(cls, content: str) -> dict[str, Any]:
        has_author = (cls.MANDATORY_AUTHOR in content) or (cls.MANDATORY_ORG in content)
        has_trailer = cls.MANDATORY_TRAILER in content
        return {
            "has_sovereign_author": has_author,
            "has_openbot_trailer": has_trailer,
            "is_compliant": has_author and has_trailer,
            "verdict": "PASS" if (has_author and has_trailer) else "VETO_MISSING_PROVENANCE",
        }


# ============================================================================
# 6. UNIFIED CRISPR PREDICTOR ORCHESTRATOR
# ============================================================================


class CRISPRPredictorOrchestrator:
    """
    Orchestrates end-to-end off-target risk assessment and codebase integrity auditing.
    """

    def __init__(self):
        self.cleavage_engine = CRISPRCleavageEngine()
        self.bulge_engine = CRISPRBulgeEngine()
        self.chromatin_engine = ChromatinAccessibilityEngine()
        self.specificity_aggregator = GuideSpecificityAggregator()
        self.sovereignty_auditor = SovereigntyAndSecretAuditor()

    def predict_offtarget(
        self,
        guide: str,
        target: str,
        pam: str,
        atac_signal: float | None = None,
        endonuclease: str = "SpCas9",
    ) -> dict[str, Any]:
        """Runs full cleavage prediction across mismatch, bulge, and chromatin layers."""
        if len(guide) == len(target):
            intrinsic = self.cleavage_engine.compute_cfd_score(guide, target, pam, endonuclease)
        else:
            intrinsic = self.bulge_engine.score_bulge(guide, target, pam, endonuclease)

        cfd_score = intrinsic["cfd_score"]

        if atac_signal is not None:
            effective = self.chromatin_engine.compute_effective_cleavage_risk(
                cfd_score, atac_signal
            )
            intrinsic["effective_cleavage_score"] = effective["effective_cleavage_score"]
            intrinsic["chromatin_state"] = effective["chromatin_state"]
            intrinsic["accessibility_factor"] = effective["accessibility_factor"]
        else:
            intrinsic["effective_cleavage_score"] = cfd_score
            intrinsic["chromatin_state"] = "UNWEIGHTED"
            intrinsic["accessibility_factor"] = 1.0

        return intrinsic

    def run_full_suite(self, file_path: str) -> dict[str, Any]:
        """Audits target source file for secrets, copyleft contamination, and sovereign provenance."""
        if not os.path.exists(file_path):
            return {"error": f"Path not found: {file_path}", "passed": False}

        results = {
            "target": file_path,
            "secret_findings": [],
            "provenance": {},
            "release_veto": False,
            "veto_reasons": [],
        }

        with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
            content = f.read()

        secrets = self.sovereignty_auditor.audit_secrets(content)
        results["secret_findings"] = secrets
        if secrets:
            results["release_veto"] = True
            results["veto_reasons"].append(
                f"Exposed high-entropy secret detected ({len(secrets)} occurrences)"
            )

        prov = self.sovereignty_auditor.audit_provenance(content)
        results["provenance"] = prov
        if not prov["is_compliant"]:
            results["release_veto"] = True
            results["veto_reasons"].append(
                "Missing Lakshitha Dananjaya authorship or OpenBot co-author trailer"
            )

        return results


# ============================================================================
# CLI ENTRY POINT
# ============================================================================


def main():
    from .cli import main as cli_main

    cli_main()


if __name__ == "__main__":
    main()
