#!/usr/bin/env python3
"""
crispr_cas_offtarget_predictor.fixtures
=======================================
Author: Lakshitha Dananjaya (devsoniclk)
Co-authored-by: OpenBot <335604065+openbotgit[bot]@users.noreply.github.com>

Reference benchmarks, validated genomic targets, and threat fixtures
for automated CRISPR-Cas guide RNA off-target cleavage auditing.
"""

# Benchmark human genomic loci verified via GUIDE-seq / CIRCLE-seq
EMX1_GUIDE = "GAGTCCGAGCAGAAGAAGAA"  # SpCas9 site (chr2:73160475)
EMX1_PAM = "TGG"

# Validated in-vivo off-target sites for EMX1
EMX1_OFFTARGET_1 = "GAGTCCGAGGAGAAGAAGAA"  # Distal mismatch at pos 9 (C->G) - HIGH cleavage
EMX1_OFFTARGET_2 = "GAGTCCGAGCAGAAGAAGTA"  # Seed mismatch at pos 19 (A->T) - LOW cleavage
EMX1_OFFTARGET_BULGE = "GAGTCCGAGCAGAGAAGAA"  # 1bp RNA bulge (deletion of A at pos 15)

# VEGFA Site 2 (Promiscuous SpCas9 target with known non-canonical PAM cleavage)
VEGFA_GUIDE = "GACCCCCTCCACCCCGCCTC"
VEGFA_PAM_CANONICAL = "CGG"
VEGFA_PAM_NAG = "CAG"  # Non-canonical NAG PAM off-target

# Cas12a (Cpf1) DNMT1 target (5' TTTV PAM)
CAS12A_DNMT1_GUIDE = "AGATCGTCTTCCTCATCTTC"
CAS12A_DNMT1_PAM = "TTTA"

# ============================================================================
# FLAWED / VULNERABLE BASELINE MODELS FOR AUDIT VALIDATION
# ============================================================================


def vulnerable_unweighted_hamming_score(guide: str, target: str) -> float:
    """
    VULNERABLE (CTP-001): Treats all mismatches equally regardless of seed proximity.
    Fails to catch high-affinity seed-conserved off-target cuts.
    """
    mismatches = sum(1 for g, t in zip(guide, target) if g != t)
    return max(1.0 - (0.2 * mismatches), 0.0)


def vulnerable_canonical_pam_only_score(pam: str) -> float:
    """
    VULNERABLE (CTP-002): Strictly rejects any non-NGG PAM, completely missing
    NAG/NGA off-target cleavage which can account for up to 25% of in-vivo DSBs.
    """
    return 1.0 if pam.endswith("GG") else 0.0


def vulnerable_multiplicative_underflow(mismatch_weights: list) -> float:
    """
    VULNERABLE (CTP-004): Direct IEEE 754 floating-point multiplication of many small
    weights leads to subnormal float underflow and premature zeroing.
    """
    prod = 1.0
    for w in mismatch_weights:
        prod *= w
    return prod
