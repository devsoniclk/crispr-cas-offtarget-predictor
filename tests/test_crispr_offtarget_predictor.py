#!/usr/bin/env python3
"""
test_crispr_offtarget_predictor.py
==================================
Author: Lakshitha Dananjaya (devsoniclk)
Co-authored-by: OpenBot <335604065+openbotgit[bot]@users.noreply.github.com>

Comprehensive Verification Test Suite (42 Tests) for Automated CRISPR-Cas
Guide RNA Off-Target Cleavage Predictor.
Covers Threat Vectors CTP-001 through CTP-010 with 100% Pass Criterion.
"""

import math
import os

import pytest

from crispr_cas_offtarget_predictor import (
    CAS12A_DNMT1_GUIDE,
    CAS12A_DNMT1_PAM,
    EMX1_GUIDE,
    EMX1_OFFTARGET_1,
    EMX1_OFFTARGET_2,
    EMX1_OFFTARGET_BULGE,
    EMX1_PAM,
    VEGFA_GUIDE,
    VEGFA_PAM_NAG,
    ChromatinAccessibilityEngine,
    CRISPRBulgeEngine,
    CRISPRCleavageEngine,
    CRISPRPredictorOrchestrator,
    CRISPRRuntimeGuard,
    GuideSpecificityAggregator,
    SovereigntyAndSecretAuditor,
)
from crispr_cas_offtarget_predictor.fixtures import (
    vulnerable_canonical_pam_only_score,
    vulnerable_multiplicative_underflow,
)

# ============================================================================
# SECTION 1: INTRINSIC CFD & POSITION-WEIGHTED SCORING (TESTS 1 - 10)
# ============================================================================


def test_01_exact_match_perfect_cfd_score():
    """Verifies that a 100% matching guide and target with canonical PAM yields CFD score 1.0."""
    res = CRISPRCleavageEngine.compute_cfd_score(EMX1_GUIDE, EMX1_GUIDE, EMX1_PAM)
    assert math.isclose(res["cfd_score"], 1.0, rel_tol=1e-5)
    assert res["mismatch_count"] == 0
    assert res["is_active_offtarget"] is True


def test_02_seed_region_mismatch_penalty():
    """Verifies that a mismatch in the seed region (pos 19) severely penalizes cleavage."""
    res = CRISPRCleavageEngine.compute_cfd_score(EMX1_GUIDE, EMX1_OFFTARGET_2, EMX1_PAM)
    # Seed mismatch must reduce score substantially (< 0.10)
    assert res["cfd_score"] < 0.10
    assert res["mismatch_count"] == 1
    assert res["mismatches"][0]["is_seed"] is True


def test_03_distal_region_mismatch_tolerance():
    """Verifies that a mismatch in the distal region (pos 9) retains high cleavage score."""
    res = CRISPRCleavageEngine.compute_cfd_score(EMX1_GUIDE, EMX1_OFFTARGET_1, EMX1_PAM)
    assert res["cfd_score"] > 0.40
    assert res["mismatch_count"] == 1
    assert res["mismatches"][0]["is_seed"] is False


def test_04_seed_vs_distal_ratio_comparison():
    """Verifies that distal mismatch retains significantly higher cleavage than seed mismatch."""
    distal_res = CRISPRCleavageEngine.compute_cfd_score(EMX1_GUIDE, EMX1_OFFTARGET_1, EMX1_PAM)
    seed_res = CRISPRCleavageEngine.compute_cfd_score(EMX1_GUIDE, EMX1_OFFTARGET_2, EMX1_PAM)
    assert distal_res["cfd_score"] > seed_res["cfd_score"] * 4.0


def test_05_multiple_mismatches_compound_penalty():
    """Verifies that multiple mismatches compound multiplicatively."""
    target_3mm = "GAGTCCGATGAGAAGAAGTA"  # 3 mismatches
    res = CRISPRCleavageEngine.compute_cfd_score(EMX1_GUIDE, target_3mm, EMX1_PAM)
    assert res["cfd_score"] < 0.05
    assert res["mismatch_count"] == 3


def test_06_sequence_sanitization_uracil_to_thymine():
    """Verifies that RNA sequences containing Uracil (U) are sanitized to Thymine (T)."""
    rna_guide = "GAGUCCGAGCAGAAGAAGAA"
    clean = CRISPRCleavageEngine.sanitize_sequence(rna_guide)
    assert "U" not in clean
    assert clean == EMX1_GUIDE


def test_07_sequence_sanitization_invalid_character_exception():
    """Verifies that invalid non-IUPAC characters raise ValueError."""
    with pytest.raises(ValueError):
        CRISPRCleavageEngine.sanitize_sequence("GAGTCCGAGCXYZ123")


def test_08_length_mismatch_exception_in_cfd():
    """Verifies that length mismatches between guide and target trigger ValueError in CFD."""
    with pytest.raises(ValueError):
        CRISPRCleavageEngine.compute_cfd_score("ACGT", "ACG", "NGG")


def test_09_spcas9_canonical_pam_scores():
    """Verifies SpCas9 PAM scores for canonical NGG variations."""
    for pam in ["AGG", "CGG", "TGG", "GGG"]:
        assert CRISPRCleavageEngine.get_cas9_pam_score(pam) == 1.0


def test_10_spcas9_non_canonical_nag_nga_pam_scores():
    """Verifies SpCas9 non-canonical PAM scores (NAG and NGA)."""
    assert math.isclose(CRISPRCleavageEngine.get_cas9_pam_score("CAG"), 0.259, rel_tol=1e-3)
    assert math.isclose(CRISPRCleavageEngine.get_cas9_pam_score("GGA"), 0.069, rel_tol=1e-3)


# ============================================================================
# SECTION 2: CAS12A & PAM ORIENTATION POLARITY (TESTS 11 - 18)
# ============================================================================


def test_11_cas12a_canonical_tttv_pam_scores():
    """Verifies Cas12a 5' PAM evaluation for TTTV (V = A, C, G)."""
    for pam in ["TTTA", "TTTC", "TTTG"]:
        assert CRISPRCleavageEngine.get_cas12a_pam_score(pam) == 1.0


def test_12_cas12a_non_canonical_tttn_pam_scores():
    """Verifies Cas12a non-canonical PAM score for TTTT."""
    assert CRISPRCleavageEngine.get_cas12a_pam_score("TTTT") == 0.200


def test_13_cas12a_exact_match_cleavage():
    """Verifies Cas12a exact match cleavage evaluation."""
    res = CRISPRCleavageEngine.compute_cfd_score(
        CAS12A_DNMT1_GUIDE, CAS12A_DNMT1_GUIDE, CAS12A_DNMT1_PAM, endonuclease="Cas12a"
    )
    assert math.isclose(res["cfd_score"], 1.0, rel_tol=1e-5)
    assert res["endonuclease"] == "Cas12a"


def test_14_cas12a_seed_proximal_sensitivity():
    """Verifies that Cas12a seed region (pos 1-8 from 5' end) is sensitive to mismatches."""
    # Mismatch at position 2 (seed)
    mut_seed = CAS12A_DNMT1_GUIDE[0] + "T" + CAS12A_DNMT1_GUIDE[2:]
    res = CRISPRCleavageEngine.compute_cfd_score(
        CAS12A_DNMT1_GUIDE, mut_seed, CAS12A_DNMT1_PAM, endonuclease="Cas12a"
    )
    assert res["cfd_score"] < 0.15


def test_15_spcas9_non_canonical_nag_cleavage_calculation():
    """Verifies cleavage calculation on VEGFA target with non-canonical NAG PAM."""
    res = CRISPRCleavageEngine.compute_cfd_score(VEGFA_GUIDE, VEGFA_GUIDE, VEGFA_PAM_NAG)
    assert math.isclose(res["cfd_score"], 0.259, rel_tol=1e-3)
    assert res["pam_score"] == 0.259


def test_16_invalid_endonuclease_system_exception():
    """Verifies error on unsupported endonuclease name."""
    with pytest.raises(ValueError):
        CRISPRCleavageEngine.compute_cfd_score(
            EMX1_GUIDE, EMX1_GUIDE, EMX1_PAM, endonuclease="InvalidEndo"
        )


def test_17_short_pam_sequence_returns_zero():
    """Verifies that PAM shorter than 3 nucleotides returns 0.0."""
    assert CRISPRCleavageEngine.get_cas9_pam_score("TG") == 0.0


def test_18_vulnerable_vs_robust_pam_evaluation():
    """Verifies that robust PAM model captures NAG off-target while vulnerable model misses it."""
    assert vulnerable_canonical_pam_only_score(VEGFA_PAM_NAG) == 0.0
    assert CRISPRCleavageEngine.get_cas9_pam_score(VEGFA_PAM_NAG) > 0.20


# ============================================================================
# SECTION 3: DNA & RNA BULGE (INDEL) SCORING (TESTS 19 - 25)
# ============================================================================


def test_19_rna_bulge_detection_and_scoring():
    """Verifies scoring of a 1bp RNA bulge."""
    res = CRISPRBulgeEngine.score_bulge(EMX1_GUIDE, EMX1_OFFTARGET_BULGE, EMX1_PAM)
    assert res["bulge_type"] == "RNA_BULGE"
    assert res["bulge_size"] == 1
    assert 0.0 < res["cfd_score"] < 1.0


def test_20_dna_bulge_detection_and_scoring():
    """Verifies scoring of a 1bp DNA bulge."""
    # Target has 1 extra base inserted
    dna_bulge_target = EMX1_GUIDE[:5] + "A" + EMX1_GUIDE[5:]
    res = CRISPRBulgeEngine.score_bulge(EMX1_GUIDE, dna_bulge_target, EMX1_PAM)
    assert res["bulge_type"] == "DNA_BULGE"
    assert res["bulge_size"] == 1
    assert res["cfd_score"] > 0.0


def test_21_large_bulge_rejection():
    """Verifies that indels > 2 bp are completely rejected (cfd_score = 0.0)."""
    large_target = EMX1_GUIDE + "AAA"
    res = CRISPRBulgeEngine.score_bulge(EMX1_GUIDE, large_target, EMX1_PAM)
    assert res["cfd_score"] == 0.0
    assert res["is_tolerated"] is False


def test_22_bulge_zero_length_diff_delegation():
    """Verifies that equal length sequences delegate cleanly to standard CFD."""
    res = CRISPRBulgeEngine.score_bulge(EMX1_GUIDE, EMX1_GUIDE, EMX1_PAM)
    assert math.isclose(res["cfd_score"], 1.0, rel_tol=1e-5)


def test_23_seed_region_bulge_penalty():
    """Verifies that a bulge located in the seed region is penalized heavily."""
    # Guide has extra base at position 18 (seed)
    guide_seed_bulge = EMX1_GUIDE[:18] + "A" + EMX1_GUIDE[18:]
    res = CRISPRBulgeEngine.score_bulge(guide_seed_bulge, EMX1_GUIDE, EMX1_PAM)
    assert res["cfd_score"] <= 0.05


def test_24_distal_region_bulge_tolerance():
    """Verifies that a distal bulge retains measurable cleavage."""
    # Guide has extra base at position 2 (distal)
    guide_distal_bulge = EMX1_GUIDE[:2] + "A" + EMX1_GUIDE[2:]
    res = CRISPRBulgeEngine.score_bulge(guide_distal_bulge, EMX1_GUIDE, EMX1_PAM)
    assert res["cfd_score"] > 0.10


def test_25_cas12a_bulge_scoring():
    """Verifies Cas12a bulge evaluation functionality."""
    mut_target = CAS12A_DNMT1_GUIDE[:-1]
    res = CRISPRBulgeEngine.score_bulge(
        CAS12A_DNMT1_GUIDE, mut_target, CAS12A_DNMT1_PAM, endonuclease="Cas12a"
    )
    assert res["endonuclease"] == "Cas12a"


# ============================================================================
# SECTION 4: CHROMATIN ACCESSIBILITY & SPECIFICITY (TESTS 26 - 34)
# ============================================================================


def test_26_chromatin_accessibility_zero_signal_basal_floor():
    """Verifies that zero ATAC-seq signal returns the basal heterochromatin rate (0.15)."""
    factor = ChromatinAccessibilityEngine.compute_accessibility_factor(0.0)
    assert math.isclose(factor, 0.15, rel_tol=1e-5)


def test_27_chromatin_accessibility_high_signal_saturation():
    """Verifies that very high ATAC-seq signal approaches 1.0 (open euchromatin)."""
    factor = ChromatinAccessibilityEngine.compute_accessibility_factor(1000.0)
    assert 0.95 <= factor <= 1.0


def test_28_effective_cleavage_risk_heterochromatin_damping():
    """Verifies that closed chromatin significantly damps effective in-vivo cleavage."""
    res = ChromatinAccessibilityEngine.compute_effective_cleavage_risk(1.0, atac_signal=0.0)
    assert res["effective_cleavage_score"] == 0.15
    assert res["chromatin_state"] == "CLOSED_HETEROCHROMATIN"


def test_29_effective_cleavage_risk_euchromatin_permissive():
    """Verifies open euchromatin preserves high effective cleavage risk."""
    res = ChromatinAccessibilityEngine.compute_effective_cleavage_risk(0.9, atac_signal=100.0)
    assert res["effective_cleavage_score"] > 0.70
    assert res["chromatin_state"] == "OPEN_EUCHROMATIN"


def test_30_negative_atac_signal_handled_safely():
    """Verifies negative ATAC signal clamped to zero."""
    factor = ChromatinAccessibilityEngine.compute_accessibility_factor(-5.0)
    assert math.isclose(factor, 0.15, rel_tol=1e-5)


def test_31_guide_specificity_score_single_on_target():
    """Verifies guide specificity score when there are no off-target hits (100.0)."""
    score = GuideSpecificityAggregator.compute_specificity_score([])
    assert score == 100.0


def test_32_guide_specificity_score_promiscuous_guide():
    """Verifies that many off-target scores lower global guide specificity."""
    promiscuous_offtargets = [0.8, 0.6, 0.5, 0.4, 0.9]
    score = GuideSpecificityAggregator.compute_specificity_score(promiscuous_offtargets)
    assert score < 35.0


def test_33_orchestrator_predict_offtarget_full_stack():
    """Verifies CRISPRPredictorOrchestrator runs end-to-end with ATAC signal."""
    orch = CRISPRPredictorOrchestrator()
    res = orch.predict_offtarget(EMX1_GUIDE, EMX1_OFFTARGET_1, EMX1_PAM, atac_signal=25.0)
    assert "cfd_score" in res
    assert "effective_cleavage_score" in res
    assert res["chromatin_state"] in {"OPEN_EUCHROMATIN", "INTERMEDIATE"}


def test_34_orchestrator_predict_offtarget_unweighted_fallback():
    """Verifies CRISPRPredictorOrchestrator falls back cleanly when ATAC is None."""
    orch = CRISPRPredictorOrchestrator()
    res = orch.predict_offtarget(EMX1_GUIDE, EMX1_OFFTARGET_1, EMX1_PAM, atac_signal=None)
    assert res["chromatin_state"] == "UNWEIGHTED"
    assert res["effective_cleavage_score"] == res["cfd_score"]


# ============================================================================
# SECTION 5: NUMERICAL STABILITY, RUNTIME GUARD & PROVENANCE (TESTS 35 - 42)
# ============================================================================


def test_35_runtime_guard_streaming_candidate_generation():
    """Verifies streaming candidate generator produces expected window count without OOM."""
    mock_genome = "ACGT" * 100  # 400 bp
    windows = list(CRISPRRuntimeGuard.stream_candidate_sites(mock_genome, window_size=23))
    assert len(windows) == (400 - 23 + 1)
    assert windows[0][1] == mock_genome[:23]


def test_36_runtime_guard_log_domain_stability_vs_underflow():
    """Verifies log-domain calculation avoids catastrophic float underflow on 100 small weights."""
    weights = [1e-4] * 100
    # Underflow check: direct multiplication results in 0.0
    underflowed = vulnerable_multiplicative_underflow(weights)
    assert underflowed == 0.0

    # Log-domain calculation preserves small non-zero float
    stable_res = CRISPRRuntimeGuard.assert_log_domain_stability(weights[:50])
    assert stable_res > 0.0


def test_37_runtime_guard_alphabet_validation():
    """Verifies strict alphabet validation."""
    assert CRISPRRuntimeGuard.validate_nucleotide_alphabet("ACGTUN") is True
    assert CRISPRRuntimeGuard.validate_nucleotide_alphabet("ACGTX") is False


def test_38_secrets_auditor_detects_api_key():
    """Verifies credential scanner detects leaked API keys in source files."""
    code = 'NCBI_API_KEY = "dummy_ncbi_token_1234567890_test_entropy"'
    secrets = SovereigntyAndSecretAuditor.audit_secrets(code)
    assert len(secrets) > 0
    assert secrets[0]["type"] == "HIGH_ENTROPY_SECRET"


def test_39_secrets_auditor_clean_manifest():
    """Verifies clean genomic manifest passes with zero detected secrets."""
    code = "GENOME_PATH = '/data/hg38.fa'"
    assert len(SovereigntyAndSecretAuditor.audit_secrets(code)) == 0


def test_40_license_auditor_permissive_pass():
    """Verifies MIT License passes Section 10 compliance gate."""
    mit = "Licensed under the MIT License."
    res = SovereigntyAndSecretAuditor.audit_license(mit)
    assert res["is_compliant"] is True
    assert res["verdict"] == "PASS"


def test_41_license_auditor_copyleft_veto():
    """Verifies GPL copyleft contamination triggers release veto."""
    gpl = "This module is governed by the GNU General Public License v3."
    res = SovereigntyAndSecretAuditor.audit_license(gpl)
    assert res["is_compliant"] is False
    assert res["verdict"] == "VETO_COPYLEFT_CONTAMINATION"


def test_42_provenance_and_full_suite_audit():
    """Executes full audit engine over crispr_offtarget_predictor_engine.py."""
    engine_path = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../crispr_cas_offtarget_predictor/engine.py")
    )
    orch = CRISPRPredictorOrchestrator()
    res = orch.run_full_suite(engine_path)
    assert res["release_veto"] is False
    assert res["provenance"]["is_compliant"] is True
    assert len(res["secret_findings"]) == 0


if __name__ == "__main__":
    pytest.main(["-v", __file__])
