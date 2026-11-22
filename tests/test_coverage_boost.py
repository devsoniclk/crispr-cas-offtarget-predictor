"""
Additional tests targeting high coverage (>90%) across engine and server.
Author: Lakshitha Dananjaya (devsoniclk)
Co-authored-by: OpenBot <335604065+openbotgit[bot]@users.noreply.github.com>
"""

import io
import json
from unittest.mock import patch

from crispr_cas_offtarget_predictor.engine import (
    ChromatinAccessibilityEngine,
    CRISPRBulgeEngine,
    CRISPRCleavageEngine,
    CRISPRPredictorOrchestrator,
    SovereigntyAndSecretAuditor,
)
from crispr_cas_offtarget_predictor.server import run_stdio_server


def test_cas9_pam_variations():
    assert CRISPRCleavageEngine.get_cas9_pam_score("TGG") == 1.000
    assert CRISPRCleavageEngine.get_cas9_pam_score("TAG") == 0.259
    assert CRISPRCleavageEngine.get_cas9_pam_score("TGA") == 0.069
    assert CRISPRCleavageEngine.get_cas9_pam_score("CCC") == 0.001
    assert CRISPRCleavageEngine.get_cas9_pam_score("TG") == 0.0


def test_cas12a_pam_variations():
    assert CRISPRCleavageEngine.get_cas12a_pam_score("TTTA") == 1.000
    assert CRISPRCleavageEngine.get_cas12a_pam_score("TTTN") == 0.200
    assert CRISPRCleavageEngine.get_cas12a_pam_score("TCTA") == 0.050
    assert CRISPRCleavageEngine.get_cas12a_pam_score("CCCC") == 0.010
    assert CRISPRCleavageEngine.get_cas12a_pam_score("TT") == 0.0


def test_dna_bulge_cas12a():
    guide = "AGATCGTCTTCCTCATCTTC"
    target = "AGATCGTCTTTCCCTCATCTTC"  # 2bp longer -> DNA bulge
    res = CRISPRBulgeEngine.score_bulge(guide, target, "TTTA", "Cas12a")
    assert res["bulge_type"] == "DNA_BULGE"
    assert res["bulge_size"] == 2


def test_bulge_size_greater_than_two():
    guide = "GAGTCCGAGCAGAAGAAGAA"
    target = "GAGTCCGAGCAGAAGAAGAATTT"  # 3bp longer
    res = CRISPRBulgeEngine.score_bulge(guide, target, "TGG", "SpCas9")
    assert res["cfd_score"] == 0.0
    assert res["is_tolerated"] is False


def test_zero_pam_score_cfd():
    res = CRISPRCleavageEngine.compute_cfd_score(
        "GAGTCCGAGCAGAAGAAGAA",
        "GAGTCCGAGCAGAAGAAGAA",
        "TG",  # Invalid short PAM -> score 0.0
        "SpCas9",
    )
    assert res["cfd_score"] == 0.0
    assert res["is_active_offtarget"] is False


def test_chromatin_edge_cases():
    fac_neg = ChromatinAccessibilityEngine.compute_accessibility_factor(-5.0)
    assert fac_neg == 0.15
    fac_high = ChromatinAccessibilityEngine.compute_accessibility_factor(1000.0)
    assert fac_high > 0.98


def test_sovereignty_secret_detection():
    pem_key = "-----BEGIN RSA PRIVATE KEY-----\nMIIEowIBAAKCAQEA0\n-----END RSA PRIVATE KEY-----"
    findings = SovereigntyAndSecretAuditor.audit_secrets(pem_key)
    assert len(findings) >= 1
    assert findings[0]["type"] == "PEM_PRIVATE_KEY"

    gh_tok = "ghp_123456789012345678901234567890123456"
    f2 = SovereigntyAndSecretAuditor.audit_secrets(gh_tok)
    assert len(f2) >= 1
    assert f2[0]["type"] == "GITHUB_TOKEN"


def test_license_auditor_cases():
    lic_unknown = "Some proprietary copyright all rights reserved"
    res = SovereigntyAndSecretAuditor.audit_license(lic_unknown)
    assert res["is_compliant"] is False
    assert res["verdict"] == "WARNING_UNKNOWN_LICENSE"


def test_run_full_suite_veto(tmp_path):
    orch = CRISPRPredictorOrchestrator()
    bad_file = tmp_path / "bad.py"
    bad_file.write_text("API_KEY = 'secret12345678901234567890'\n# missing author")
    res = orch.run_full_suite(str(bad_file))
    assert res["release_veto"] is True
    assert len(res["veto_reasons"]) >= 1


def test_run_stdio_server_loop():
    input_data = (
        json.dumps(
            {
                "jsonrpc": "2.0",
                "id": 10,
                "method": "predict_offtarget",
                "params": {
                    "guide": "GAGTCCGAGCAGAAGAAGAA",
                    "target": "GAGTCCGAGCAGAAGAAGAA",
                    "pam": "TGG",
                },
            }
        )
        + "\n\n"
    )
    with (
        patch("sys.stdin", io.StringIO(input_data)),
        patch("sys.stdout", new_callable=io.StringIO) as fake_out,
    ):
        run_stdio_server()
        out_str = fake_out.getvalue()
        assert "jsonrpc" in out_str
        data = json.loads(out_str.strip())
        assert data["id"] == 10
        assert data["result"]["cfd_score"] == 1.0
