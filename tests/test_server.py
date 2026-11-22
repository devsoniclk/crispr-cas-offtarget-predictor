"""
Tests for FastMCP JSON-RPC stdio Server.
Author: Lakshitha Dananjaya (devsoniclk)
Co-authored-by: OpenBot <335604065+openbotgit[bot]@users.noreply.github.com>
"""

import pytest

from crispr_cas_offtarget_predictor.fixtures import (
    EMX1_GUIDE,
    EMX1_OFFTARGET_1,
    EMX1_OFFTARGET_BULGE,
    EMX1_PAM,
)
from crispr_cas_offtarget_predictor.server import CRISPROffTargetMCPServer


@pytest.fixture
def server():
    return CRISPROffTargetMCPServer()


def test_server_predict_offtarget(server):
    res = server.predict_offtarget(
        guide=EMX1_GUIDE,
        target=EMX1_OFFTARGET_1,
        pam=EMX1_PAM,
        atac_signal=25.0,
        endonuclease="SpCas9",
    )
    assert "cfd_score" in res
    assert "effective_cleavage_score" in res
    assert res["cfd_score"] > 0.0
    assert res["accessibility_factor"] > 0.5


def test_server_score_bulge(server):
    res = server.score_bulge(
        guide=EMX1_GUIDE,
        target=EMX1_OFFTARGET_BULGE,
        pam=EMX1_PAM,
        endonuclease="SpCas9",
    )
    assert "bulge_type" in res
    assert res["bulge_type"] == "RNA_BULGE"
    assert res["bulge_size"] == 1


def test_server_modulate_chromatin(server):
    res = server.modulate_chromatin(intrinsic_cfd=0.8, atac_signal=0.0)
    assert res["accessibility_factor"] == 0.15
    assert res["effective_cleavage_score"] == pytest.approx(0.8 * 0.15)


def test_server_compute_guide_specificity(server):
    res = server.compute_guide_specificity([0.05, 0.02, 0.01])
    assert "specificity_score" in res
    assert res["specificity_score"] > 90.0
    assert res["tier"] == "HIGH_SPECIFICITY"


def test_server_audit_codebase(server):
    import os

    import crispr_cas_offtarget_predictor

    pkg_dir = os.path.dirname(crispr_cas_offtarget_predictor.__file__)
    engine_file = os.path.join(pkg_dir, "engine.py")
    res = server.audit_codebase(engine_file)
    assert res["release_veto"] is False
    assert res["provenance"]["is_compliant"] is True


def test_server_jsonrpc_dispatcher(server):
    req = {
        "jsonrpc": "2.0",
        "id": 1,
        "method": "predict_offtarget",
        "params": {
            "guide": EMX1_GUIDE,
            "target": EMX1_OFFTARGET_1,
            "pam": EMX1_PAM,
            "atac_signal": 15.0,
        },
    }
    resp = server.handle_request(req)
    assert resp["id"] == 1
    assert "result" in resp
    assert resp["result"]["cfd_score"] > 0.0

    # Test tools/list
    list_req = {"jsonrpc": "2.0", "id": 2, "method": "tools/list"}
    list_resp = server.handle_request(list_req)
    assert len(list_resp["result"]["tools"]) == 5

    # Test unknown method
    err_req = {"jsonrpc": "2.0", "id": 3, "method": "nonexistent_method"}
    err_resp = server.handle_request(err_req)
    assert "error" in err_resp
    assert err_resp["error"]["code"] == -32601


def test_server_all_jsonrpc_methods(server):
    # Test score_bulge
    b_req = {
        "jsonrpc": "2.0",
        "id": 11,
        "method": "score_bulge",
        "params": {"guide": "GAGTCCGAGCAGAAGAAGAA", "target": "GAGTCCGAGCAGAGAAGAA", "pam": "TGG"},
    }
    b_resp = server.handle_request(b_req)
    assert b_resp["result"]["bulge_type"] == "RNA_BULGE"

    # Test modulate_chromatin
    m_req = {
        "jsonrpc": "2.0",
        "id": 12,
        "method": "modulate_chromatin",
        "params": {"intrinsic_cfd": 0.5, "atac_signal": 20.0, "kd": 10.0},
    }
    m_resp = server.handle_request(m_req)
    assert "effective_cleavage_score" in m_resp["result"]

    # Test compute_guide_specificity
    s_req = {
        "jsonrpc": "2.0",
        "id": 13,
        "method": "compute_guide_specificity",
        "params": {"off_target_scores": [0.1, 0.05]},
    }
    s_resp = server.handle_request(s_req)
    assert "specificity_score" in s_resp["result"]

    # Test audit_codebase
    import os

    import crispr_cas_offtarget_predictor

    pkg_dir = os.path.dirname(crispr_cas_offtarget_predictor.__file__)
    a_req = {
        "jsonrpc": "2.0",
        "id": 14,
        "method": "audit_codebase",
        "params": {"file_path": os.path.join(pkg_dir, "engine.py")},
    }
    a_resp = server.handle_request(a_req)
    assert a_resp["result"]["release_veto"] is False
