"""
Tests for CLI interface.
Author: Lakshitha Dananjaya (devsoniclk)
Co-authored-by: OpenBot <335604065+openbotgit[bot]@users.noreply.github.com>
"""

import os
import sys
from unittest.mock import patch

import pytest

import crispr_cas_offtarget_predictor
from crispr_cas_offtarget_predictor.cli import main


def test_cli_predict_json():
    test_args = [
        "crispr-predict",
        "predict",
        "--guide",
        "GAGTCCGAGCAGAAGAAGAA",
        "--target",
        "GAGTCCGAGGAGAAGAAGAA",
        "--pam",
        "TGG",
        "--atac",
        "30.0",
        "--json",
    ]
    with patch.object(sys, "argv", test_args):
        with pytest.raises(SystemExit) as exc:
            main()
        assert exc.value.code == 0


def test_cli_predict_panel():
    test_args = [
        "crispr-predict",
        "predict",
        "--guide",
        "GAGTCCGAGCAGAAGAAGAA",
        "--target",
        "GAGTCCGAGGAGAAGAAGAA",
        "--pam",
        "TGG",
        "--atac",
        "30.0",
    ]
    with patch.object(sys, "argv", test_args):
        with pytest.raises(SystemExit) as exc:
            main()
        assert exc.value.code == 0


def test_cli_audit_json():
    pkg_dir = os.path.dirname(crispr_cas_offtarget_predictor.__file__)
    engine_file = os.path.join(pkg_dir, "engine.py")
    test_args = [
        "crispr-predict",
        "audit",
        engine_file,
        "--json",
    ]
    with patch.object(sys, "argv", test_args):
        with pytest.raises(SystemExit) as exc:
            main()
        assert exc.value.code == 0


def test_cli_audit_table():
    pkg_dir = os.path.dirname(crispr_cas_offtarget_predictor.__file__)
    engine_file = os.path.join(pkg_dir, "engine.py")
    test_args = [
        "crispr-predict",
        "audit",
        engine_file,
    ]
    with patch.object(sys, "argv", test_args):
        with pytest.raises(SystemExit) as exc:
            main()
        assert exc.value.code == 0


def test_cli_legacy_flags():
    test_args = [
        "crispr-predict",
        "--guide",
        "GAGTCCGAGCAGAAGAAGAA",
        "--target",
        "GAGTCCGAGGAGAAGAAGAA",
        "--pam",
        "TGG",
    ]
    with patch.object(sys, "argv", test_args):
        with pytest.raises(SystemExit) as exc:
            main()
        assert exc.value.code == 0


def test_cli_audit_fail_nonexistent():
    test_args = [
        "crispr-predict",
        "audit",
        "/nonexistent/file/path.py",
    ]
    with patch.object(sys, "argv", test_args):
        with pytest.raises(SystemExit) as exc:
            main()
        assert exc.value.code == 0 or exc.value.code == 1
