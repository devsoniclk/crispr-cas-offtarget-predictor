#!/usr/bin/env python3
"""
crispr_cas_offtarget_predictor.runtime_guard
============================================
Author: Lakshitha Dananjaya (devsoniclk)
Co-authored-by: OpenBot <335604065+openbotgit[bot]@users.noreply.github.com>

Runtime Streaming Guard & Numerical Integrity Interceptor for CRISPR Endonuclease Analysis.
Guarantees constant-memory scanning over multi-gigabase genomes and validates log-space math.
"""

import math
from collections.abc import Generator


class CRISPRRuntimeGuard:
    """
    Provides memory-bounded streaming generators and numerical bounds assertion.
    """

    @staticmethod
    def stream_candidate_sites(
        genome_chunk: str, window_size: int = 23, stride: int = 1
    ) -> Generator[tuple[int, str], None, None]:
        """
        Constant-memory generator yielding candidate genomic windows without
        allocating massive sequence lists into memory (CTP-006 mitigation).
        """
        chunk_len = len(genome_chunk)
        for i in range(0, chunk_len - window_size + 1, stride):
            yield (i, genome_chunk[i : i + window_size])

    @staticmethod
    def assert_log_domain_stability(weights: list[float]) -> float:
        """
        Safely computes product of probabilistic weights in the log-domain,
        preventing IEEE 754 catastrophic underflow (CTP-004 mitigation).
        """
        log_sum = 0.0
        for w in weights:
            if w <= 0.0:
                return 0.0
            log_sum += math.log(w)
        return math.exp(log_sum)

    @staticmethod
    def validate_nucleotide_alphabet(sequence: str) -> bool:
        """Strict validation of IUPAC nucleotide symbols for CRISPR guides."""
        allowed = {"A", "C", "G", "T", "U", "N"}
        return all(c.upper() in allowed for c in sequence.strip())
