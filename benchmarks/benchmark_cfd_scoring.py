#!/usr/bin/env python3
"""
Benchmark: Computational Throughput of Log-Space CFD Scoring Engine.
"""

import time

from crispr_cas_offtarget_predictor import EMX1_GUIDE, EMX1_PAM, CRISPRCleavageEngine

ITERATIONS = 50_000
target = "GAGTCCGAGGAGAAGAAGAA"

start = time.perf_counter()
for _ in range(ITERATIONS):
    CRISPRCleavageEngine.compute_cfd_score(EMX1_GUIDE, target, EMX1_PAM, "SpCas9")
elapsed = time.perf_counter() - start

ops_sec = ITERATIONS / elapsed
print(f"Scored {ITERATIONS:,} target sequences in {elapsed:.3f}s ({ops_sec:,.1f} scores/sec)")
