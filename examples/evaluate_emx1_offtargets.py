#!/usr/bin/env python3
"""
Example: EMX1 Locus Guide RNA Off-Target Cleavage & Chromatin Accessibility Evaluation.
"""

from crispr_cas_offtarget_predictor import (
    EMX1_GUIDE,
    EMX1_OFFTARGET_1,
    EMX1_OFFTARGET_2,
    EMX1_OFFTARGET_BULGE,
    EMX1_PAM,
    CRISPRPredictorOrchestrator,
)

orchestrator = CRISPRPredictorOrchestrator()

sites = [
    ("On-Target Exact Match", EMX1_GUIDE, 50.0),
    ("Distal Mismatch (Pos 9 C->G, Active)", EMX1_OFFTARGET_1, 45.0),
    ("Seed Mismatch (Pos 19 A->T, Tolerable)", EMX1_OFFTARGET_2, 45.0),
    ("RNA Bulge (1bp Indel)", EMX1_OFFTARGET_BULGE, 20.0),
]

print("=== CRISPR-Cas Off-Target Cleavage Evaluation for EMX1 ===")
print(f"Guide: {EMX1_GUIDE} (PAM: {EMX1_PAM})\n")

for label, target_seq, atac_signal in sites:
    res = orchestrator.predict_offtarget(
        guide=EMX1_GUIDE,
        target=target_seq,
        pam=EMX1_PAM,
        atac_signal=atac_signal,
        endonuclease="SpCas9",
    )
    print(f"[{label}]")
    print(f"  Target: {target_seq}")
    print(f"  CFD Score: {res['cfd_score']:.5f}")
    print(f"  ATAC Signal: {atac_signal:.1f} ({res['chromatin_state']})")
    print(f"  Effective Risk: {res['effective_cleavage_score']:.5f}")
    print(f"  Active Risk: {'YES' if res.get('is_active_offtarget') else 'NO'}\n")
