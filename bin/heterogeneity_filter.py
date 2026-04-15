#!/usr/bin/env python3
"""
heterogeneity_filter.py

Author: Nadeem Khan, INRS-CAFSB
"""

import argparse
import pandas as pd
import numpy as np


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--input',       required=True)
    p.add_argument('--trait',       required=True)
    p.add_argument('--type',        default='qt')
    p.add_argument('--p_sig',       type=float, default=5e-8)
    p.add_argument('--p_sug',       type=float, default=1e-5)
    p.add_argument('--i2_thresh',   type=float, default=75)
    p.add_argument('--clump_p1',    type=float, default=5e-8)
    p.add_argument('--clump_p2',    type=float, default=1e-2)
    p.add_argument('--clump_r2',    type=float, default=0.05)
    p.add_argument('--clump_kb',    type=int,   default=500)
    p.add_argument('--out_loci',    required=True)
    p.add_argument('--out_summary', required=True)
    p.add_argument('--out_het',     required=True)
    return p.parse_args()


def positional_clump(df: pd.DataFrame, p_col: str, p1: float,
                     clump_kb: int) -> pd.DataFrame:
    """
    Greedy positional clumping: sort by p-value, keep lead SNP,
    remove all SNPs within clump_kb on same chromosome.
    Returns DataFrame of independent lead SNPs.
    """
    sig = df[df[p_col] <= p1].copy()
    if sig.empty:
        return sig
    sig = sig.sort_values(p_col)
    leads = []
    used  = set()
    for _, row in sig.iterrows():
        snp = row['SNP']
        if snp in used:
            continue
        leads.append(row)
        # Mark window
        window = sig[
            (sig['CHR'] == row['CHR']) &
            (sig['POS'] >= row['POS'] - clump_kb * 1000) &
            (sig['POS'] <= row['POS'] + clump_kb * 1000)
        ]['SNP']
        used.update(window)
    return pd.DataFrame(leads)


def main():
    args = parse_args()

    print(f"[filter] Loading: {args.input}")
    df = pd.read_csv(args.input, sep='\t', low_memory=False,
                     compression='gzip' if args.input.endswith('.gz') else 'infer')

    # Choose primary p-value column (GC-corrected RE preferred)
    p_col = 'p_re_gc' if 'p_re_gc' in df.columns else \
            'p_re'    if 'p_re'    in df.columns else 'p_fe'

    # ── High-heterogeneity loci ────────────────────────────────
    hi_het = df[df['i2'] > args.i2_thresh].copy() if 'i2' in df.columns \
             else pd.DataFrame()
    hi_het['trait'] = args.trait
    hi_het.to_csv(args.out_het, sep='\t', index=False)
    print(f"[filter] High-het loci (I²>{args.i2_thresh}%): {len(hi_het)}")

    # ── Positional clumping ────────────────────────────────────
    top_loci = positional_clump(df, p_col, args.clump_p1, args.clump_kb)
    top_loci['trait'] = args.trait
    top_loci['type']  = args.type

    # Direction-of-effect annotation
    beta_cols = [c for c in df.columns if c.startswith('BETA_')]
    if beta_cols and len(top_loci) > 0:
        signs = top_loci[beta_cols].apply(
            lambda row: ''.join(['+' if v > 0 else '-' if v < 0 else '?'
                                 for v in row.dropna()]), axis=1)
        top_loci['direction'] = signs

    # Heterogeneity flag
    if 'i2' in top_loci.columns:
        top_loci['hi_het'] = top_loci['i2'] > args.i2_thresh

    top_loci.to_csv(args.out_loci, sep='\t', index=False)
    print(f"[filter] Independent lead SNPs: {len(top_loci)}")

    # ── Summary ───────────────────────────────────────────────
    summ = pd.DataFrame([{
        'trait':         args.trait,
        'type':          args.type,
        'n_total':       len(df),
        'n_sig':         (df[p_col] < args.p_sig).sum(),
        'n_sug':         ((df[p_col] < args.p_sug) & (df[p_col] >= args.p_sig)).sum(),
        'n_independent': len(top_loci),
        'n_hi_het':      len(hi_het),
        'p_col_used':    p_col,
    }])
    summ.to_csv(args.out_summary, sep='\t', index=False)
    print(f"[filter] Summary written: {args.out_summary}")


if __name__ == '__main__':
    main()
