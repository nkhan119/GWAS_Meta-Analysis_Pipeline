#!/usr/bin/env python3
"""
validate_sumstats.py

Author: Nadeem Khan, INRS-CAFSB
"""

import argparse
import gzip
import sys
import os
import pandas as pd
import numpy as np


REQUIRED_COLS = ['SNP', 'CHR', 'POS', 'EA', 'OA', 'BETA', 'SE', 'P', 'N']


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--trait',       required=True)
    p.add_argument('--type',        default='qt')
    p.add_argument('--cohort',      action='append', default=[],
                   help='cohort:filepath pairs, can repeat')
    p.add_argument('--min_n',       type=int,   default=2)
    p.add_argument('--p_col',       default='P')
    p.add_argument('--beta_col',    default='BETA')
    p.add_argument('--se_col',      default='SE')
    p.add_argument('--snp_col',     default='SNP')
    p.add_argument('--chr_col',     default='CHR')
    p.add_argument('--pos_col',     default='POS')
    p.add_argument('--ea_col',      default='EA')
    p.add_argument('--oa_col',      default='OA')
    p.add_argument('--eaf_col',     default='EAF')
    p.add_argument('--n_col',       default='N')
    p.add_argument('--out_prefix',  required=True)
    p.add_argument('--report',      required=True)
    return p.parse_args()


def load_sumstats(filepath: str, cohort: str) -> pd.DataFrame:
    print(f"[validate] Loading {cohort}: {filepath}")
    df = pd.read_csv(filepath, sep='\t', low_memory=False,
                     compression='gzip' if filepath.endswith('.gz') else 'infer')
    # Standardise column names (upper-case)
    df.columns = [c.upper() for c in df.columns]

    # ── Alias non-standard column names to the pipeline standard ──
    # Covers common PLINK/REGENIE/BOLT/SAIGE output variants.
    ALIASES = {
        # Position
        'BP':              'POS',
        'BASE_PAIR_LOCATION': 'POS',
        # Effect allele / other allele
        'A1':             'EA',
        'ALLELE1':        'EA',
        'ALT':            'EA',
        'A2':             'OA',
        'ALLELE0':        'OA',
        'REF':            'OA',
        # P-value
        'PVALUE':         'P',
        'P_BOLT_LMM':     'P',
        'P_BOLT_LMM_INF': 'P',
        # Effect size
        'OR':             'BETA',
        'Z':              'BETA',
        # Standard error
        'SE_BOLT_LMM':    'SE',
        # Sample size
        'N_OBS':          'N',
        'NOBS':           'N',
        'OBS_CT':         'N',
        # SNP identifier
        'RSID':           'SNP',
        'ID':             'SNP',
        'SNPID':          'SNP',
        'MARKERNAME':     'SNP',
        # Effect-allele frequency
        'A1FREQ':         'EAF',
        'MAF':            'EAF',
        'AF1':            'EAF',
    }
    for src, tgt in ALIASES.items():
        if src in df.columns and tgt not in df.columns:
            df = df.rename(columns={src: tgt})

    # Derive P from LOG10P if P is absent
    if 'P' not in df.columns and 'LOG10P' in df.columns:
        df['P'] = 10.0 ** (-df['LOG10P'].astype(float))
        print(f"[validate] {cohort}: derived P from LOG10P")

    return df


def qc_sumstats(df: pd.DataFrame, cohort: str) -> tuple[pd.DataFrame, dict]:
    stats = {'cohort': cohort}
    stats['n_raw'] = len(df)

    # Check required columns
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        print(f"[validate] WARNING: {cohort} missing columns: {missing}")
        for c in missing:
            df[c] = np.nan

    # Remove duplicates
    df = df.drop_duplicates(subset=['SNP'])
    stats['n_after_dedup'] = len(df)

    # Remove missing BETA/SE/P
    df = df.dropna(subset=['BETA', 'SE', 'P'])
    stats['n_after_missing'] = len(df)

    # Filter SE ≤ 0
    df = df[df['SE'] > 0]
    stats['n_after_se'] = len(df)

    # Filter extreme BETA (|BETA| < 10)
    df = df[np.abs(df['BETA']) < 10]
    stats['n_after_beta'] = len(df)

    # Filter P in (0, 1]
    df = df[(df['P'] > 0) & (df['P'] <= 1)]
    stats['n_final'] = len(df)

    # Flag palindromic SNPs
    # Cast EA/OA to string first — columns may be NaN-filled object dtype
    pal = {('A','T'),('T','A'),('C','G'),('G','C')}
    ea_str = df['EA'].astype(str).str.upper()
    oa_str = df['OA'].astype(str).str.upper()
    df = df.assign(EA_u=ea_str, OA_u=oa_str)
    df['palindromic'] = df.apply(lambda r: (r['EA_u'], r['OA_u']) in pal, axis=1)
    stats['n_palindromic'] = int(df['palindromic'].sum())
    df = df.drop(columns=['EA_u','OA_u'])

    # Genomic lambda
    chi2 = (df['BETA'] / df['SE'])**2
    stats['lambda_gc'] = float(np.round(np.median(chi2) / 0.4549, 4))

    print(f"[validate] {cohort}: {stats['n_raw']} → {stats['n_final']} variants  "
          f"(lambda={stats['lambda_gc']})")
    return df, stats


def main():
    args = parse_args()
    reports = []

    for cf in args.cohort:
        cohort, fpath = cf.split(':', 1)
        df, stats = qc_sumstats(load_sumstats(fpath, cohort), cohort)
        out = f"{args.out_prefix}_{cohort}.tsv.gz"
        df.to_csv(out, sep='\t', index=False, compression='gzip')
        print(f"[validate] Written: {out}")
        reports.append(stats)

    rpt = pd.DataFrame(reports)
    rpt.to_csv(args.report, sep='\t', index=False)
    print(f"[validate] Report: {args.report}")


if __name__ == '__main__':
    main()
