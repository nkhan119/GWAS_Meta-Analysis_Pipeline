#!/usr/bin/env python3
"""
harmonise_meta.py
Aligns alleles across cohorts, filters to common variants,
applies genomic-control correction if lambda > 1.1.
Author: Nadeem Khan, INRS-CAFSB
"""

import argparse
import sys
import pandas as pd
import numpy as np
from scipy import stats


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--trait',       required=True)
    p.add_argument('--type',        default='qt')
    p.add_argument('--input',       action='append', default=[],
                   help='cohort:filepath, repeatable')
    p.add_argument('--min_cohorts', type=int, default=2)
    p.add_argument('--ref_build',   default='GRCh37')
    p.add_argument('--out',         required=True)
    p.add_argument('--lambda_out',  required=True)
    return p.parse_args()


COMPLEMENT = {'A':'T','T':'A','C':'G','G':'C'}
PALINDROMIC = {('A','T'),('T','A'),('C','G'),('G','C')}


def complement_flip(df, cohort):
    """Flip alleles to reference orientation."""
    flip = (df['EA_ref'] != df[f'EA_{cohort}']) & \
           (df['EA_ref'] == df[f'OA_{cohort}'])
    df.loc[flip, f'BETA_{cohort}'] = -df.loc[flip, f'BETA_{cohort}']
    df.loc[flip, f'EAF_{cohort}']  =  1 - df.loc[flip, f'EAF_{cohort}'].fillna(0.5)
    return df


def load_cohort(cf_str):
    cohort, fpath = cf_str.split(':', 1)
    print(f"[harmonise] Loading {cohort}: {fpath}")
    df = pd.read_csv(fpath, sep='\t', low_memory=False,
                     compression='gzip' if fpath.endswith('.gz') else 'infer')
    df.columns = [c.upper() for c in df.columns]

    # Alias non-standard column names (mirrors validate_sumstats.py)
    ALIASES = {
        'BP': 'POS', 'BASE_PAIR_LOCATION': 'POS',
        'A1': 'EA',  'ALLELE1': 'EA', 'ALT': 'EA',
        'A2': 'OA',  'ALLELE0': 'OA', 'REF': 'OA',
        'PVALUE': 'P', 'P_BOLT_LMM': 'P', 'P_BOLT_LMM_INF': 'P',
        'OR': 'BETA', 'Z': 'BETA',
        'SE_BOLT_LMM': 'SE',
        'N_OBS': 'N', 'NOBS': 'N', 'OBS_CT': 'N',
        'RSID': 'SNP', 'ID': 'SNP', 'SNPID': 'SNP', 'MARKERNAME': 'SNP',
        'A1FREQ': 'EAF', 'MAF': 'EAF', 'AF1': 'EAF',
    }
    for src, tgt in ALIASES.items():
        if src in df.columns and tgt not in df.columns:
            df = df.rename(columns={src: tgt})
    if 'P' not in df.columns and 'LOG10P' in df.columns:
        df['P'] = 10.0 ** (-df['LOG10P'].astype(float))

    # Rename to cohort-specific columns
    rename = {
        'BETA': f'BETA_{cohort}',
        'SE':   f'SE_{cohort}',
        'P':    f'P_{cohort}',
        'EAF':  f'EAF_{cohort}',
        'N':    f'N_{cohort}',
        'EA':   f'EA_{cohort}',
        'OA':   f'OA_{cohort}',
    }
    df = df.rename(columns=rename)
    return cohort, df


def gc_correct(df, cohort):
    beta_col = f'BETA_{cohort}'
    se_col   = f'SE_{cohort}'
    if beta_col not in df.columns or se_col not in df.columns:
        return df, 1.0
    z2     = (df[beta_col] / df[se_col]) ** 2
    lam    = float(np.median(z2.dropna()) / 0.4549)
    if lam > 1.1:
        print(f"[harmonise] GC correction {cohort}: lambda={lam:.4f}")
        df[se_col] = df[se_col] * np.sqrt(lam)
    return df, round(lam, 4)


def main():
    args = parse_args()
    cohort_dfs = {}
    for cf in args.input:
        cohort, df = load_cohort(cf)
        cohort_dfs[cohort] = df
    cohorts = list(cohort_dfs.keys())

    # ── Merge on SNP ──────────────────────────────────────────
    base_cols = ['SNP', 'CHR', 'POS']
    first_cohort = cohorts[0]
    merged = cohort_dfs[first_cohort][
        base_cols + [f'EA_{first_cohort}', f'OA_{first_cohort}',
                     f'BETA_{first_cohort}', f'SE_{first_cohort}',
                     f'P_{first_cohort}', f'EAF_{first_cohort}', f'N_{first_cohort}']
    ].copy()
    merged = merged.rename(columns={f'EA_{first_cohort}': 'EA_ref',
                                     f'OA_{first_cohort}': 'OA_ref'})

    for cohort in cohorts[1:]:
        df = cohort_dfs[cohort]
        keep = ['SNP', f'EA_{cohort}', f'OA_{cohort}',
                f'BETA_{cohort}', f'SE_{cohort}',
                f'P_{cohort}', f'EAF_{cohort}', f'N_{cohort}']
        keep = [c for c in keep if c in df.columns]
        merged = merged.merge(df[keep], on='SNP', how='inner')
        merged = complement_flip(merged, cohort)

    print(f"[harmonise] Common variants after merge: {len(merged)}")

    # ── Remove palindromic SNPs ────────────────────────────────
    pal = merged.apply(
        lambda r: (str(r['EA_ref']).upper(), str(r['OA_ref']).upper()) in PALINDROMIC, axis=1)
    merged = merged[~pal]
    print(f"[harmonise] After palindromic filter: {len(merged)}")

    # ── Require min_cohorts non-missing ───────────────────────
    beta_cols = [f'BETA_{c}' for c in cohorts]
    n_present = merged[beta_cols].notna().sum(axis=1)
    merged = merged[n_present >= args.min_cohorts]
    print(f"[harmonise] After min_cohorts={args.min_cohorts} filter: {len(merged)}")

    # ── GC correction + lambda report ─────────────────────────
    lambda_rows = []
    for cohort in cohorts:
        merged, lam = gc_correct(merged, cohort)
        lambda_rows.append({'cohort': cohort, 'trait': args.trait, 'lambda_gc': lam})

    lambda_df = pd.DataFrame(lambda_rows)
    lambda_df.to_csv(args.lambda_out, sep='\t', index=False)

    # ── Rename ref alleles for output ─────────────────────────
    merged = merged.rename(columns={'EA_ref': 'EA', 'OA_ref': 'OA'})

    # Average EAF across cohorts
    eaf_cols = [f'EAF_{c}' for c in cohorts if f'EAF_{c}' in merged.columns]
    if eaf_cols:
        merged['EAF'] = merged[eaf_cols].mean(axis=1)

    merged.to_csv(args.out, sep='\t', index=False, compression='gzip')
    print(f"[harmonise] Written: {args.out}")


if __name__ == '__main__':
    main()
