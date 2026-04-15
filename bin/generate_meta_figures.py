#!/usr/bin/env python3
"""
generate_meta_figures.py

Author: Nadeem Khan, INRS-CAFSB
"""

import argparse
import glob
import os
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D

# ── Typography & palette ──────────────────────────────────────
plt.rcParams.update({
    'font.family':       'DejaVu Sans',
    'font.size':         11,
    'axes.titlesize':    13,
    'axes.titleweight':  'normal',
    'axes.labelsize':    11,
    'axes.spines.top':   False,
    'axes.spines.right': False,
    'axes.grid':         False,
    'figure.facecolor':  'white',
    'axes.facecolor':    'white',
    'savefig.facecolor': 'white',
    'savefig.dpi':       300,
    'savefig.bbox':      'tight',
    'legend.frameon':    True,
    'legend.framealpha': 1.0,
    'legend.edgecolor':  '#e2e8f0',
    'legend.fontsize':   9,
})

CHR_EVEN       = '#3B82F6'
CHR_ODD        = '#6366F1'
SIG_RED        = '#DC2626'
SUG_ORG        = '#EA580C'
COHORT_PALETTE = ['#2563EB', '#7C3AED', '#059669', '#D97706', '#DC2626', '#0891B2']
META_COLOR     = '#0f172a'


# ── Utilities ─────────────────────────────────────────────────
def ensure_dir(p):
    os.makedirs(p, exist_ok=True)
    return p


def save(fig, path_png):
    fig.savefig(path_png, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'[figures] Written: {path_png}')


def get_p_col(df):
    for c in ['p_re_gc', 'p_re', 'p_fe', 'P']:
        if c in df.columns:
            return c
    return None


def chr_offsets(df):
    df = df.copy()
    df['CHR'] = pd.to_numeric(df['CHR'], errors='coerce')
    df = df.dropna(subset=['CHR', 'POS'])
    df['CHR'] = df['CHR'].astype(int)
    chrom_max = df.groupby('CHR')['POS'].max().sort_index()
    offsets, cum = {}, 0
    for ch in range(1, 23):
        offsets[ch] = cum
        cum += int(chrom_max.get(ch, 250_000_000)) + 8_000_000
    df['x_pos'] = df.apply(lambda r: r['POS'] + offsets.get(r['CHR'], 0), axis=1)
    return df, offsets


def load_meta_files(d):
    out = {}
    for f in glob.glob(os.path.join(d, '**', '*_meta_results.tsv.gz'), recursive=True):
        t = os.path.basename(f).replace('_meta_results.tsv.gz', '')
        try:
            out[t] = pd.read_csv(f, sep='\t', low_memory=False, compression='gzip')
        except Exception as e:
            print(f'[figures] Warning: {f}: {e}')
    return out


def load_qq_files(d):
    out = {}
    for f in glob.glob(os.path.join(d, '**', '*_qq_data.tsv.gz'), recursive=True):
        t = os.path.basename(f).replace('_qq_data.tsv.gz', '')
        try:
            out[t] = pd.read_csv(f, sep='\t', compression='gzip')
        except Exception as e:
            print(f'[figures] Warning: {f}: {e}')
    return out


def load_het_files(d):
    out = {}
    for f in glob.glob(os.path.join(d, '**', '*_het_stats.tsv.gz'), recursive=True):
        t = os.path.basename(f).replace('_het_stats.tsv.gz', '')
        try:
            out[t] = pd.read_csv(f, sep='\t', compression='gzip')
        except Exception as e:
            print(f'[figures] Warning: {f}: {e}')
    return out


def load_loci_files(d):
    out = {}
    for f in glob.glob(os.path.join(d, '**', '*_annotated_loci.tsv'), recursive=True):
        t = os.path.basename(f).replace('_annotated_loci.tsv', '')
        try:
            out[t] = pd.read_csv(f, sep='\t', low_memory=False)
        except Exception as e:
            print(f'[figures] Warning: {f}: {e}')
    return out


# ═══════════════════════════════════════════════════════════════
# 1. MANHATTAN PLOT
# ═══════════════════════════════════════════════════════════════
def make_manhattan(df, trait, p_sig, p_sug, out_dir):
    p_col = get_p_col(df)
    if p_col is None:
        print(f'[figures] Manhattan: no P column for {trait}')
        return

    df = df[df[p_col] > 0].dropna(subset=[p_col, 'CHR', 'POS']).copy()
    df, offsets = chr_offsets(df)
    df['logp'] = -np.log10(df[p_col])

    # Thin non-significant points
    bg  = df[df[p_col] >= p_sug]
    sug = df[(df[p_col] < p_sug) & (df[p_col] >= p_sig)]
    sig = df[df[p_col] < p_sig]
    if len(bg) > 100_000:
        bg = bg.sample(100_000, random_state=42)

    fig, ax = plt.subplots(figsize=(16, 4.5))

    # Background points by chromosome colour
    for ch in df['CHR'].unique():
        col = CHR_EVEN if ch % 2 == 0 else CHR_ODD
        sub = bg[bg['CHR'] == ch]
        if not sub.empty:
            ax.scatter(sub['x_pos'], sub['logp'], c=col, s=2.5,
                       alpha=0.60, linewidths=0, rasterized=True)

    if not sug.empty:
        ax.scatter(sug['x_pos'], sug['logp'], c=SUG_ORG, s=5,
                   alpha=0.88, linewidths=0, zorder=3,
                   label=f'Suggestive  P < {p_sug:.0e}')
    if not sig.empty:
        ax.scatter(sig['x_pos'], sig['logp'], c=SIG_RED, s=7,
                   alpha=0.92, linewidths=0, zorder=4,
                   label=f'Significant  P < {p_sig:.0e}')

    # Threshold lines
    ax.axhline(-np.log10(p_sig), color=SIG_RED, lw=1.2, ls=':', zorder=2)
    ax.axhline(-np.log10(p_sug), color=SUG_ORG, lw=1.0, ls='--', zorder=2)

    # Chromosome ticks
    chroms = sorted(offsets.keys())
    tick_x, tick_l = [], []
    for i, ch in enumerate(chroms):
        nxt  = chroms[i + 1] if i + 1 < len(chroms) else ch
        mid  = (offsets[ch] + offsets.get(nxt, offsets[ch] + 250_000_000)) / 2
        tick_x.append(mid)
        tick_l.append(str(ch))
    ax.set_xticks(tick_x)
    ax.set_xticklabels(tick_l, fontsize=8)
    ax.set_xlim(0, max(offsets.values()) + 250_000_000)
    ax.set_xlabel('Chromosome', labelpad=8)
    ax.set_ylabel('−log₁₀(P)', labelpad=8)
    ax.set_title(f'{trait.replace("_", " ").title()}  ·  Meta-analysis Manhattan',
                 pad=12, loc='left')

    if not sug.empty or not sig.empty:
        ax.legend(loc='upper right', bbox_to_anchor=(1.01, 1),
                  borderaxespad=0, framealpha=1)

    fig.tight_layout()
    ensure_dir(os.path.join(out_dir, 'manhattan'))
    save(fig, os.path.join(out_dir, 'manhattan', f'{trait}_manhattan.png'))


# ═══════════════════════════════════════════════════════════════
# 2. QQ PLOT
# ═══════════════════════════════════════════════════════════════
def make_qq(qq_data, trait, out_dir):
    if qq_data is None or qq_data.empty:
        return
    df = qq_data.dropna()
    if df.empty:
        return

    fig, ax = plt.subplots(figsize=(5, 5))

    max_val = max(float(df['expected'].max()), float(df['observed'].max())) + 0.3

    # 95 % CI band
    n     = len(df)
    n_ci  = min(n, 5000)
    ranks = np.arange(1, n_ci + 1)
    x_ci  = -np.log10(ranks / (n + 1))[::-1]
    ci_up = -np.log10(np.array([np.percentile(
        np.random.beta(r, n - r + 1, 500), 97.5) for r in ranks]))
    ci_lo = -np.log10(np.array([np.percentile(
        np.random.beta(r, n - r + 1, 500), 2.5)  for r in ranks]))
    ax.fill_between(x_ci, ci_lo[:len(x_ci)], ci_up[:len(x_ci)],
                    color='#3B82F6', alpha=0.12, label='95% CI')

    # Diagonal
    ax.plot([0, max_val], [0, max_val], color='#94a3b8', lw=1.4, ls='--',
            label='Expected (null)')

    # Points
    ax.scatter(df['expected'], df['observed'], c='#3B82F6', s=4,
               alpha=0.60, linewidths=0, label='Observed')

    ax.set_xlabel('Expected −log₁₀(P)', labelpad=8)
    ax.set_ylabel('Observed −log₁₀(P)', labelpad=8)
    ax.set_title(f'{trait.replace("_", " ").title()}  ·  QQ Plot',
                 pad=12, loc='left')
    ax.set_xlim(0, max_val)
    ax.set_ylim(0, max_val)
    ax.legend(loc='upper left', framealpha=1)

    fig.tight_layout()
    ensure_dir(os.path.join(out_dir, 'qqplot'))
    save(fig, os.path.join(out_dir, 'qqplot', f'{trait}_qqplot.png'))


# ═══════════════════════════════════════════════════════════════
# 3. FOREST PLOT  (top 5 lead SNPs)
# ═══════════════════════════════════════════════════════════════
def make_forest(loci_df, trait, out_dir):
    if loci_df is None or loci_df.empty:
        return
    beta_cols = sorted([c for c in loci_df.columns if c.startswith('BETA_')])
    se_cols   = sorted([c for c in loci_df.columns if c.startswith('SE_')])
    cohorts   = [c.replace('BETA_', '') for c in beta_cols]
    if not cohorts:
        return

    p_col = get_p_col(loci_df)
    top   = loci_df.nsmallest(min(5, len(loci_df)), p_col) if p_col else loci_df.head(5)
    if top.empty:
        return

    # Build row list
    rows = []
    for _, snp_row in top.iterrows():
        snp = str(snp_row.get('SNP', ''))
        rows.append(dict(label=snp, beta=None, lo=None, hi=None,
                         color=None, is_header=True))
        for j, (coh, bc, sc) in enumerate(zip(cohorts, beta_cols, se_cols)):
            b = snp_row.get(bc, np.nan)
            s = snp_row.get(sc, np.nan)
            if pd.isna(b) or pd.isna(s):
                continue
            rows.append(dict(label=f'  {coh}', beta=b, lo=b - 1.96*s,
                             hi=b + 1.96*s,
                             color=COHORT_PALETTE[j % len(COHORT_PALETTE)],
                             is_header=False, is_meta=False, coh=coh))
        mb = snp_row.get('beta_re', np.nan)
        ms = snp_row.get('se_re',   np.nan)
        if not (pd.isna(mb) or pd.isna(ms)):
            rows.append(dict(label='  ◆ Meta RE', beta=mb, lo=mb - 1.96*ms,
                             hi=mb + 1.96*ms, color=META_COLOR,
                             is_header=False, is_meta=True, coh='Meta RE'))

    n_rows = len(rows)
    fig_h  = max(4.0, n_rows * 0.38)
    fig, ax = plt.subplots(figsize=(8, fig_h))

    ax.axvline(0, color='#94a3b8', lw=1.0, ls=':', zorder=0)

    y_pos    = n_rows - 1
    y_ticks  = []
    y_labels = []
    legend_seen = set()
    legend_handles = []

    for r in rows:
        if r['is_header']:
            ax.text(-0.02, y_pos, r['label'], transform=ax.get_yaxis_transform(),
                    ha='right', va='center', fontsize=10, fontweight='bold',
                    color='#0f172a')
            y_ticks.append(y_pos)
            y_labels.append('')
        else:
            ax.plot([r['lo'], r['hi']], [y_pos, y_pos],
                    color=r['color'], lw=1.8, solid_capstyle='round')
            mk = 'D' if r.get('is_meta') else 's'
            ms = 8  if r.get('is_meta') else 6
            ax.plot(r['beta'], y_pos, marker=mk, ms=ms,
                    color=r['color'], zorder=3)
            y_ticks.append(y_pos)
            y_labels.append(r['label'])
            if r['coh'] not in legend_seen:
                legend_seen.add(r['coh'])
                legend_handles.append(
                    Line2D([0], [0], marker='D' if r.get('is_meta') else 's',
                           color=r['color'], ms=6, lw=0, label=r['coh']))
        y_pos -= 1

    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels, fontsize=9)
    ax.set_xlabel('Effect size (β)', labelpad=8)
    ax.set_title(f'{trait.replace("_", " ").title()}  ·  Forest Plot',
                 pad=12, loc='left')
    ax.spines['left'].set_visible(False)

    if legend_handles:
        ax.legend(handles=legend_handles, loc='upper right',
                  bbox_to_anchor=(1.18, 1), borderaxespad=0, framealpha=1)

    fig.tight_layout()
    ensure_dir(os.path.join(out_dir, 'forest'))
    save(fig, os.path.join(out_dir, 'forest', f'{trait}_forest.png'))


# ═══════════════════════════════════════════════════════════════
# 4. HETEROGENEITY HEATMAP
# ═══════════════════════════════════════════════════════════════
def make_het_heatmap(het_dict, p_sig, out_dir):
    if not het_dict:
        return
    rows = []
    for trait, df in het_dict.items():
        p_col = get_p_col(df)
        sub   = df[df[p_col] < p_sig].head(20) if (p_col and p_col in df.columns) \
                else df.head(20)
        for _, row in sub.iterrows():
            rows.append({'trait': trait.replace('_', ' ').title(),
                         'SNP':   str(row.get('SNP', '')),
                         'i2':    row.get('i2', np.nan)})
    if not rows:
        return
    df_h = pd.DataFrame(rows).dropna(subset=['i2'])
    if df_h.empty:
        return
    pivot = df_h.pivot_table(index='trait', columns='SNP',
                              values='i2', aggfunc='mean')
    pivot = pivot.dropna(axis=1, how='all').dropna(axis=0, how='all').iloc[:, :25]
    if pivot.empty:
        return

    n_row, n_col = pivot.shape
    fig, ax = plt.subplots(figsize=(max(6, n_col * 0.8 + 2),
                                    max(3, n_row * 0.7 + 1.5)))

    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list(
        'het', ['#f0fdf4', '#fef9c3', '#fee2e2'])

    im = ax.imshow(pivot.values, aspect='auto', cmap=cmap, vmin=0, vmax=100)

    ax.set_xticks(range(n_col))
    ax.set_xticklabels([s[:14] for s in pivot.columns],
                       rotation=45, ha='right', fontsize=8)
    ax.set_yticks(range(n_row))
    ax.set_yticklabels(pivot.index, fontsize=9)
    ax.set_xlabel('Lead SNP', labelpad=8)
    ax.set_title('Heterogeneity  ·  I² across traits and loci', pad=12, loc='left')

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label('I² (%)', fontsize=10)
    cbar.set_ticks([0, 25, 50, 75, 100])

    # Annotate cells
    for i in range(n_row):
        for j in range(n_col):
            v = pivot.values[i, j]
            if not np.isnan(v):
                ax.text(j, i, f'{v:.0f}', ha='center', va='center',
                        fontsize=7, color='#374151')

    fig.tight_layout()
    ensure_dir(os.path.join(out_dir, 'global'))
    save(fig, os.path.join(out_dir, 'global', 'heterogeneity_heatmap.png'))


# ═══════════════════════════════════════════════════════════════
# 5. EFFECT-SIZE SCATTER  (Cohort A vs Cohort B)
# ═══════════════════════════════════════════════════════════════
def make_effect_scatter(meta_dict, p_sig, out_dir):
    for trait, df in meta_dict.items():
        b_cols = [c for c in df.columns if c.startswith('BETA_')]
        if len(b_cols) < 2:
            continue
        c1, c2 = b_cols[0], b_cols[1]
        n1, n2 = c1.replace('BETA_', ''), c2.replace('BETA_', '')
        p_col  = get_p_col(df)
        sub    = df.dropna(subset=[c1, c2])
        if sub.empty:
            continue
        lim = float(max(sub[c1].abs().max(), sub[c2].abs().max())) * 1.10

        if p_col and p_col in sub.columns:
            sig    = sub[sub[p_col] < p_sig]
            nonsig = sub[sub[p_col] >= p_sig].sample(
                min(50_000, int((sub[p_col] >= p_sig).sum())), random_state=42)
        else:
            sig    = pd.DataFrame()
            nonsig = sub.sample(min(50_000, len(sub)), random_state=42)

        fig, ax = plt.subplots(figsize=(5.5, 5.5))

        ax.plot([-lim, lim], [-lim, lim], color='#94a3b8', lw=1.2, ls='--',
                label='Identity (β₁ = β₂)', zorder=1)
        ax.scatter(nonsig[c1], nonsig[c2], c='#93c5fd', s=3, alpha=0.30,
                   linewidths=0, rasterized=True, label='Non-significant', zorder=2)
        if not sig.empty:
            ax.scatter(sig[c1], sig[c2], c=SIG_RED, s=6, alpha=0.85,
                       linewidths=0, zorder=3,
                       label=f'GW significant  P < {p_sig:.0e}')

        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xlabel(f'β  {n1}', labelpad=8)
        ax.set_ylabel(f'β  {n2}', labelpad=8)
        ax.set_title(f'{trait.replace("_", " ").title()}  ·  Effect Concordance',
                     pad=12, loc='left')
        ax.legend(loc='upper left', framealpha=1)

        fig.tight_layout()
        ensure_dir(os.path.join(out_dir, 'scatter'))
        save(fig, os.path.join(out_dir, 'scatter', f'{trait}_effect_scatter.png'))


# ═══════════════════════════════════════════════════════════════
# 6. LAMBDA SUMMARY BAR
# ═══════════════════════════════════════════════════════════════
def make_lambda_bar(meta_dict, out_dir):
    rows = []
    for trait, df in meta_dict.items():
        if 'BETA' in df.columns and 'SE' in df.columns:
            v    = df.dropna(subset=['BETA', 'SE'])
            v    = v[v['SE'] > 0]
            chi2 = (v['BETA'] / v['SE']) ** 2
        else:
            p_col = get_p_col(df)
            if p_col is None:
                continue
            p = df[p_col].dropna()
            p = p[(p > 0) & (p < 1)]
            if len(p) < 100:
                continue
            from scipy.stats import norm
            chi2 = pd.Series(norm.ppf(p.values / 2) ** 2)
        lam = float(chi2.dropna().median() / 0.4549)
        rows.append({'trait': trait.replace('_', ' ').title(), 'lambda': round(lam, 4)})

    if not rows:
        return
    df_l = pd.DataFrame(rows).sort_values('lambda', ascending=False)

    fig, ax = plt.subplots(figsize=(max(5, len(df_l) * 0.9), 4))

    colours = [SIG_RED if v > 1.1 else '#3B82F6' for v in df_l['lambda']]
    bars = ax.bar(df_l['trait'], df_l['lambda'], color=colours, width=0.55,
                  edgecolor='none')

    for bar, val in zip(bars, df_l['lambda']):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.01, f'{val:.3f}',
                ha='center', va='bottom', fontsize=9)

    ax.axhline(1.0, color='#94a3b8', lw=1.2, ls=':', label='λ = 1.0')
    ax.axhline(1.1, color=SIG_RED,   lw=1.0, ls='--', label='λ = 1.1 threshold')

    ax.set_ylim(0, max(float(df_l['lambda'].max()) * 1.22, 1.28))
    ax.set_ylabel('λ GC', labelpad=8)
    ax.set_title('Genomic Inflation  ·  λ GC per trait', pad=12, loc='left')
    ax.tick_params(axis='x', rotation=20)
    ax.legend(loc='upper right', framealpha=1)

    fig.tight_layout()
    ensure_dir(os.path.join(out_dir, 'global'))
    save(fig, os.path.join(out_dir, 'global', 'lambda_summary.png'))


# ═══════════════════════════════════════════════════════════════
# CLI + main
# ═══════════════════════════════════════════════════════════════
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--meta_dir',  required=True)
    p.add_argument('--het_dir',   required=True)
    p.add_argument('--qq_dir',    required=True)
    p.add_argument('--loci_dir',  required=True)
    p.add_argument('--author',    default='')
    p.add_argument('--institute', default='')
    p.add_argument('--p_sig',     type=float, default=5e-8)
    p.add_argument('--p_sug',     type=float, default=1e-5)
    p.add_argument('--i2_thresh', type=float, default=75)
    p.add_argument('--out_dir',   default='.')
    return p.parse_args()


def main():
    args = parse_args()
    print('[figures] Loading data files...')
    meta_dict = load_meta_files(args.meta_dir)
    qq_dict   = load_qq_files(args.qq_dir)
    het_dict  = load_het_files(args.het_dir)
    loci_dict = load_loci_files(args.loci_dir)
    traits    = sorted(meta_dict.keys())
    print(f'[figures] Traits found: {traits}')

    for trait in traits:
        print(f'\n[figures] Processing trait: {trait}')
        df = meta_dict[trait]
        make_manhattan(df, trait, args.p_sig, args.p_sug, args.out_dir)
        make_qq(qq_dict.get(trait), trait, args.out_dir)
        make_forest(loci_dict.get(trait), trait, args.out_dir)
        make_effect_scatter({trait: df}, args.p_sig, args.out_dir)

    print('\n[figures] Generating global figures...')
    make_het_heatmap(het_dict, args.p_sig, args.out_dir)
    make_lambda_bar(meta_dict, args.out_dir)
    print('\n[figures] All figures complete.')


if __name__ == '__main__':
    main()
