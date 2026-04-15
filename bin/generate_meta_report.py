#!/usr/bin/env python3
"""
generate_meta_report.py
Generates the interactive HTML report for GWAS meta-analysis.
White background, clean typography, publication-quality embedded Plotly figures.
Author: Nadeem Khan, INRS-CAFSB
"""

import argparse
import glob
import os
import json
import warnings
warnings.filterwarnings('ignore')
from datetime import datetime

import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.io as pio


# ── Design tokens ─────────────────────────────────────────────
FONT    = 'Inter, -apple-system, Arial, sans-serif'
BG      = '#ffffff'
C_BLUE  = '#2563EB'
C_IND   = '#4F46E5'
C_RED   = '#DC2626'
C_ORG   = '#EA580C'
C_GRN   = '#16A34A'
C_SLT   = '#64748B'
C_LBLUE = '#DBEAFE'

CHR_EVEN = '#3B82F6'
CHR_ODD  = '#818CF8'

BASE_LAYOUT = dict(
    paper_bgcolor = BG,
    plot_bgcolor  = BG,
    font          = dict(family=FONT, size=13, color='#0f172a'),
    margin        = dict(l=70, r=40, t=65, b=60),
    showlegend    = False,
)
AX = dict(
    showgrid   = False, zeroline  = False,
    showline   = True,  linecolor = '#e2e8f0', linewidth = 1.5,
    ticks      = 'outside', ticklen = 4, tickcolor = '#94a3b8',
    title_standoff = 12,
)


# ── Data loaders ──────────────────────────────────────────────
def glob_load(pattern, sep='\t', gz=True):
    out = {}
    for f in glob.glob(pattern, recursive=True):
        key = os.path.basename(f).split('_')[0]
        # better key: strip known suffixes
        for sfx in ['_meta_results.tsv.gz','_qq_data.tsv.gz',
                    '_het_stats.tsv.gz','_annotated_loci.tsv',
                    '_filter_summary.tsv','_meta_summary.tsv']:
            key = os.path.basename(f).replace(sfx, '')
        try:
            comp = 'gzip' if f.endswith('.gz') else 'infer'
            out[key] = pd.read_csv(f, sep=sep, low_memory=False, compression=comp)
        except Exception as e:
            print(f'[report] Warning loading {f}: {e}')
    return out


def get_p_col(df):
    for c in ['p_re_gc','p_re','p_fe','P']:
        if c in df.columns: return c
    return None


# ── Plotly figure builders ─────────────────────────────────────
def chr_offsets(df):
    df = df.copy()
    df['CHR'] = pd.to_numeric(df['CHR'], errors='coerce')
    df = df.dropna(subset=['CHR','POS'])
    df['CHR'] = df['CHR'].astype(int)
    chrom_max = df.groupby('CHR')['POS'].max().sort_index()
    offsets = {}
    cum = 0
    for chrom in range(1, 23):
        offsets[chrom] = cum
        cum += int(chrom_max.get(chrom, 3e8)) + 8_000_000
    df['x_pos'] = df.apply(lambda r: r['POS'] + offsets.get(r['CHR'],0), axis=1)
    df['chrom_col'] = df['CHR'].apply(lambda c: CHR_EVEN if c%2==0 else CHR_ODD)
    return df, offsets


def manhattan_fig(df, trait, p_sig=5e-8, p_sug=1e-5):
    p_col = get_p_col(df)
    if p_col is None or df.empty: return go.Figure()
    df = df[df[p_col] > 0].dropna(subset=[p_col,'CHR','POS'])
    df, offsets = chr_offsets(df)
    df['logp'] = -np.log10(df[p_col])
    # Subsample
    nonsig = df[df[p_col] >= p_sug]
    sig    = df[df[p_col] <  p_sug]
    if len(nonsig) > 80_000:
        nonsig = nonsig.sample(80_000, random_state=42)
    pdata = pd.concat([nonsig, sig])

    fig = go.Figure()
    for col_val, grp in pdata.groupby('chrom_col'):
        is_sig = grp[p_col] < p_sig
        is_sug = (grp[p_col] < p_sug) & ~is_sig
        for mask, c, sz in [(~is_sig & ~is_sug, col_val, 3.5),
                             (is_sug, C_ORG, 5),
                             (is_sig, C_RED, 6)]:
            sub = grp[mask]
            if sub.empty: continue
            fig.add_trace(go.Scattergl(
                x=sub['x_pos'], y=sub['logp'],
                mode='markers',
                marker=dict(color=c, size=sz, opacity=0.75),
                hovertext=sub.get('SNP','').astype(str) + '<br>' +
                          'P=' + sub[p_col].map('{:.2e}'.format),
                hoverinfo='text',
            ))
    sig_lp = -np.log10(p_sig)
    sug_lp = -np.log10(p_sug)
    fig.add_hline(y=sig_lp, line=dict(color=C_RED, width=1.2, dash='dot'))
    fig.add_hline(y=sug_lp, line=dict(color=C_ORG, width=1.0, dash='dash'))
    # Chr ticks
    tv, tt = [], []
    for ch in range(1,23):
        if ch in offsets and ch+1 in offsets:
            tv.append((offsets[ch] + offsets.get(ch+1, offsets[ch]+2e8))/2)
            tt.append(str(ch))
    fig.update_layout(**BASE_LAYOUT,
        xaxis=dict(**AX, tickvals=tv, ticktext=tt, title='Chromosome'),
        yaxis=dict(**AX, title='−log₁₀(P)'),
        title=dict(text=f'{trait.replace("_"," ").title()} — Manhattan',
                   font=dict(size=14, weight='normal'), x=0.01),
        height=360, width=None)
    return fig


def qq_fig(qq_df, trait):
    if qq_df is None or qq_df.empty: return go.Figure()
    df = qq_df.dropna()
    max_v = max(df['expected'].max(), df['observed'].max()) + 0.3
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=[0, max_v], y=[0, max_v],
        mode='lines', line=dict(color='#cbd5e1', width=1.5, dash='dash'),
        hoverinfo='skip'))
    fig.add_trace(go.Scattergl(
        x=df['expected'], y=df['observed'],
        mode='markers',
        marker=dict(color=C_BLUE, size=3.5, opacity=0.55),
        hovertemplate='Exp: %{x:.3f}<br>Obs: %{y:.3f}<extra></extra>'))
    fig.update_layout(**BASE_LAYOUT,
        xaxis=dict(**AX, title='Expected −log₁₀(P)'),
        yaxis=dict(**AX, title='Observed −log₁₀(P)'),
        title=dict(text=f'{trait.replace("_"," ").title()} — QQ',
                   font=dict(size=14, weight='normal'), x=0.01),
        height=340, width=None)
    return fig


def forest_fig(loci_df, trait):
    if loci_df is None or loci_df.empty: return go.Figure()
    beta_cols = sorted([c for c in loci_df.columns if c.startswith('BETA_')])
    se_cols   = sorted([c for c in loci_df.columns if c.startswith('SE_')])
    cohorts   = [c.replace('BETA_','') for c in beta_cols]
    if not cohorts: return go.Figure()

    p_col = get_p_col(loci_df)
    top   = loci_df.nsmallest(min(6, len(loci_df)), p_col) if p_col else loci_df.head(6)

    rows, y_out, betas, los, his, cols = [], [], [], [], [], []
    pal = [C_BLUE, C_IND, C_GRN, C_ORG, C_RED]
    y_cursor = 0
    for _, snp_row in top.iterrows():
        snp = str(snp_row.get('SNP',''))
        rows.append(dict(y=y_cursor, label=f'<b>{snp}</b>', is_header=True))
        y_cursor -= 1
        for j,(coh,bc,sc) in enumerate(zip(cohorts,beta_cols,se_cols)):
            b = snp_row.get(bc, np.nan)
            s = snp_row.get(sc, np.nan)
            if pd.isna(b) or pd.isna(s): continue
            rows.append(dict(y=y_cursor, label=f'  {coh}',
                              beta=b, lo=b-1.96*s, hi=b+1.96*s,
                              col=pal[j % len(pal)], is_header=False))
            y_cursor -= 0.7
        mb = snp_row.get('beta_re', np.nan)
        ms = snp_row.get('se_re', np.nan)
        if not pd.isna(mb) and not pd.isna(ms):
            rows.append(dict(y=y_cursor, label='  Meta RE',
                              beta=mb, lo=mb-1.96*ms, hi=mb+1.96*ms,
                              col='#0f172a', is_header=False))
        y_cursor -= 1.5

    fig = go.Figure()
    fig.add_vline(x=0, line=dict(color='#e2e8f0', width=1, dash='dot'))
    for r in rows:
        if r.get('is_header'): continue
        fig.add_trace(go.Scatter(
            x=[r['lo'], r['hi']], y=[r['y'], r['y']],
            mode='lines', line=dict(color=r['col'], width=1.8),
            hoverinfo='skip'))
        sz = 12 if r['label'].strip()=='Meta RE' else 9
        fig.add_trace(go.Scatter(
            x=[r['beta']], y=[r['y']],
            mode='markers',
            marker=dict(symbol='square', size=sz, color=r['col']),
            hovertemplate=f'{r["label"].strip()}: β={r["beta"]:.4f} [{r["lo"]:.4f},{r["hi"]:.4f}]<extra></extra>'))
    y_ticks = [r['y'] for r in rows]
    y_labels = [r['label'] for r in rows]
    fig.update_layout(**BASE_LAYOUT,
        xaxis=dict(**AX, title='Effect size (β)', zeroline=True,
                   zerolinecolor='#e2e8f0', zerolinewidth=1),
        yaxis=dict(showgrid=False, zeroline=False, showline=False,
                   tickvals=y_ticks, ticktext=y_labels,
                   tickfont=dict(size=11)),
        title=dict(text=f'{trait.replace("_"," ").title()} — Forest Plot',
                   font=dict(size=14, weight='normal'), x=0.01),
        height=max(300, abs(y_cursor)*22 + 80), width=None)
    return fig


def het_heatmap_fig(het_dict, p_sig):
    rows = []
    for trait, df in het_dict.items():
        p_col = get_p_col(df)
        sig = df[df[p_col] < p_sig] if p_col and p_col in df.columns \
              else df.head(15)
        for _, r in sig.head(15).iterrows():
            rows.append({'trait': trait.replace('_',' ').title(),
                         'SNP':   r.get('SNP',''),
                         'i2':    r.get('i2', np.nan)})
    if not rows: return go.Figure()
    df_h = pd.DataFrame(rows).dropna(subset=['i2'])
    if df_h.empty: return go.Figure()
    pivot = df_h.pivot_table(index='trait',columns='SNP',values='i2',aggfunc='mean')
    pivot = pivot.dropna(axis=1,how='all').dropna(axis=0,how='all').iloc[:,:25]
    if pivot.empty: return go.Figure()
    fig = go.Figure(go.Heatmap(
        z=pivot.values,
        x=[s[:14] for s in pivot.columns],
        y=pivot.index.tolist(),
        colorscale=[[0,'#f0fdf4'],[0.5,'#fef9c3'],[1,'#fee2e2']],
        zmin=0, zmax=100,
        colorbar=dict(title=dict(text='I² (%)'), tickvals=[0,25,50,75,100]),
        hovertemplate='%{y} · %{x}<br>I²=%{z:.1f}%<extra></extra>',
    ))
    fig.update_layout(**BASE_LAYOUT,
        xaxis=dict(showgrid=False, tickangle=40, tickfont=dict(size=10)),
        yaxis=dict(showgrid=False),
        title=dict(text='Heterogeneity — I² across traits',
                   font=dict(size=14, weight='normal'), x=0.01),
        height=max(280, len(pivot)*40+100),
        width=None)
    return fig


def effect_scatter_fig(df, trait, p_sig):
    b_cols = [c for c in df.columns if c.startswith('BETA_')]
    if len(b_cols) < 2: return go.Figure()
    c1, c2 = b_cols[0], b_cols[1]
    n1 = c1.replace('BETA_',''); n2 = c2.replace('BETA_','')
    p_col = get_p_col(df)
    sub = df.dropna(subset=[c1,c2])
    if p_col:
        sig = sub[sub[p_col] < p_sig]
        ns  = sub[sub[p_col] >= p_sig].sample(min(40_000,len(sub)),random_state=42)
    else:
        sig = pd.DataFrame(); ns = sub.sample(min(40_000,len(sub)),random_state=42)
    lim = max(abs(sub[c1]).max(), abs(sub[c2]).max()) * 1.08
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=[-lim,lim],y=[-lim,lim],
        mode='lines',line=dict(color='#cbd5e1',width=1.2,dash='dash'),
        hoverinfo='skip'))
    fig.add_hline(y=0,line=dict(color='#f1f5f9',width=0.8))
    fig.add_vline(x=0,line=dict(color='#f1f5f9',width=0.8))
    fig.add_trace(go.Scattergl(x=ns[c1],y=ns[c2],mode='markers',
        marker=dict(color='#93c5fd',size=3,opacity=0.35),hoverinfo='skip'))
    if not sig.empty:
        fig.add_trace(go.Scattergl(x=sig[c1],y=sig[c2],mode='markers',
            marker=dict(color=C_RED,size=5,opacity=0.8),
            hovertext=sig.get('SNP','').astype(str),hoverinfo='text'))
    fig.update_layout(**BASE_LAYOUT,
        xaxis=dict(**AX,title=f'β  {n1}',range=[-lim,lim]),
        yaxis=dict(**AX,title=f'β  {n2}',range=[-lim,lim]),
        title=dict(text=f'{trait.replace("_"," ").title()} — Effect Concordance',
                   font=dict(size=14,weight='normal'),x=0.01),
        height=370, width=None)
    return fig


def lambda_bar_fig(meta_dict):
    rows = []
    for trait, df in meta_dict.items():
        p_col = get_p_col(df)
        if p_col is None: continue
        p = df[p_col].dropna(); p = p[(p>0)&(p<1)]
        if len(p) < 100: continue
        # Correct λ GC: median(χ²) / 0.4549  where χ² = z²
        from scipy.stats import norm as _norm
        chi2_vals = pd.Series(_norm.ppf(p.values / 2) ** 2)
        lam = float(chi2_vals.dropna().median() / 0.4549)
        rows.append({'trait': trait.replace('_',' ').title(), 'lambda': round(lam, 4)})
    if not rows: return go.Figure()
    df_l = pd.DataFrame(rows).sort_values('lambda', ascending=False)
    colours = [C_RED if v > 1.1 else C_BLUE for v in df_l['lambda']]
    fig = go.Figure(go.Bar(
        x=df_l['trait'], y=df_l['lambda'],
        marker_color=colours,
        text=df_l['lambda'].map('{:.3f}'.format),
        textposition='outside',
    ))
    fig.add_hline(y=1.0,line=dict(color='#94a3b8',width=1.2,dash='dot'))
    fig.add_hline(y=1.1,line=dict(color=C_RED,width=1,dash='dash'))
    fig.update_layout(**BASE_LAYOUT,
        xaxis=dict(**AX,title=''),
        yaxis=dict(**AX,title='λ GC',
                   range=[0, max(df_l['lambda'].max()*1.18, 1.25)]),
        title=dict(text='Genomic Inflation — λ GC per trait',
                   font=dict(size=14,weight='normal'),x=0.01),
        height=340, width=None)
    return fig


# ── Figure → HTML div ─────────────────────────────────────────
def fig_div(fig, div_id):
    if fig is None or (hasattr(fig,'data') and len(fig.data)==0):
        return '<div class="empty-fig">No data available for this figure.</div>'
    return pio.to_html(fig, full_html=False, include_plotlyjs=False,
                       div_id=div_id,
                       config={'displayModeBar': True,
                               'modeBarButtonsToRemove': ['lasso2d','select2d'],
                               'toImageButtonOptions': {'format':'png','scale':3}})


# ── Summary table HTML ────────────────────────────────────────
def summary_table_html(summ_dict):
    rows = []
    for trait, df in summ_dict.items():
        for _, r in df.iterrows():
            rows.append(r.to_dict())
    if not rows:
        return '<p style="color:#64748b">No summary data found.</p>'
    df = pd.DataFrame(rows)
    # format numeric
    for c in ['top_p']:
        if c in df.columns:
            df[c] = df[c].apply(lambda x: f'{x:.2e}' if pd.notna(x) else '—')
    for c in ['lambda_gc','top_beta']:
        if c in df.columns:
            df[c] = df[c].apply(lambda x: f'{x:.4f}' if pd.notna(x) else '—')
    cols_show = [c for c in ['trait','type','n_variants','n_cohorts',
                              'n_sig','n_sug','lambda_gc',
                              'top_snp','top_p','top_beta'] if c in df.columns]
    df = df[cols_show]
    headers = {
        'trait':'Trait','type':'Type','n_variants':'Variants',
        'n_cohorts':'Cohorts','n_sig':'GW Sig','n_sug':'Suggestive',
        'lambda_gc':'λ GC','top_snp':'Lead SNP',
        'top_p':'Lead P','top_beta':'Lead β'
    }
    th = ''.join(f'<th>{headers.get(c,c)}</th>' for c in cols_show)
    trs = ''
    for _, row in df.iterrows():
        tds = ''.join(f'<td>{row[c]}</td>' for c in cols_show)
        trs += f'<tr>{tds}</tr>\n'
    return f'''
    <div class="table-wrap">
      <table class="data-table">
        <thead><tr>{th}</tr></thead>
        <tbody>{trs}</tbody>
      </table>
    </div>'''


def loci_table_html(loci_dict):
    all_rows = []
    for trait, df in loci_dict.items():
        df['trait'] = trait
        all_rows.append(df)
    if not all_rows:
        return '<p style="color:#64748b">No annotated loci available.</p>'
    df = pd.concat(all_rows, ignore_index=True)
    p_col = get_p_col(df)
    show_cols = [c for c in ['trait','SNP','CHR','POS','EA','OA','EAF',
                               p_col,'beta_re','se_re','i2',
                               'nearest_gene','consequence','direction']
                 if c and c in df.columns]
    if p_col and p_col in df.columns:
        df = df.sort_values(p_col)
    df = df[show_cols].head(200)
    if p_col in df.columns:
        df[p_col] = df[p_col].apply(lambda x: f'{x:.2e}' if pd.notna(x) else '—')
    for c in ['beta_re','se_re','EAF']:
        if c in df.columns:
            df[c] = df[c].apply(lambda x: f'{x:.4f}' if pd.notna(x) else '—')
    if 'i2' in df.columns:
        df['i2'] = df['i2'].apply(lambda x: f'{x:.1f}' if pd.notna(x) else '—')
    th = ''.join(f'<th>{c}</th>' for c in show_cols)
    trs = ''
    for _, row in df.iterrows():
        tds = ''.join(f'<td>{row[c]}</td>' for c in show_cols)
        trs += f'<tr>{tds}</tr>\n'
    return f'''
    <div class="table-wrap">
      <table class="data-table" id="loci-table">
        <thead><tr>{th}</tr></thead>
        <tbody>{trs}</tbody>
      </table>
    </div>'''


# ── Full HTML ─────────────────────────────────────────────────
def build_html(args, meta_dict, qq_dict, het_dict, loci_dict, summ_dict, filter_dict):
    timestamp = datetime.now().strftime('%d %B %Y, %H:%M')
    traits    = sorted(meta_dict.keys())

    # Build per-trait figure blocks
    trait_sections = ''
    for trait in traits:
        label = trait.replace('_',' ').title()
        df    = meta_dict.get(trait, pd.DataFrame())
        man   = fig_div(manhattan_fig(df, trait, args.p_sig, args.p_sug), f'man_{trait}')
        qq    = fig_div(qq_fig(qq_dict.get(trait), trait), f'qq_{trait}')
        fst   = fig_div(forest_fig(loci_dict.get(trait), trait), f'fst_{trait}')
        esc   = fig_div(effect_scatter_fig(df, trait, args.p_sig), f'esc_{trait}')

        # Mini stats
        p_col = get_p_col(df)
        n_sig = int((df[p_col] < args.p_sig).sum()) if p_col and not df.empty else 0
        n_sug = int(((df[p_col] < args.p_sug) & (df[p_col] >= args.p_sig)).sum()) \
                if p_col and not df.empty else 0
        n_var = len(df)

        trait_sections += f'''
      <section class="trait-section" id="trait-{trait}">
        <h2 class="trait-heading">{label}</h2>
        <div class="stat-row">
          <div class="stat-card"><span class="stat-num">{n_var:,}</span><span class="stat-lbl">Variants</span></div>
          <div class="stat-card accent-red"><span class="stat-num">{n_sig:,}</span><span class="stat-lbl">GW Significant</span></div>
          <div class="stat-card accent-org"><span class="stat-num">{n_sug:,}</span><span class="stat-lbl">Suggestive</span></div>
        </div>
        <div class="fig-grid fig-grid-wide">
          <div class="fig-card span-2">
            <p class="fig-label">Manhattan Plot</p>
            {man}
          </div>
        </div>
        <div class="fig-grid fig-grid-3col">
          <div class="fig-card">
            <p class="fig-label">QQ Plot</p>
            {qq}
          </div>
          <div class="fig-card">
            <p class="fig-label">Forest Plot — Lead SNPs</p>
            {fst}
          </div>
          <div class="fig-card">
            <p class="fig-label">Effect Concordance</p>
            {esc}
          </div>
        </div>
      </section>'''

    # Global figures
    het_fig  = fig_div(het_heatmap_fig(het_dict, args.p_sig), 'het_heatmap')
    lam_fig  = fig_div(lambda_bar_fig(meta_dict), 'lambda_bar')
    summ_tbl = summary_table_html(summ_dict)
    loci_tbl = loci_table_html(loci_dict)

    # Sidebar nav
    nav_links = '\n'.join(
        f'<a href="#trait-{t}" class="nav-link">{t.replace("_"," ").title()}</a>'
        for t in traits
    )

    n_cohorts_global = 2  # from config
    n_traits_total   = len(traits)
    n_sig_total = 0
    for t in traits:
        df_t  = meta_dict.get(t, pd.DataFrame())
        p_col = get_p_col(df_t)
        if p_col and not df_t.empty:
            n_sig_total += int((df_t[p_col] < args.p_sig).sum())

    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>GWAS Meta-Analysis Report — CDC 1.0.0</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>
/* ── Reset & Base ── */
*, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
html {{ scroll-behavior: smooth; }}
body {{
  font-family: {FONT};
  font-size: 14px;
  line-height: 1.65;
  color: #0f172a;
  background: #ffffff;
  display: flex;
  min-height: 100vh;
}}

/* ── Sidebar ── */
.sidebar {{
  position: fixed;
  top: 0; left: 0;
  width: 220px;
  height: 100vh;
  background: #f8fafc;
  border-right: 1px solid #e2e8f0;
  padding: 28px 0 20px;
  overflow-y: auto;
  z-index: 100;
}}
.sidebar-logo {{
  padding: 0 20px 24px;
  font-size: 12px;
  font-weight: 600;
  letter-spacing: 0.08em;
  text-transform: uppercase;
  color: {C_SLT};
  border-bottom: 1px solid #e2e8f0;
  margin-bottom: 16px;
}}
.sidebar-logo span {{ display: block; font-size: 16px; font-weight: 700;
  color: #0f172a; letter-spacing: 0; text-transform: none; margin-bottom: 2px; }}
.nav-section {{ padding: 8px 20px 4px;
  font-size: 10.5px; font-weight: 600; letter-spacing: 0.06em;
  text-transform: uppercase; color: #94a3b8; }}
.nav-link {{
  display: block; padding: 6px 20px; color: #334155;
  text-decoration: none; font-size: 13px; border-radius: 0;
  transition: background 0.15s, color 0.15s;
}}
.nav-link:hover {{ background: #eff6ff; color: {C_BLUE}; }}
.nav-link.active {{ background: #eff6ff; color: {C_BLUE}; font-weight: 500; }}

/* ── Main content ── */
.main {{
  margin-left: 220px;
  padding: 0 52px 80px;
  max-width: 1420px;
  width: 100%;
}}

/* ── Header ── */
.page-header {{
  border-bottom: 2px solid #f1f5f9;
  padding: 42px 0 28px;
  margin-bottom: 40px;
}}
.page-header h1 {{
  font-size: 24px;
  font-weight: 600;
  color: #0f172a;
  letter-spacing: -0.02em;
  margin-bottom: 6px;
}}
.page-header .meta {{
  font-size: 13px;
  color: {C_SLT};
}}
.page-header .meta span {{ margin-right: 20px; }}

/* ── Overview stats ── */
.overview-stats {{
  display: flex;
  gap: 18px;
  margin: 32px 0 44px;
  flex-wrap: wrap;
}}
.ov-card {{
  background: #f8fafc;
  border: 1px solid #e2e8f0;
  border-radius: 8px;
  padding: 20px 28px;
  min-width: 140px;
  flex: 1;
}}
.ov-card .num {{
  font-size: 28px;
  font-weight: 700;
  color: #0f172a;
  letter-spacing: -0.03em;
  line-height: 1;
  margin-bottom: 5px;
}}
.ov-card .lbl {{
  font-size: 12px;
  color: {C_SLT};
  font-weight: 500;
  text-transform: uppercase;
  letter-spacing: 0.06em;
}}
.ov-card.blue  .num {{ color: {C_BLUE}; }}
.ov-card.red   .num {{ color: {C_RED};  }}
.ov-card.green .num {{ color: {C_GRN};  }}

/* ── Section headings ── */
h2.section-heading {{
  font-size: 18px;
  font-weight: 600;
  color: #0f172a;
  margin-bottom: 20px;
  padding-bottom: 10px;
  border-bottom: 1px solid #f1f5f9;
}}
h2.trait-heading {{
  font-size: 18px;
  font-weight: 600;
  color: #0f172a;
  margin: 48px 0 18px;
  padding-bottom: 10px;
  border-bottom: 1px solid #f1f5f9;
}}

/* ── Per-trait stat row ── */
.stat-row {{
  display: flex; gap: 14px; margin-bottom: 22px; flex-wrap: wrap;
}}
.stat-card {{
  background: #f8fafc;
  border: 1px solid #e2e8f0;
  border-radius: 7px;
  padding: 14px 22px;
  display: flex; flex-direction: column; align-items: center;
}}
.stat-num {{
  font-size: 22px; font-weight: 700; color: #0f172a;
  letter-spacing: -0.02em; line-height: 1;
}}
.stat-lbl {{
  font-size: 11px; color: {C_SLT}; font-weight: 500;
  text-transform: uppercase; letter-spacing: 0.06em; margin-top: 3px;
}}
.stat-card.accent-red .stat-num  {{ color: {C_RED};  }}
.stat-card.accent-org .stat-num  {{ color: {C_ORG};  }}

/* ── Figure grid ── */
.fig-grid {{
  display: grid; gap: 20px; margin-bottom: 20px;
}}
.fig-grid-wide {{ grid-template-columns: 1fr; }}
.fig-grid-2col {{ grid-template-columns: 1fr 1fr; }}
.fig-grid-3col {{ grid-template-columns: 1fr 1fr 1fr; }}
@media (max-width: 1100px) {{
  .fig-grid-3col {{ grid-template-columns: 1fr 1fr; }}
}}
@media (max-width: 780px) {{
  .fig-grid-3col, .fig-grid-2col {{ grid-template-columns: 1fr; }}
}}
.span-2 {{ grid-column: span 2; }}

.fig-card {{
  background: #ffffff;
  border: 1px solid #e2e8f0;
  border-radius: 10px;
  padding: 18px 16px 10px;
}}
.fig-label {{
  font-size: 11.5px;
  font-weight: 600;
  text-transform: uppercase;
  letter-spacing: 0.07em;
  color: {C_SLT};
  margin-bottom: 10px;
}}
.empty-fig {{
  padding: 36px 20px;
  text-align: center;
  color: #94a3b8;
  font-size: 13px;
}}

/* ── Tables ── */
.table-wrap {{
  overflow-x: auto;
  border: 1px solid #e2e8f0;
  border-radius: 10px;
  margin-bottom: 32px;
}}
.data-table {{
  width: 100%;
  border-collapse: collapse;
  font-size: 12.5px;
}}
.data-table thead tr {{
  background: #f8fafc;
}}
.data-table th {{
  padding: 11px 14px;
  text-align: left;
  font-weight: 600;
  font-size: 11px;
  text-transform: uppercase;
  letter-spacing: 0.05em;
  color: {C_SLT};
  border-bottom: 1px solid #e2e8f0;
  white-space: nowrap;
}}
.data-table td {{
  padding: 9px 14px;
  border-bottom: 1px solid #f1f5f9;
  color: #334155;
  font-size: 12.5px;
  white-space: nowrap;
}}
.data-table tbody tr:hover {{ background: #f8fafc; }}
.data-table tbody tr:last-child td {{ border-bottom: none; }}

/* ── Section spacing ── */
.trait-section {{ margin-bottom: 16px; }}
section {{ margin-bottom: 48px; }}

/* ── Divider ── */
.divider {{
  height: 1px;
  background: #f1f5f9;
  margin: 44px 0;
}}

/* ── Footer ── */
.page-footer {{
  padding: 28px 0 0;
  border-top: 1px solid #f1f5f9;
  font-size: 12px;
  color: #94a3b8;
  margin-top: 56px;
}}

/* ── Tab nav ── */
.tab-bar {{
  display: flex; gap: 4px; margin-bottom: 28px;
  border-bottom: 1px solid #e2e8f0;
}}
.tab-btn {{
  padding: 9px 18px;
  font-size: 13px;
  font-weight: 500;
  color: {C_SLT};
  background: none;
  border: none;
  border-bottom: 2px solid transparent;
  cursor: pointer;
  transition: color 0.15s, border-color 0.15s;
}}
.tab-btn:hover {{ color: {C_BLUE}; }}
.tab-btn.active {{ color: {C_BLUE}; border-bottom-color: {C_BLUE}; }}
.tab-panel {{ display: none; }}
.tab-panel.active {{ display: block; }}
</style>
</head>
<body>

<!-- ── Sidebar ── -->
<nav class="sidebar">
  <div class="sidebar-logo">
    <span>CDC 1.0.0</span>
    GWAS Meta-Analysis
  </div>
  <div class="nav-section">Overview</div>
  <a href="#overview"  class="nav-link">Summary</a>
  <a href="#global"    class="nav-link">Global Figures</a>
  <a href="#loci"      class="nav-link">Top Loci</a>
  <div class="nav-section" style="margin-top:12px">Per Trait</div>
  {nav_links}
</nav>

<!-- ── Main ── -->
<div class="main">

  <!-- Header -->
  <div class="page-header" id="overview">
    <h1>GWAS Meta-Analysis Report</h1>
    <div class="meta">
      <span>{args.author} · {args.affiliation}</span>
      <span>{args.institute}</span>
      <span style="color:#94a3b8">{timestamp}</span>
    </div>
  </div>

  <!-- Overview cards -->
  <div class="overview-stats">
    <div class="ov-card blue"><div class="num">{n_traits_total}</div>
      <div class="lbl">Traits</div></div>
    <div class="ov-card"><div class="num">2</div>
      <div class="lbl">Cohorts</div></div>
    <div class="ov-card red"><div class="num">{n_sig_total:,}</div>
      <div class="lbl">GW Significant</div></div>
    <div class="ov-card green">
      <div class="num">{sum(len(v) for v in loci_dict.values()):,}</div>
      <div class="lbl">Independent Loci</div></div>
    <div class="ov-card"><div class="num">IVW RE</div>
      <div class="lbl">Method</div></div>
    <div class="ov-card"><div class="num">5×10⁻⁸</div>
      <div class="lbl">Significance</div></div>
  </div>

  <!-- Summary table -->
  <section id="summary-table">
    <h2 class="section-heading">Meta-analysis Summary</h2>
    {summ_tbl}
  </section>

  <!-- Global figures -->
  <section id="global">
    <h2 class="section-heading">Global Figures</h2>
    <div class="fig-grid fig-grid-2col">
      <div class="fig-card">
        <p class="fig-label">Heterogeneity — I² Heatmap</p>
        {het_fig}
      </div>
      <div class="fig-card">
        <p class="fig-label">Genomic Inflation — λ GC</p>
        {lam_fig}
      </div>
    </div>
  </section>

  <!-- Per-trait sections -->
  {trait_sections}

  <div class="divider"></div>

  <!-- Top loci table -->
  <section id="loci">
    <h2 class="section-heading">Top Loci — Annotated</h2>
    {loci_tbl}
  </section>

  <!-- Footer -->
  <div class="page-footer">
    CDC 1.0.0 · GWAS Meta-Analysis Pipeline ·
    {args.author}, {args.institute} ·
    Generated {timestamp} ·
    <a href="https://github.com/nkhan119" style="color:#94a3b8">github.com/nkhan119</a>
  </div>

</div>

<script>
// Active nav on scroll
const sections = document.querySelectorAll('[id]');
const navLinks  = document.querySelectorAll('.nav-link');
const io = new IntersectionObserver(entries => {{
  entries.forEach(e => {{
    if (e.isIntersecting) {{
      navLinks.forEach(l => l.classList.remove('active'));
      const link = document.querySelector(`.nav-link[href="#${{e.target.id}}"]`);
      if (link) link.classList.add('active');
    }}
  }});
}}, {{ rootMargin: '-20% 0px -70% 0px' }});
sections.forEach(s => io.observe(s));
</script>
</body>
</html>'''
    return html


# ── CLI ───────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--meta_dir',   required=True)
    p.add_argument('--het_dir',    required=True)
    p.add_argument('--qq_dir',     required=True)
    p.add_argument('--loci_dir',   required=True)
    p.add_argument('--filter_dir', default='')
    p.add_argument('--author',     default='Nadeem Khan')
    p.add_argument('--affiliation',default='Bioinformatician')
    p.add_argument('--institute',  default='INRS-CAFSB')
    p.add_argument('--p_sig',      type=float, default=5e-8)
    p.add_argument('--p_sug',      type=float, default=1e-5)
    p.add_argument('--i2_thresh',  type=float, default=75)
    p.add_argument('--out',        default='Meta_Report.html')
    return p.parse_args()


def main():
    args = parse_args()
    print('[report] Loading data...')
    meta_dict   = glob_load(os.path.join(args.meta_dir,  '**/*_meta_results.tsv.gz'))
    qq_dict     = glob_load(os.path.join(args.qq_dir,    '**/*_qq_data.tsv.gz'))
    het_dict    = glob_load(os.path.join(args.het_dir,   '**/*_het_stats.tsv.gz'))
    loci_dict   = glob_load(os.path.join(args.loci_dir,  '**/*_annotated_loci.tsv'))
    summ_dict   = glob_load(os.path.join(args.meta_dir,  '**/*_meta_summary.tsv'))
    filter_dict = glob_load(os.path.join(args.filter_dir,'**/*_filter_summary.tsv')) \
                  if args.filter_dir else {}

    print(f'[report] Traits found: {list(meta_dict.keys())}')
    print('[report] Building HTML...')
    html = build_html(args, meta_dict, qq_dict, het_dict,
                      loci_dict, summ_dict, filter_dict)
    with open(args.out, 'w', encoding='utf-8') as fh:
        fh.write(html)
    print(f'[report] Written: {args.out}  ({os.path.getsize(args.out)//1024} KB)')


if __name__ == '__main__':
    main()
