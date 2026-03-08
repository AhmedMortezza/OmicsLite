"""
Step 8 — Build HTML Report
============================
Membaca semua file JSON output dari tahap sebelumnya,
membuat visualisasi Plotly, dan merakit satu HTML report
dengan CSS terpisah (style.css).

Input:
  output_GSE15852/metadata.json
  output_GSE15852/stats.csv
  output_GSE15852/pca_umap.json
  output_GSE15852/heatmap.json
  output_GSE15852/boxplot.json
  output_GSE15852/enrichment.json
  output_GSE15852/network.json
  report/style.css

Output:
  output_GSE15852/GSE15852_Report.html
"""

import json
import os
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio

OUTPUT_DIR = "output_GSE15852"
REPORT_DIR = os.path.dirname(os.path.abspath(__file__))

# ── Load semua data ──────────────────────────────
print("Membaca semua output JSON...")

def load_json(name):
    path = os.path.join(OUTPUT_DIR, name)
    if not os.path.exists(path):
        print(f"  [!] File tidak ditemukan: {path}")
        return None
    with open(path) as f:
        return json.load(f)

meta       = load_json("metadata.json")
pca_umap   = load_json("pca_umap.json")
hm_data    = load_json("heatmap.json")
box_data   = load_json("boxplot.json")
enr_data   = load_json("enrichment.json")
net_data   = load_json("network.json")
stats      = pd.read_csv(os.path.join(OUTPUT_DIR, "stats.csv"))

# Load CSS
css_path = os.path.join(REPORT_DIR, "style.css")
with open(css_path, encoding='utf-8') as f:
    css_content = f.read()

# ── Helper ───────────────────────────────────────
PLOTLY_CONFIG = dict(
    template='plotly_white',
    paper_bgcolor='rgba(0,0,0,0)',
    plot_bgcolor='rgba(0,0,0,0)',
    font=dict(family='IBM Plex Mono, monospace', color='#3d3830', size=11),
    margin=dict(l=60, r=40, t=50, b=60),
)
COLORS = {'Cancer': '#b85c5c', 'Normal': '#7b8c7a', 'Up': '#b85c5c', 'Down': '#4f7fa8', 'NS': '#c4bdb4'}

def to_html(fig):
    return pio.to_html(fig, full_html=False, include_plotlyjs=False,
                       config={'responsive': True, 'displayModeBar': True,
                               'modeBarButtonsToRemove': ['select2d','lasso2d']})

# ═════════════════════════════════════════════════
# VOLCANO PLOT
# ═════════════════════════════════════════════════
print("Membuat Volcano Plot...")
stats_v = stats.copy()
stats_v['neg_log10_padj'] = -np.log10(stats_v['p_adj'].clip(1e-300, 1))
stats_v['label'] = stats_v.apply(
    lambda r: r['gene_symbol'] if (r['DEG'] != 'NS' and r['neg_log10_padj'] > 6
                                   and pd.notna(r['gene_symbol'])) else '', axis=1)

fig_vol = go.Figure()
for deg, color in [('NS', COLORS['NS']), ('Down', COLORS['Down']), ('Up', COLORS['Up'])]:
    sub = stats_v[stats_v['DEG'] == deg]
    fig_vol.add_trace(go.Scatter(
        x=sub['logFC'], y=sub['neg_log10_padj'],
        mode='markers',
        name=deg,
        marker=dict(color=color, size=4, opacity=0.7 if deg=='NS' else 0.9,
                    line=dict(width=0)),
        text=sub['gene_symbol'].fillna(sub['probe_id']),
        customdata=np.stack([sub['logFC'], sub['p_adj'], sub['gene_symbol'].fillna('—')], axis=-1),
        hovertemplate='<b>%{customdata[2]}</b><br>logFC: %{customdata[0]:.3f}<br>adj.p: %{customdata[1]:.2e}<extra></extra>',
    ))

# Label untuk top genes
label_sub = stats_v[stats_v['label'] != '']
for _, r in label_sub.iterrows():
    fig_vol.add_annotation(x=r['logFC'], y=r['neg_log10_padj'],
        text=r['label'], showarrow=False, yshift=8,
        font=dict(size=8, color='#c9d4e8'), bgcolor='rgba(8,11,18,.7)')

fig_vol.add_vline(x=1.5,  line_dash='dash', line_color='#2a3a5c', line_width=1)
fig_vol.add_vline(x=-1.5, line_dash='dash', line_color='#2a3a5c', line_width=1)
fig_vol.add_hline(y=-np.log10(0.05), line_dash='dash', line_color='#2a3a5c', line_width=1)

fig_vol.update_layout(
    **PLOTLY_CONFIG,
    height=520,
    xaxis_title='Log₂ Fold Change (Cancer vs Normal)',
    yaxis_title='-log₁₀(adj. p-value)',
    legend=dict(orientation='h', y=1.02, x=0),
)
volcano_html = to_html(fig_vol)

# ═════════════════════════════════════════════════
# PCA PLOT
# ═════════════════════════════════════════════════
print("Membuat PCA Plot...")
fig_pca = go.Figure()

if pca_umap and 'pca' in pca_umap:
    pca_d  = pca_umap['pca']
    coords = pca_d['coords']
    samps  = pca_d['samples']
    grps   = pca_d['groups']
    ve     = pca_d['var_exp']

    for grp in ['Cancer', 'Normal']:
        idx = [i for i, g in enumerate(grps) if g == grp]
        fig_pca.add_trace(go.Scatter(
            x=[coords[i][0] for i in idx],
            y=[coords[i][1] for i in idx],
            mode='markers+text',
            name=grp,
            text=[samps[i] for i in idx],
            textposition='top center',
            textfont=dict(size=7),
            marker=dict(color=COLORS[grp], size=10,
                        line=dict(width=1, color='rgba(255,255,255,.3)')),
            hovertemplate=f'<b>%{{text}}</b><br>PC1: %{{x:.2f}}<br>PC2: %{{y:.2f}}<extra>{grp}</extra>',
        ))
    fig_pca.update_layout(
        **PLOTLY_CONFIG, height=480,
        xaxis_title=f'PC1 ({ve[0]:.1f}%)',
        yaxis_title=f'PC2 ({ve[1]:.1f}%)',
        legend=dict(orientation='h', y=1.02),
    )
pca_html = to_html(fig_pca)

# ═════════════════════════════════════════════════
# UMAP PLOT
# ═════════════════════════════════════════════════
print("Membuat UMAP Plot...")
fig_umap = go.Figure()
umap_available = pca_umap and pca_umap.get('umap')

if umap_available:
    um   = pca_umap['umap']
    coords_u = um['coords']
    samps_u  = um['samples']
    grps_u   = um['groups']

    for grp in ['Cancer', 'Normal']:
        idx = [i for i, g in enumerate(grps_u) if g == grp]
        fig_umap.add_trace(go.Scatter(
            x=[coords_u[i][0] for i in idx],
            y=[coords_u[i][1] for i in idx],
            mode='markers+text',
            name=grp,
            text=[samps_u[i] for i in idx],
            textposition='top center',
            textfont=dict(size=7),
            marker=dict(color=COLORS[grp], size=10,
                        line=dict(width=1, color='rgba(255,255,255,.3)')),
            hovertemplate=f'<b>%{{text}}</b><br>UMAP1: %{{x:.2f}}<br>UMAP2: %{{y:.2f}}<extra>{grp}</extra>',
        ))
    fig_umap.update_layout(
        **PLOTLY_CONFIG, height=480,
        xaxis_title='UMAP 1', yaxis_title='UMAP 2',
        legend=dict(orientation='h', y=1.02),
    )
    umap_html = to_html(fig_umap)
else:
    umap_html = "<div class='notice'>UMAP tidak tersedia — install <code>umap-learn</code> lalu jalankan ulang.</div>"

# ═════════════════════════════════════════════════
# HEATMAP
# ═════════════════════════════════════════════════
print("Membuat Heatmap...")
if hm_data:
    fig_hm = go.Figure(go.Heatmap(
        z=hm_data['z_matrix'],
        x=hm_data['samples'],
        y=hm_data['genes'],
        colorscale='RdBu_r',
        zmid=0, zmin=-3, zmax=3,
        colorbar=dict(title='Z-score', thickness=12, len=.8),
        hovertemplate='Gene: %{y}<br>Sample: %{x}<br>Z: %{z:.2f}<extra></extra>',
    ))
    # Annotasi group bar
    n_c = hm_data['n_cancer']
    n_n = hm_data['n_normal']
    fig_hm.add_shape(type='rect',
        x0=-.5, x1=n_c-.5, y0=len(hm_data['genes'])-.5, y1=len(hm_data['genes'])+.5,
        fillcolor='#f05d5d', line_width=0)
    fig_hm.add_shape(type='rect',
        x0=n_c-.5, x1=n_c+n_n-.5, y0=len(hm_data['genes'])-.5, y1=len(hm_data['genes'])+.5,
        fillcolor='#4fafed', line_width=0)
    fig_hm.update_layout(**PLOTLY_CONFIG, height=780,
                         xaxis=dict(showticklabels=False, title='Samples'),
                         yaxis=dict(tickfont=dict(size=9)))
    hm_html = to_html(fig_hm)
else:
    hm_html = "<div class='notice'>Data heatmap tidak tersedia.</div>"

# ═════════════════════════════════════════════════
# BOX PLOT
# ═════════════════════════════════════════════════
print("Membuat Box Plot...")
if box_data and box_data.get('genes'):
    fig_box = go.Figure()
    for gd in box_data['genes']:
        for grp in ['Cancer', 'Normal']:
            vals = gd[grp.lower()]
            fig_box.add_trace(go.Box(
                y=vals, name=gd['gene'],
                legendgroup=grp,
                legendgrouptitle_text=grp if gd == box_data['genes'][0] else None,
                showlegend=(gd == box_data['genes'][0]),
                marker_color=COLORS[grp],
                line_color=COLORS[grp],
                fillcolor=COLORS[grp].replace(')', ', 0.15)').replace('rgb', 'rgba') if 'rgb' in COLORS[grp] else COLORS[grp],
                boxpoints='all', jitter=.35, pointpos=0,
                marker=dict(size=4, opacity=.6),
                offsetgroup=grp,
                hovertemplate=f'{grp}: %{{y:.2f}}<extra>{gd["gene"]}</extra>',
            ))
    fig_box.update_layout(
        **PLOTLY_CONFIG, height=520,
        boxmode='group',
        xaxis_tickangle=-30,
        yaxis_title='Normalized Expression',
        legend=dict(orientation='h', y=1.02),
    )
    box_html = to_html(fig_box)
else:
    box_html = "<div class='notice'>Data box plot tidak tersedia.</div>"

# ═════════════════════════════════════════════════
# ENRICHMENT
# ═════════════════════════════════════════════════
print("Membuat Enrichment Chart...")
enr_chart_html  = "<div class='notice'>Enrichment tidak tersedia.</div>"
enr_table_html  = ""

if enr_data and enr_data.get('status') == 'ok' and enr_data.get('pathways'):
    pw = enr_data['pathways']
    df_pw = pd.DataFrame(pw)
    fig_enr = px.bar(
        df_pw.sort_values('combined_score'),
        x='combined_score', y='term',
        color='p_adj',
        color_continuous_scale='Blues_r',
        orientation='h',
        labels={'combined_score': 'Combined Score', 'term': '', 'p_adj': 'adj. p-value'},
        template='plotly_dark',
    )
    fig_enr.update_layout(
        **PLOTLY_CONFIG, height=580,
        yaxis=dict(tickfont=dict(size=9)),
        coloraxis_colorbar=dict(title='adj.p', thickness=12),
    )
    enr_chart_html = to_html(fig_enr)

    # Tabel
    rows_html = ''
    for r in enr_data.get('table_rows', []):
        gs = r.get('Gene_set', '—')
        term = r.get('Term', '—')
        overlap = r.get('Overlap', '—')
        padj = r.get('Adjusted P-value', '—')
        genes = r.get('Genes', '—')
        genes_str = str(genes)[:80] + '…' if len(str(genes)) > 80 else str(genes)
        padj_str = f'{padj:.3e}' if isinstance(padj, float) else str(padj)
        rows_html += f'''<tr>
          <td><span class="geneset-chip">{gs}</span></td>
          <td>{term}</td>
          <td>{overlap}</td>
          <td>{padj_str}</td>
          <td style="color:var(--text-dim);font-size:.7rem">{genes_str}</td>
        </tr>'''
    enr_table_html = f'''
    <table class="enr-table">
      <thead><tr>
        <th>Gene Set</th><th>Term</th><th>Overlap</th><th>adj. p-value</th><th>Genes</th>
      </tr></thead>
      <tbody>{rows_html}</tbody>
    </table>'''

# ═════════════════════════════════════════════════
# NETWORK
# ═════════════════════════════════════════════════
print("Membuat Network Plot...")
net_html = "<div class='notice'>STRING network tidak tersedia (koneksi gagal atau tidak ada edge signifikan).</div>"

if net_data and net_data.get('status') == 'ok':
    nodes = net_data['nodes']
    edges = net_data['edges']

    edge_x, edge_y, edge_w = [], [], []
    node_map = {n['id']: n for n in nodes}
    for e in edges:
        a, b = node_map.get(e['source']), node_map.get(e['target'])
        if a and b:
            edge_x += [a['x'], b['x'], None]
            edge_y += [a['y'], b['y'], None]

    fig_net = go.Figure()
    fig_net.add_trace(go.Scatter(
        x=edge_x, y=edge_y, mode='lines',
        line=dict(color='rgba(42,58,92,.6)', width=1.5),
        hoverinfo='none', showlegend=False,
    ))
    for deg_status in ['Up', 'Down', 'NS']:
        sub = [n for n in nodes if n['deg'] == deg_status]
        if not sub: continue
        fig_net.add_trace(go.Scatter(
            x=[n['x'] for n in sub], y=[n['y'] for n in sub],
            mode='markers+text',
            name=deg_status,
            text=[n['id'] for n in sub],
            textposition='top center',
            textfont=dict(size=8),
            marker=dict(color=COLORS[deg_status], size=13,
                        line=dict(width=1.5, color='rgba(255,255,255,.25)')),
            hovertemplate='<b>%{text}</b><br>logFC: %{customdata:.3f}<extra></extra>',
            customdata=[n['logfc'] for n in sub],
        ))
    fig_net.update_layout(
        **PLOTLY_CONFIG, height=600,
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        legend=dict(orientation='h', y=1.02),
    )
    net_html = to_html(fig_net)

# ═════════════════════════════════════════════════
# DEG TABLE
# ═════════════════════════════════════════════════
deg_show = stats[stats['DEG'] != 'NS'].sort_values('p_adj').head(200)
table_rows = ''
for _, r in deg_show.iterrows():
    fc_class = 'fc-up' if r['logFC'] > 0 else 'fc-down'
    gene     = r['gene_symbol'] if pd.notna(r.get('gene_symbol')) else '—'
    badge    = f'<span class="badge badge-{r["DEG"]}">{r["DEG"]}</span>'
    table_rows += f'''<tr>
      <td class="probe-id">{r["probe_id"]}</td>
      <td class="gene-name">{gene}</td>
      <td class="{fc_class}">{r["logFC"]:+.3f}</td>
      <td>{r["p_value"]:.3e}</td>
      <td>{r["p_adj"]:.3e}</td>
      <td>{badge}</td>
    </tr>'''

# ─── Summary numbers ─────────────────────────────
n_probes   = meta['n_probes']
n_cancer   = meta['n_cancer']
n_normal   = meta['n_normal']
n_up       = (stats['DEG'] == 'Up').sum()
n_down     = (stats['DEG'] == 'Down').sum()
n_sig      = n_up + n_down
gen_ts     = pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')

# ═════════════════════════════════════════════════
# BUILD HTML
# ═════════════════════════════════════════════════
print("Merakit HTML report...")

html = f"""<!DOCTYPE html>
<html lang="id">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1.0"/>
<title>GSE15852 — Analysis Report</title>
<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
<style>
{css_content}
</style>
</head>
<body>

<!-- TOP BAR -->
<div class="topbar">
  <div class="topbar-inner">
    <div class="topbar-label">Transcriptomic Analysis Report</div>
    <h1>GSE<span>15852</span> — Cancer vs Normal</h1>
    <p>Differential expression analysis · FDR correction (Benjamini-Hochberg)<br>
       PCA · UMAP · STRING PPI Network · KEGG/GO Enrichment</p>
  </div>
</div>

<!-- SUMMARY CARDS -->
<div class="cards">
  <div class="card">         <div class="num">{n_probes:,}</div><div class="lbl">Total Probes</div></div>
  <div class="card c-sig">   <div class="num">{n_sig:,}</div>  <div class="lbl">Sig. DEGs (FDR &lt; 0.05)</div></div>
  <div class="card c-up">    <div class="num">{n_up:,}</div>   <div class="lbl">Up-regulated</div></div>
  <div class="card c-down">  <div class="num">{n_down:,}</div> <div class="lbl">Down-regulated</div></div>
  <div class="card">         <div class="num">{n_cancer}</div> <div class="lbl">Cancer Samples</div></div>
  <div class="card">         <div class="num">{n_normal}</div> <div class="lbl">Normal Samples</div></div>
</div>

<!-- TABS -->
<div class="tabs" role="tablist">
  <button class="tab-btn active" onclick="showTab('volcano',this)" role="tab">
    <span class="tab-icon">🌋</span> Volcano
  </button>
  <button class="tab-btn" onclick="showTab('dimred',this)" role="tab">
    <span class="tab-icon">⚬</span> PCA / UMAP
  </button>
  <button class="tab-btn" onclick="showTab('heatmap',this)" role="tab">
    <span class="tab-icon">▦</span> Heatmap
  </button>
  <button class="tab-btn" onclick="showTab('boxplot',this)" role="tab">
    <span class="tab-icon">□</span> Box Plot
  </button>
  <button class="tab-btn" onclick="showTab('network',this)" role="tab">
    <span class="tab-icon">◎</span> Network
  </button>
  <button class="tab-btn" onclick="showTab('enrichment',this)" role="tab">
    <span class="tab-icon">⬡</span> Enrichment
  </button>
  <button class="tab-btn" onclick="showTab('table',this)" role="tab">
    <span class="tab-icon">≡</span> DEG Table
  </button>
</div>

<!-- ══ PANELS ══ -->

<!-- Volcano -->
<div id="tab-volcano" class="panel active" role="tabpanel">
  <div class="plot-card">
    <div class="plot-card-header">
      <h3>Volcano Plot</h3>
      <span class="subtitle">FDR-adjusted p-value · threshold: |logFC| &gt; 1.5, adj.p &lt; 0.05</span>
    </div>
    <div class="legend">
      <div class="legend-item"><div class="legend-dot" style="background:#f05d5d"></div>Up-regulated ({n_up})</div>
      <div class="legend-item"><div class="legend-dot" style="background:#4fafed"></div>Down-regulated ({n_down})</div>
      <div class="legend-item"><div class="legend-dot" style="background:#4a5568"></div>Not Significant</div>
    </div>
    {volcano_html}
  </div>
</div>

<!-- PCA / UMAP -->
<div id="tab-dimred" class="panel" role="tabpanel">
  <div class="two-col">
    <div class="plot-card">
      <div class="plot-card-header">
        <h3>PCA</h3>
        <span class="subtitle">Principal Component Analysis</span>
      </div>
      {pca_html}
    </div>
    <div class="plot-card">
      <div class="plot-card-header">
        <h3>UMAP</h3>
        <span class="subtitle">Uniform Manifold Approximation</span>
      </div>
      {umap_html}
    </div>
  </div>
</div>

<!-- Heatmap -->
<div id="tab-heatmap" class="panel" role="tabpanel">
  <div class="plot-card">
    <div class="plot-card-header">
      <h3>Heatmap — Top 50 DEGs</h3>
      <span class="subtitle">Z-score normalized per gene · Red bar = Cancer · Blue bar = Normal</span>
    </div>
    {hm_html}
  </div>
</div>

<!-- Box Plot -->
<div id="tab-boxplot" class="panel" role="tabpanel">
  <div class="plot-card">
    <div class="plot-card-header">
      <h3>Box Plot — Top 12 DEGs</h3>
      <span class="subtitle">Individual data points overlay</span>
    </div>
    {box_html}
  </div>
</div>

<!-- Network -->
<div id="tab-network" class="panel" role="tabpanel">
  <div class="plot-card">
    <div class="plot-card-header">
      <h3>STRING PPI Network</h3>
      <span class="subtitle">Confidence score ≥ 700 · Red = Up · Blue = Down</span>
    </div>
    {net_html}
  </div>
</div>

<!-- Enrichment -->
<div id="tab-enrichment" class="panel" role="tabpanel">
  <div class="plot-card">
    <div class="plot-card-header">
      <h3>Pathway Enrichment</h3>
      <span class="subtitle">KEGG 2021 Human + GO Biological Process 2023 · via Enrichr</span>
    </div>
    {enr_chart_html}
  </div>
  <div class="plot-card">
    <div class="plot-card-header"><h3>Pathway Table</h3></div>
    <div class="table-wrap">{enr_table_html}</div>
  </div>
</div>

<!-- DEG Table -->
<div id="tab-table" class="panel" role="tabpanel">
  <div class="plot-card">
    <div class="plot-card-header">
      <h3>DEG Results</h3>
      <span class="subtitle">Top 200 · FDR &lt; 0.05, |logFC| &gt; 1.5 · sorted by adj. p-value</span>
    </div>
    <div class="table-wrap">
      <table class="deg-table">
        <thead><tr>
          <th>Probe ID</th>
          <th>Gene Symbol</th>
          <th>logFC</th>
          <th>p-value</th>
          <th>adj. p-value</th>
          <th>Status</th>
        </tr></thead>
        <tbody>{table_rows}</tbody>
      </table>
    </div>
  </div>
</div>

<!-- FOOTER -->
<footer>
  GSE15852 Analysis Pipeline · FDR: Benjamini-Hochberg ·
  <a href="https://string-db.org" target="_blank">STRING</a> PPI ·
  <a href="https://maayanlab.cloud/Enrichr/" target="_blank">Enrichr</a> ·
  Generated {gen_ts}
</footer>

<script>
function showTab(name, btn) {{
  document.querySelectorAll('.panel').forEach(p => p.classList.remove('active'));
  document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
  document.getElementById('tab-' + name).classList.add('active');
  btn.classList.add('active');
}}
</script>
</body>
</html>"""

out_path = os.path.join(OUTPUT_DIR, "GSE15852_Report.html")
with open(out_path, 'w', encoding='utf-8') as f:
    f.write(html)

print(f"  ✅ Report tersimpan: {out_path}")
print(f"  Size: {os.path.getsize(out_path) / 1024:.0f} KB")

# =================================================
# EXPORT PNG
# =================================================
print("\nExport JPG...")

JPG_SCALE  = 2
JPG_WIDTH  = 1400

figures = [
    ("volcano",    fig_vol,  600),
    ("pca",        fig_pca,  560),
    ("heatmap",    fig_hm if hm_data else None,   860),
    ("boxplot",    fig_box if box_data and box_data.get('genes') else None, 600),
    ("network",    fig_net if net_data and net_data.get('status') == 'ok' else None, 700),
    ("enrichment", fig_enr if enr_data and enr_data.get('status') == 'ok' and enr_data.get('pathways') else None, 650),
]

if umap_available:
    figures.append(("umap", fig_umap, 560))

try:
    import kaleido  # noqa
    exported = []
    skipped  = []
    for name, fig, h in figures:
        if fig is None:
            skipped.append(name)
            continue
        # Set background putih untuk export PNG
        fig.update_layout(
            paper_bgcolor='white',
            plot_bgcolor='white',
            font_color='#3d3830',
        )
        jpg_path = os.path.join(OUTPUT_DIR, f"{name}.jpg")
        fig.write_image(jpg_path, format="jpg", width=JPG_WIDTH, height=h, scale=JPG_SCALE)
        # Kembalikan ke transparan supaya tidak merusak HTML
        fig.update_layout(
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
        )
        exported.append(name)
        print(f"  Tersimpan: {jpg_path}")

    if skipped:
        print(f"  Dilewati (data tidak tersedia): {', '.join(skipped)}")
    print(f"\n  Total JPG: {len(exported)} file")

except ImportError:
    print("  [!] kaleido tidak terinstall.")
    print("      Jalankan: pip install kaleido")
except Exception as e:
    print(f"  [!] Export JPG gagal: {e}")