#!/usr/bin/env python3.7
"""
Visualize statistics of all supersingular cases across g=1..12.
Parses pre-computed compact log files. Produces one PNG figure.
Run with: python3.7 visualize.py
"""

import re
import ast
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# ─────────────────────────────────────────────────────────────
# Parse log files
# ─────────────────────────────────────────────────────────────

LOG_FILES = [
    ("output_g123456.log", range(1, 7)),
    *((f"output_g{g}.log", [g]) for g in range(7, 13)),
]

SEQ_RE = re.compile(
    r'^\s+\d+\.\s+a-seq=\[.*?\]\s+a-number=(\d+)\s+type-seq=(\[.*?\])\s+sum=(\d+)\s+Isog=(\[.*?\])'
)
G_RE = re.compile(r'^g = (\d+)')


def parse_logs():
    data = {}
    for filename, _ in LOG_FILES:
        if not os.path.exists(filename):
            print(f"Warning: {filename} not found, skipping")
            continue
        current_g = None
        with open(filename) as f:
            for line in f:
                m = G_RE.match(line)
                if m:
                    current_g = int(m.group(1))
                    if current_g not in data:
                        data[current_g] = []
                    continue
                m = SEQ_RE.match(line)
                if m and current_g is not None:
                    an        = int(m.group(1))
                    type_seq  = ast.literal_eval(m.group(2))
                    sum_val   = int(m.group(3))
                    isog      = ast.literal_eval(m.group(4))
                    data[current_g].append((sum_val, type_seq, isog, an))
    return data


# ─────────────────────────────────────────────────────────────
# Figure: Sum frequency bar charts per g, color-coded by a-number
# ─────────────────────────────────────────────────────────────

# Palette — 12 distinct, warm/earthy tones
AN_COLORS = {
    1:  '#7B8B3A',  # olive green
    2:  '#C03530',  # brick red
    3:  '#E07828',  # burnt orange
    4:  '#5080B0',  # slate blue
    5:  '#C89A18',  # amber gold
    6:  '#7868A8',  # periwinkle
    7:  '#3D8C68',  # teal
    8:  '#D45838',  # coral
    9:  '#A87898',  # dusty mauve
    10: '#688028',  # yellow-olive
    11: '#804020',  # dark brown
    12: '#60A0C0',  # steel blue
}


def fig_sum_frequency(data):
    from collections import Counter, defaultdict
    from matplotlib.ticker import MaxNLocator
    import matplotlib.patches as mpatches
    import matplotlib.lines as mlines

    gs = [g for g in sorted(data) if g >= 2]
    ncols = 4
    nrows = (len(gs) + ncols - 1) // ncols   # 3 rows: g=2..12 fills 11 of 12 slots

    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(ncols * 4.2, nrows * 3.2 + 0.6),
                             gridspec_kw={'hspace': 0.56, 'wspace': 0.35})
    axes_flat = axes.flatten()

    for i, g in enumerate(gs):
        ax = axes_flat[i]
        lo = g * (g - 1) // 2
        hi = g * g * (g - 1) // 2

        by_an = defaultdict(Counter)
        for sum_val, _ts, _isog, an in data[g]:
            by_an[an][sum_val] += 1

        all_xs = sorted(set(row[0] for row in data[g]))
        bottoms = np.zeros(len(all_xs))
        xs_arr = np.array(all_xs, dtype=float)

        for an in sorted(by_an.keys()):
            heights = np.array([by_an[an].get(x, 0) for x in all_xs], dtype=float)
            ax.bar(xs_arr, heights, bottom=bottoms,
                   color=AN_COLORS[an], width=1.0, zorder=2, linewidth=0)
            bottoms += heights

        ax.axvline(lo, color='#333333', linestyle=':', linewidth=1.3, zorder=4)
        ax.axvline(hi, color='#333333', linestyle='--', linewidth=1.3, zorder=4)

        ax.set_title(f'g = {g}', fontsize=11, fontweight='bold', pad=4)
        ax.set_xlabel(r'$\sigma$', fontsize=10)
        ax.set_ylabel('count', fontsize=8)
        ax.tick_params(labelsize=7)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=5, prune='both'))
        ax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
        ax.set_xlim(lo - (hi - lo) * 0.04, hi + (hi - lo) * 0.04)
        ax.grid(axis='y', alpha=0.25, zorder=0, linestyle=':')
        ax.set_axisbelow(True)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # ── Legend box in the 12th slot (where g=13 would go) ───────
    legend_ax = axes_flat[len(gs)]   # slot 11 (0-indexed), bottom-right
    legend_ax.set_visible(True)
    legend_ax.set_xlim(0, 1)
    legend_ax.set_ylim(0, 1)
    legend_ax.set_xticks([])
    legend_ax.set_yticks([])
    legend_ax.set_facecolor('#F7F5F2')
    for spine in legend_ax.spines.values():
        spine.set_linewidth(0.8)
        spine.set_color('#AAAAAA')

    max_an = max(max(row[3] for row in data[g]) for g in gs)
    patches = [mpatches.Patch(facecolor=AN_COLORS[k], edgecolor='none',
                               label=f'a-number = {k}')
               for k in range(1, max_an + 1)]
    line_ss = mlines.Line2D([], [], color='#333333', linestyle=':',
                             linewidth=1.4, label='superspecial')
    line_sg = mlines.Line2D([], [], color='#333333', linestyle='--',
                             linewidth=1.4, label='supergeneral')

    leg = legend_ax.legend(handles=patches + [line_ss, line_sg],
                           loc='center', ncol=2,
                           fontsize=8, frameon=False,
                           handlelength=1.3, handleheight=1.0,
                           columnspacing=0.8, labelspacing=0.45,
                           borderpad=0.5, title='Color code')
    leg.get_title().set_fontsize(9)
    leg.get_title().set_fontweight('bold')

    # Hide any remaining empty slots
    for j in range(len(gs) + 1, len(axes_flat)):
        axes_flat[j].set_visible(False)

    fig.suptitle(r'Statistics of generalized Artin invariant $\sigma$',
                 fontsize=14, fontweight='bold', y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig('sum_frequency.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("Saved sum_frequency.png")


# ─────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print("Parsing log files...")
    data = parse_logs()
    print(f"Loaded g values: {sorted(data.keys())}")
    for g in sorted(data):
        print(f"  g={g}: {len(data[g])} sequences")

    print("\nGenerating figures...")
    fig_sum_frequency(data)
    print("\nDone.")
