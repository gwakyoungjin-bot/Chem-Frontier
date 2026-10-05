"""분리 영역도 스타일 후보 — 비교용으로 여러 벌 뽑는다.

배경:
  처음 그림(guideline_map.png)은 색 = Ni 의 f_Td 라는 **연속 sequential** 값이라
  viridis 가 정확히 맞는 용도였다. 지금 영역도는 색이 **이산 범주**(어느 조합이
  음이온인가)라 역할이 다르다. 그래서 viridis 를 그대로 옮기면 용도가 어긋난다.

  두 갈래로 나눠 각각 제대로 만든다:
    A  viridis small multiples — 금속별 f_Td 연속장 3패널. viridis 본래 용도.
    B  영역도 + viridis 에서 뽑은 이산 색 — 영역도 구조에 그 색감만 가져옴.
    C  영역도 + 검증된 범주형 팔레트 (plot_regime.py, 현재본)

  viridis 는 흔히 말하는 '무지개(jet)' 가 아니다. 명도가 단조 증가하도록 설계된
  perceptually-uniform 맵이라 sequential 용으로 권장되는 쪽이다. 다만 **이산 3단계**
  에 쓰면 단일 색상 램프가 아니게 되므로, B 는 그 점을 감수한 선택이다.

실행: anchored_map.py 가 호출한다 (직접 실행용 아님).
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

INK, INK2, INK3 = '#0b0b0b', '#52514e', '#9c9b96'
OBS = {'1:19': '#ffffff', '3:17': '#ffffff', '1:4': '#ffffff',
       '3:1': '#a9ded8', '19:1': '#6fc9bf'}
MEAS_X = {'1:19': 0.05, '3:17': 0.15, '1:4': 0.20, '3:1': 0.75, '19:1': 0.95}


def _style():
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 10.5,
        'axes.edgecolor': INK3, 'axes.linewidth': 0.9,
        'xtick.color': INK2, 'ytick.color': INK2,
        'axes.labelcolor': INK, 'text.color': INK,
        'figure.facecolor': 'white', 'axes.facecolor': 'white',
    })


# ---------------------------------------------------------------- A
def small_multiples(xs, Tc, fields, out):
    """금속별 f_Td 연속장 3패널. viridis 를 본래 용도(연속 sequential)로 쓴다."""
    _style()
    syms = [s for s in ('Ni', 'Co', 'Mn') if s in fields]
    fig, axes = plt.subplots(1, len(syms), figsize=(4.0 * len(syms) + 1.2, 4.5),
                             sharey=True)
    if len(syms) == 1:
        axes = [axes]
    im = None
    for ax, s in zip(axes, syms):
        im = ax.pcolormesh(xs, Tc, np.array(fields[s]), cmap='viridis',
                           vmin=0, vmax=1, shading='gouraud', rasterized=True)
        cs = ax.contour(xs, Tc, np.array(fields[s]), levels=[0.5],
                        colors='white', linewidths=2.0)
        ax.clabel(cs, fmt={0.5: '50%'}, fontsize=8, inline=True)
        for lab, x in MEAS_X.items():
            ax.plot(x, 60, 'o', ms=8, mfc=OBS[lab], mec='w', mew=1.4, zorder=6)
        ax.axhline(60, color='w', ls=(0, (1, 2.5)), lw=1.0, alpha=0.8)
        ax.axvline(0.333, color='w', ls=(0, (2, 3)), lw=1.0, alpha=0.8)
        ax.set_title(s, fontsize=13, fontweight='bold', color=INK, pad=8)
        ax.set_xlabel('x(ChCl)')
        pd = (max(xs) - min(xs)) * 0.04          # 마커가 축에 안 잘리게
        ax.set_xlim(min(xs) - pd, max(xs) + pd); ax.set_ylim(min(Tc), max(Tc))
        ax.tick_params(length=3.5, width=0.9)
    axes[0].set_ylabel('temperature  /  °C')
    cb = fig.colorbar(im, ax=list(axes), fraction=0.022, pad=0.035, aspect=28)
    cb.set_label('tetrahedral (anionic) fraction', fontsize=10)
    cb.outline.set_visible(False)
    fig.text(0.055, 0.985, 'Tetrahedral fraction of each metal', fontsize=15,
             fontweight='bold', ha='left', va='top')
    fig.text(0.055, 0.935,
             'white line = 50 % contour · circles = our samples at 60 °C '
             '(fill = observed solution colour)',
             fontsize=9.5, color=INK2, ha='left', va='top')
    fig.subplots_adjust(top=0.78)
    fig.savefig(out, dpi=220, bbox_inches='tight', facecolor='white')
    print('wrote', out)


# ---------------------------------------------------------------- B
def regime_viridis(xs, Tc, lines, out, subtitle=''):
    """영역도 구조 + viridis 에서 뽑은 이산 색. C 와 구조는 같고 색만 다르다."""
    _style()
    ymin, ymax = min(Tc), max(Tc)
    pad = (max(xs) - min(xs)) * 0.035
    x0, x1 = min(xs) - pad, max(xs) + pad
    vir = plt.get_cmap('viridis')
    fill = [vir(0.06), vir(0.52), vir(0.92)]      # 보라 -> 청록 -> 노랑
    lcol = ['#ffffff', '#ffffff', '#ffffff']      # 진한 바탕이라 선은 흰색

    fig, ax = plt.subplots(figsize=(9.2, 5.9))
    fig.subplots_adjust(top=0.80, bottom=0.26)
    X = np.array(xs)
    co = np.array([v if isinstance(v, float) else (ymin if v is None else ymax)
                   for v in lines.get('Co', [ymax] * len(xs))], float)
    mn = np.maximum(np.array([v if isinstance(v, float) else (ymin if v is None else ymax)
                              for v in lines.get('Mn', [ymax] * len(xs))], float), co)
    for lo, hi, c in ((np.full_like(co, ymin), co, fill[0]), (co, mn, fill[1]),
                      (mn, np.full_like(mn, ymax), fill[2])):
        ax.fill_between(X, lo, hi, color=c, zorder=0)
    ax.axvspan(x0, X[0], color=fill[0], zorder=0)
    for lo, hi, c in ((ymin, co[-1], fill[0]), (co[-1], mn[-1], fill[1]),
                      (mn[-1], ymax, fill[2])):
        ax.fill_between([X[-1], x1], lo, hi, color=c, zorder=0)

    for sym, dash in (('Co', 'solid'), ('Mn', (0, (7, 3.5)))):
        ys = lines.get(sym)
        if not ys:
            continue
        xv = [x for x, y in zip(xs, ys) if isinstance(y, float)]
        yv = [y for y in ys if isinstance(y, float)]
        if not xv:
            continue
        ax.plot(xv, yv, ls=dash, lw=2.8, color='w', zorder=4,
                solid_capstyle='round', dash_capstyle='round')
        # 두 라벨이 서로 겹치지 않게: Co 는 선 왼쪽 위, Mn 은 선 오른쪽 아래
        k = max(1, int(len(xv) * (0.12 if sym == 'Co' else 0.55)))
        off = (6, 11) if sym == 'Co' else (10, -34)
        ax.annotate('%s  Oh→Td%s' % (sym, '\n(uncalibrated)' if sym == 'Mn' else ''),
                    (xv[k], yv[k]), textcoords='offset points', xytext=off,
                    color='w', fontsize=10, fontweight='bold', zorder=6,
                    linespacing=1.3)
    ax.annotate('Ni  Oh→Td  lies above 125 °C\nNi never becomes anionic here',
                (0.06, 0.93), xycoords='axes fraction', ha='left', va='top',
                fontsize=10, fontweight='bold', color='w', zorder=6, linespacing=1.35)
    ax.axvline(0.333, color='w', ls=(0, (2, 3)), lw=1.1, alpha=0.75, zorder=2)
    ax.annotate('conventional 1:2', (0.341, ymin + 2), rotation=90, va='bottom',
                ha='left', color='w', fontsize=9, zorder=6)
    ax.axhline(60, color='w', ls=(0, (1, 2.6)), lw=1.1, alpha=0.75, zorder=2)
    # 배경이 진해 마커 채움색이 묻힌다. 흰 링을 두껍게 둘러 어떤 배경에서도
    # '관찰된 색' 이 읽히게 한다 — 이 마커가 예측-실측 대조의 핵심이므로.
    for lab, x in MEAS_X.items():
        ax.plot(x, 60, 'o', ms=15.5, mfc='none', mec='w', mew=3.4, zorder=6)
        ax.plot(x, 60, 'o', ms=12, mfc=OBS[lab], mec=INK, mew=1.0, zorder=7)
        ax.annotate(lab, (x, 60), textcoords='offset points',
                    xytext=(0, 16 if x < 0.5 else -25), ha='center',
                    fontsize=9.5, color='w', zorder=7)
    ax.annotate('our experiment · 60 °C\nmarker fill = observed colour',
                (0.032, 0.05), xycoords='axes fraction', ha='left', va='bottom',
                fontsize=9.5, color='w', zorder=6, linespacing=1.4)
    ax.set_xlim(x0, x1); ax.set_ylim(ymin, ymax)
    ax.set_xlabel('x(ChCl)     —     ChCl : citric acid  (water-free basis)', labelpad=9)
    ax.set_ylabel('temperature  /  °C', labelpad=9)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    fig.text(0.055, 0.955, 'Where chloride speciation splits Ni, Co and Mn',
             fontsize=15, fontweight='bold', ha='left', va='top')
    if subtitle:
        fig.text(0.055, 0.905, subtitle, fontsize=9.5, color=INK2, ha='left',
                 va='top', linespacing=1.45)
    h = [Patch(facecolor=fill[0], label='all cationic  —  no separation'),
         Patch(facecolor=fill[1], label='Co anionic / Ni+Mn cationic  →  Co separable'),
         Patch(facecolor=fill[2], label='Co+Mn anionic / Ni cationic  →  Ni separable'),
         Line2D([], [], marker='o', ls='', mfc=OBS['19:1'], mec=INK, mew=1.2,
                ms=10, label='measured sample (fill = observed colour)')]
    ax.legend(handles=h, loc='upper center', bbox_to_anchor=(0.5, -0.155), ncol=2,
              frameon=False, fontsize=9.5, handlelength=1.8, columnspacing=2.2,
              labelspacing=0.75, handletextpad=0.9)
    fig.savefig(out, dpi=220, bbox_inches='tight', facecolor='white')
    print('wrote', out)
