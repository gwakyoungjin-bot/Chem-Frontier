"""포스터 양식에 맞춘 그림 2종 — 흰 배경 + 경희대 크림슨/카키.

왜 다시 만드나:
  기존 D_poster_ko 는 어두운 배경(#0e1413)이었다. 실제 포스터 양식은 흰 배경에
  크림슨(#a40031) 상단바와 카키(#c1ad61) 섹션 라벨을 쓰는 밝은 학술 포스터다.
  어두운 그림을 그대로 올리면 페이지에서 홀로 튄다. 양식 색을 추출해 맞춘다.

  fig1  조성-온도 화학종 지도 + 바이알 10개   -> Results 칸
  fig2  배위 사다리 분자구조 (Ni 계열)         -> Method 칸

fig2 를 Ni 로만 그리는 이유:
  Co/Mn 의 물 많은 구조 5종은 잔여 허수진동수가 있는 안장점이다(CO_REPORT.txt).
  Ni 5종은 전부 통과했으므로 그림에 쓰기 안전하다.

실행: python dft_ni_complexes/plot_poster2.py
"""
import csv
import math
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from PIL import Image

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'dft_results', 'figure_styles')
GEO = os.path.join(ROOT, 'dft_results', 'geometries', 'optimized')

# 포스터 양식에서 추출한 색
CRIMSON, KHAKI = '#a40031', '#c1ad61'
INK, MUTED, LINE = '#1a1a1a', '#5e5e5e', '#d8d3c4'
# 영역 채움: 카키 계열 단색 램프 (순차형). 무지개 금지.
FILL = ['#f6f3ea', '#e2d9bd', '#c1ad61']

# 60 C 사면체 분율 — ftd_60C.csv + ftd_scan (검증됨)
DATA = [('1:19', 0.05, 0.0000, '무색'), ('1:9', 0.10, 0.0000, '무색'),
        ('3:17', 0.15, 0.0000, '무색'), ('1:4', 0.20, 0.0000, '무색'),
        ('1:3', 0.25, 0.0000, '무색'), ('2:3', 0.40, 0.0000, '무색'),
        ('3:2', 0.60, 0.0053, '청록'), ('3:1', 0.75, 0.2424, '청록'),
        ('17:3', 0.85, 0.7668, '청록'), ('19:1', 0.95, 0.9614, '청록')]
# 사진에서 바이알을 잘라낼 위치 (라벨 -> (사진, 중심 x))
VIALS = {'3:17': ('A', 968), '1:4': ('A', 1524), '1:19': ('A', 1767),
         '3:1': ('A', 2374), '19:1': ('A', 2865),
         '17:3': ('B', 1150), '3:2': ('B', 1650), '1:3': ('B', 2180),
         '2:3': ('B', 2620), '1:9': ('B', 3020)}
PHOTO = {'A': os.path.join(ROOT, 'photos', 'vials_260922.png'),
         'B': os.path.join(ROOT, 'photos', 'vials_B_260930.png')}
# 두 사진의 액체 높이가 달라 크롭을 사진별로 잡는다
CROP = {'A': ((1880, 2130), 105), 'B': ((1800, 2060), 115)}


def style():
    plt.rcParams.update({
        'font.family': 'Malgun Gothic', 'font.size': 11,
        'axes.unicode_minus': False, 'figure.facecolor': 'white',
        'axes.facecolor': 'white', 'text.color': INK,
        'axes.labelcolor': INK, 'axes.edgecolor': LINE,
        'xtick.color': MUTED, 'ytick.color': MUTED,
        'svg.fonttype': 'none',        # PPT 에서 글자 편집 가능
    })


def load_window():
    xs, co, mn = [], [], []
    with open(os.path.join(ROOT, 'dft_results', 'window.csv'), encoding='utf-8') as f:
        for r in csv.DictReader(f):
            xs.append(float(r['x_ChCl']))
            for k, d in (('T_Co_C', co), ('T_Mn_C', mn)):
                v = r[k]
                d.append(np.nan if v.startswith('>') else
                         (-99.0 if v.startswith('<') else float(v)))
    return np.array(xs), np.array(co), np.array(mn)


# ---------------------------------------------------------------- fig 1
def fig_map():
    style()
    xs, co, mn = load_window()
    ymin, ymax = 25.0, 125.0
    co = np.where(np.isnan(co), ymax, np.where(co < 0, ymin, co))
    mn = np.maximum(np.where(np.isnan(mn), ymax, np.where(mn < 0, ymin, mn)), co)

    fig = plt.figure(figsize=(10.4, 8.0))
    ax = fig.add_axes([0.095, 0.355, 0.875, 0.455])
    axp = fig.add_axes([0.095, 0.155, 0.875, 0.165], sharex=ax)

    pad = 0.045
    x0, x1 = xs.min() - pad, xs.max() + pad
    ax.fill_between(xs, ymin, co, color=FILL[0], lw=0)
    ax.fill_between(xs, co, mn, color=FILL[1], lw=0)
    ax.fill_between(xs, mn, ymax, color=FILL[2], lw=0)
    ax.axvspan(x0, xs[0], color=FILL[0], lw=0)
    for lo, hi, c in ((ymin, co[-1], FILL[0]), (co[-1], mn[-1], FILL[1]),
                      (mn[-1], ymax, FILL[2])):
        ax.fill_between([xs[-1], x1], lo, hi, color=c, lw=0)

    ok = co < ymax - 0.1
    ax.plot(xs[ok], co[ok], '-', lw=2.4, color=CRIMSON, solid_capstyle='round')
    okm = (mn < ymax - 0.1) & (mn > co + 0.1)
    ax.plot(xs[okm], mn[okm], ls=(0, (6, 3)), lw=2.0, color=CRIMSON, alpha=0.55,
            dash_capstyle='round')

    ax.text(0.085, 100, '세 금속 모두 양이온\n화학종으로 구분되지 않음',
            color=MUTED, fontsize=12, fontweight='bold', ha='left', va='center',
            linespacing=1.5)
    ax.text(0.845, 68, 'Co → 음이온 (사면체)\nNi → 여전히 양이온',
            color=INK, fontsize=10.5, fontweight='bold', ha='center', va='center',
            linespacing=1.4)
    ax.text(0.878, 108, 'Mn 도 합류\n(미보정)', color='#6b5f33', fontsize=9.5,
            fontweight='bold', ha='center', va='center', linespacing=1.3)

    ax.axvline(0.333, color=CRIMSON, lw=2.6, zorder=5)
    ax.annotate('', xy=(0.345, 42), xytext=(0.60, 42), zorder=6,
                arrowprops=dict(arrowstyle='<|-', color=CRIMSON, lw=2.0,
                                mutation_scale=18))
    ax.text(0.615, 42, '관행 몰비 1:2 는\n세 금속이 같은 화학종인\n영역 한복판',
            color=CRIMSON, fontsize=11, fontweight='bold', va='center', zorder=6,
            linespacing=1.35)
    ax.axhline(60, color=MUTED, ls=(0, (1, 2.6)), lw=1.2)
    ax.text(x0 + 0.012, 62.5, '실험 조건 · 60 °C', color=MUTED, fontsize=9.5)

    ax.set_xlim(x0, x1); ax.set_ylim(ymin, ymax)
    ax.set_ylabel('온도  /  °C', labelpad=10, fontsize=11.5)
    ax.tick_params(labelbottom=False, length=3.5)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)

    # 바이알 10개
    src = {k: Image.open(p).convert('RGB') for k, p in PHOTO.items()
           if os.path.exists(p)}
    placed = 0
    for lab, x, f, col in DATA:
        if lab not in VIALS:
            continue
        which, cx = VIALS[lab]
        if which not in src:
            continue
        (cy0, cy1), cw = CROP[which]
        im = src[which].crop((cx - cw, cy0, cx + cw, cy1))
        im = im.resize((62, 74), Image.LANCZOS)
        hit = f > 0
        axp.add_artist(AnnotationBbox(
            OffsetImage(np.asarray(im), zoom=0.92), (x, 0.60), frameon=True,
            pad=0.10, bboxprops=dict(edgecolor=CRIMSON if hit else '#bdbdbd',
                                     lw=2.6 if hit else 0.9)))
        axp.text(x, 0.02, lab, ha='center', va='bottom', fontsize=9.5,
                 fontweight='bold', color=INK)
        placed += 1

    axp.set_ylim(0, 1); axp.set_yticks([])
    axp.set_xlabel('x(ChCl)          ←  시트르산 많음            염(ChCl) 많음  →',
                   labelpad=22, fontsize=11.5)
    for s in ('top', 'right', 'left'):
        axp.spines[s].set_visible(False)
    axp.axvline(0.333, color=CRIMSON, lw=2.6, alpha=0.45, zorder=0)

    fig.text(0.095, 0.965,
             '조성비가 어느 금속을 음이온으로 만들지 결정한다',
             fontsize=17, fontweight='bold', color=CRIMSON, ha='left', va='top')
    fig.text(0.095, 0.918,
             '계산이 사면체 코발트를 예측한 조성과 실제로 착색된 조성이 10개 전부 일치',
             fontsize=11, color=INK, ha='left', va='top')
    fig.text(0.095, 0.885,
             'DFT(r2SCAN-3c) + COSMO-RS · 실험과 동일한 물 함량 · Ni·Co 는 각각 문헌 실측값에 고정',
             fontsize=8.8, color=MUTED, ha='left', va='top')
    fig.text(0.095, 0.030,
             '붉은 테두리 = 계산이 사면체 Co 를 예측한 조성 (분율 0.005 / 0.24 / 0.77 / 0.96)   ·   '
             '나머지 6조성은 계산값이 정확히 0',
             fontsize=9, color=MUTED, ha='left')

    for ext, dpi in (('png', 220), ('svg', None)):
        p = os.path.join(OUT, 'P1_map_10vials.%s' % ext)
        fig.savefig(p, facecolor='white', bbox_inches='tight',
                    **({'dpi': dpi} if dpi else {}))
        print('wrote', p)
    plt.close(fig)
    return placed


# ---------------------------------------------------------------- fig 2
def read_xyz(p):
    L = [l.split() for l in open(p, encoding='utf-8') if l.strip()]
    n = int(L[0][0])
    return [(r[0], tuple(float(v) for v in r[1:4])) for r in L[2:2 + n]]


COLOR = {'Ni': '#4a9d8f', 'O': '#c0392b', 'Cl': '#4a8c3f', 'H': '#f2f2f2'}
RAD = {'Ni': 300, 'O': 150, 'Cl': 230, 'H': 55}


def fig_ladder():
    style()
    names = [('Ni_aq6', r'$[\mathrm{Ni(H_2O)_6}]^{2+}$', '팔면체'),
             ('NiCl1_aq5', r'$[\mathrm{NiCl(H_2O)_5}]^{+}$', '팔면체'),
             ('NiCl2_aq4', r'$[\mathrm{NiCl_2(H_2O)_4}]$', '팔면체'),
             ('NiCl3_aq1', r'$[\mathrm{NiCl_3(H_2O)}]^{-}$', '사면체'),
             ('NiCl4', r'$[\mathrm{NiCl_4}]^{2-}$', '사면체')]
    fig, axes = plt.subplots(1, 5, figsize=(13.2, 3.5))
    for ax, (nm, formula, geom) in zip(axes, names):
        p = os.path.join(GEO, '%s_opt.xyz' % nm)
        at = read_xyz(p)
        m = next(q for e, q in at if e == 'Ni')
        at = [(e, tuple(a - b for a, b in zip(q, m))) for e, q in at]
        # 리간드(O·Cl) 배치의 주성분 2축에 투영 — 그냥 xy 로 자르면 기하가 겹쳐 보인다
        L = np.array([q for e, q in at if e in ('O', 'Cl')])
        if len(L) >= 3:
            u, sv, vt = np.linalg.svd(L - L.mean(0), full_matrices=False)
            at = [(e, tuple(np.array(q) @ vt[:3].T)) for e, q in at]
        # 금속-리간드 결합선
        for e, q in at:
            if e in ('O', 'Cl'):
                ax.plot([0, q[0]], [0, q[1]], '-', color='#b9b9b9', lw=1.6, zorder=1)
        order = sorted(at, key=lambda t: t[1][2])       # z 순으로 겹침 처리
        for e, q in order:
            ax.scatter(q[0], q[1], s=RAD.get(e, 60), c=COLOR.get(e, '#999'),
                       edgecolors='white', linewidths=1.1, zorder=2)
        ax.set_xlim(-3.2, 3.2); ax.set_ylim(-3.2, 3.2)
        ax.set_aspect('equal'); ax.axis('off')
        ax.set_title(formula, fontsize=11.5, fontweight='bold', color=INK, pad=6)
        ax.text(0, -3.0, geom, ha='center', fontsize=10, color=CRIMSON,
                fontweight='bold')

    # 팔면체 / 사면체 경계
    fig.lines.append(plt.Line2D([0.617, 0.617], [0.10, 0.80], color=CRIMSON,
                                lw=1.6, ls=(0, (5, 4)),
                                transform=fig.transFigure))
    fig.text(0.617, 0.055, '배위수 6 → 4', ha='center', fontsize=10,
             color=CRIMSON, fontweight='bold')
    fig.text(0.02, 0.965, '염화물이 물을 밀어내며 배위 구조가 바뀐다',
             fontsize=15, fontweight='bold', color=CRIMSON, va='top')
    fig.text(0.02, 0.905,
             'ORCA 6.1.1 · r2SCAN-3c/CPCM(SMD water) 최적화 구조.  Ni–O 2.08 Å,  Ni–Cl 2.29 Å',
             fontsize=9.2, color=MUTED, va='top')
    lax = fig.add_axes([0.02, 0.005, 0.30, 0.075]); lax.axis('off')
    lax.set_xlim(0, 10); lax.set_ylim(0, 1)
    for i, (e, t) in enumerate((('Ni', 'Ni'), ('O', 'O (물)'),
                                ('Cl', 'Cl'), ('H', 'H'))):
        lax.scatter(i * 2.4 + 0.2, 0.5, s=110, c=COLOR[e],
                    edgecolors='#999999', linewidths=0.9)
        lax.text(i * 2.4 + 0.65, 0.5, t, va='center', fontsize=9.5, color=MUTED)
    fig.text(0.34, 0.035, 'Co·Mn 도 같은 사다리로 계산했다 (총 15종)',
             fontsize=9, color=MUTED)
    fig.subplots_adjust(top=0.80, bottom=0.14, left=0.02, right=0.99, wspace=0.05)

    for ext, dpi in (('png', 220), ('svg', None)):
        p = os.path.join(OUT, 'P2_ladder_Ni.%s' % ext)
        fig.savefig(p, facecolor='white', bbox_inches='tight',
                    **({'dpi': dpi} if dpi else {}))
        print('wrote', p)
    plt.close(fig)


def main():
    os.makedirs(OUT, exist_ok=True)
    n = fig_map()
    print('  바이알 %d개 배치' % n)
    fig_ladder()
    assert n >= 8, '바이알이 %d개뿐 — 사진 좌표 확인 필요' % n
    print('self-check OK')


if __name__ == '__main__':
    main()
