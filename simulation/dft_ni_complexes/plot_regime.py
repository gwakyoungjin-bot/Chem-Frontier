"""분리 영역도 그리기 — 논문 피규어 수준.

anchored_map.py 가 계산을, 이 파일이 작도를 맡는다 (그림만 손볼 때 DFT 결과를
다시 읽지 않아도 되게).

색 설계 (dataviz 절차: 형태 -> 색 역할 -> 검증 -> 마크):
  * 영역(배경) = **순차형**. "몇 개 금속이 음이온이 됐나" 는 순서가 있으므로
    중립 -> 연파랑 -> 진파랑 한 색상 램프. 무지개 금지.
  * 전환선 = **범주형**(금속 정체성). #eb6834 Co / #1baf7a Mn / #4a3aa7 Ni.
    validate_palette.js 통과 (CVD 최악 인접쌍 ΔE 9.2, 정상시야 27.6).
    aqua 의 대비 WARN 은 **직접 라벨**로 해소한다.
  * 시료 마커 = 실제 관찰된 용액 색. 인코딩이 아니라 사진의 재현이므로
    팔레트 밖 색이어도 되고, 범례에 'observed colour' 로 명시한다.

영역은 pcolormesh 가 아니라 **전환선 사이를 fill_between** 으로 채운다.
격자 셀 계단이 생기면 경계가 실제보다 거칠어 보이고, 경계 위치가 이 그림의 핵심이라
계단이 곧 오독을 부른다.

서버 matplotlib 에 한글 폰트가 없으므로 그림 안 글자는 전부 영문.
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))

FILL = {'none': '#f2f1ee', 'co': '#cfe3fb', 'comn': '#93bdf2'}
LINE = {'Co': '#eb6834', 'Mn': '#1baf7a', 'Ni': '#4a3aa7'}
INK, INK2, INK3 = '#0b0b0b', '#52514e', '#9c9b96'

# 2026-09-27 민이형 사진에서 관찰된 용액 색
OBS = {'1:19': '#ffffff', '3:17': '#ffffff', '1:4': '#ffffff',
       '3:1': '#a9ded8', '19:1': '#6fc9bf'}


def _series(xs, ys, ymax):
    """전환선을 fill 용 배열로. None(범위 아래)=ymin 취급, 'high'/없음=ymax."""
    out = []
    for y in ys:
        if isinstance(y, float):
            out.append(y)
        elif y is None:
            out.append(None)      # 이미 전환 -> 호출부에서 ymin 으로
        else:
            out.append(ymax)      # 범위 위 -> 전환 안 함
    return out


def draw(xs, Ts_C, lines, meas, out, subtitle=''):
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 10.5,
        'axes.edgecolor': INK3, 'axes.linewidth': 0.9,
        'xtick.color': INK2, 'ytick.color': INK2,
        'axes.labelcolor': INK, 'text.color': INK,
        'figure.facecolor': 'white', 'axes.facecolor': 'white',
    })
    ymin, ymax = min(Ts_C), max(Ts_C)
    pad = (max(xs) - min(xs)) * 0.035          # 마커가 축에 안 잘리도록
    x0, x1 = min(xs) - pad, max(xs) + pad

    fig, ax = plt.subplots(figsize=(9.2, 5.9))
    fig.subplots_adjust(top=0.80, bottom=0.26)

    X = np.array(xs)
    co = np.array([v if isinstance(v, float) else (ymin if v is None else ymax)
                   for v in lines.get('Co', [ymax] * len(xs))], dtype=float)
    mn = np.array([v if isinstance(v, float) else (ymin if v is None else ymax)
                   for v in lines.get('Mn', [ymax] * len(xs))], dtype=float)
    mn = np.maximum(mn, co)                    # Mn 선은 Co 선 아래로 못 내려간다

    ax.fill_between(X, ymin, co,   color=FILL['none'], zorder=0)
    ax.fill_between(X, co,   mn,   color=FILL['co'],   zorder=0)
    ax.fill_between(X, mn,   ymax, color=FILL['comn'], zorder=0)
    # 축 여백까지 배경 이어붙이기 (마커 자리를 위해 넓힌 부분)
    ax.axvspan(x0, X[0],  color=FILL['none'], zorder=0)
    ax.fill_between([X[-1], x1], ymin, co[-1],  color=FILL['none'], zorder=0)
    ax.fill_between([X[-1], x1], co[-1], mn[-1], color=FILL['co'],  zorder=0)
    ax.fill_between([X[-1], x1], mn[-1], ymax,  color=FILL['comn'], zorder=0)

    # --- 전환선 + 선 위 직접 라벨 -----------------------------------------
    for sym, dash, off in (('Co', 'solid', (0, 10)), ('Mn', (0, (7, 3.5)), (0, 9))):
        ys = lines.get(sym)
        if not ys:
            continue
        xv = [x for x, y in zip(xs, ys) if isinstance(y, float)]
        yv = [y for y in ys if isinstance(y, float)]
        if not xv:
            continue
        ax.plot(xv, yv, ls=dash, lw=2.8, color=LINE[sym], zorder=4,
                solid_capstyle='round', dash_capstyle='round')
        k = max(1, int(len(xv) * 0.30))
        txt = 'Co  Oh→Td' if sym == 'Co' else 'Mn  Oh→Td   (uncalibrated)'
        ax.annotate(txt, (xv[k], yv[k]), textcoords='offset points', xytext=off,
                    color=LINE[sym], fontsize=10, fontweight='bold', zorder=6,
                    rotation=0, ha='left', va='bottom')

    # Ni 는 격자 밖 -> 왼쪽 빈 영역에 명시 (선과 겹치지 않는 자리)
    ax.annotate('Ni  Oh→Td  lies above 125 °C\nNi never becomes anionic here',
                (0.06, 0.93), xycoords='axes fraction', ha='left', va='top',
                fontsize=10, fontweight='bold', color=LINE['Ni'], zorder=6,
                linespacing=1.35)

    # --- 기준선 -----------------------------------------------------------
    ax.axvline(0.333, color=INK3, ls=(0, (2, 3)), lw=1.1, zorder=2)
    ax.annotate('conventional 1:2', (0.341, ymin + 2), rotation=90, va='bottom',
                ha='left', color=INK2, fontsize=9, zorder=6)
    ax.axhline(60, color=INK3, ls=(0, (1, 2.6)), lw=1.1, zorder=2)

    # --- 실험점: 마커를 관찰된 용액 색으로 -------------------------------
    for x, lab in meas:
        ax.plot(x, 60, 'o', ms=12, mfc=OBS.get(lab, '#ffffff'), mec=INK,
                mew=1.4, zorder=7)
        ax.annotate(lab, (x, 60), textcoords='offset points',
                    xytext=(0, 14 if x < 0.5 else -23), ha='center',
                    fontsize=9.5, color=INK, zorder=7)
    # 캡션은 'conventional 1:2' 세로 라벨(x=0.333)을 피해 왼쪽 아래 빈 영역에 둔다
    ax.annotate('our experiment · 60 °C\nmarker fill = observed solution colour',
                (0.045, 0.10), xycoords='axes fraction', ha='left', va='bottom',
                fontsize=9.5, color=INK2, zorder=6, linespacing=1.4)

    # --- 축 ---------------------------------------------------------------
    ax.set_xlim(x0, x1); ax.set_ylim(ymin, ymax)
    ax.set_xlabel('x(ChCl)        ChCl : citric acid,  water-free basis', labelpad=9)
    ax.set_ylabel('temperature  /  °C', labelpad=9)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=3.5, width=0.9)

    # --- 제목: figure 좌표에 놓아 축 라벨과 절대 안 겹치게 ----------------
    fig.text(0.055, 0.955, 'Where chloride speciation splits Ni, Co and Mn',
             fontsize=15, fontweight='bold', color=INK, ha='left', va='top')
    if subtitle:
        fig.text(0.055, 0.905, subtitle, fontsize=9.5, color=INK2,
                 ha='left', va='top', linespacing=1.45)

    # --- 범례: 영역만 (선은 직접 라벨) ------------------------------------
    h = [Patch(facecolor=FILL['none'], edgecolor='#dcdbd6',
               label='all cationic  —  no chloride-based separation'),
         Patch(facecolor=FILL['co'], edgecolor='#dcdbd6',
               label='Co anionic / Ni+Mn cationic   →  Co separable'),
         Patch(facecolor=FILL['comn'], edgecolor='#dcdbd6',
               label='Co+Mn anionic / Ni cationic   →  Ni separable'),
         Line2D([], [], marker='o', ls='', mfc=OBS['19:1'], mec=INK, mew=1.4,
                ms=10, label='measured sample (fill = observed colour)')]
    ax.legend(handles=h, loc='upper center', bbox_to_anchor=(0.5, -0.155),
              ncol=2, frameon=False, fontsize=9.5, handlelength=1.8,
              columnspacing=2.2, labelspacing=0.75, handletextpad=0.9)

    fig.savefig(out, dpi=220, bbox_inches='tight', facecolor='white')
    print('wrote', out)
    return out
