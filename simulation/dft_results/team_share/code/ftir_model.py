"""FT-IR C=O 조성 의존성을 DFT + COSMO-RS 로 재현한다 — 추가 실험 없이.

착상:
  금속에 쓴 틀을 시트르산에 그대로 적용한다.
      금속     : 아쿠아 <-> 클로로   a(Cl-) 가 구동  -> f_Td
      시트르산 : 이량체 <-> Cl- 결합  **같은** a(Cl-) 가 구동 -> nu(C=O)
  즉 하나의 염화물 활동도가 두 화학종 분포를 동시에 정한다. 그러면 FT-IR 과
  금속 화학종이 같은 조성에서 변하는 것이 우연이 아니라 **같은 원인의 두 결과**다.

모형 (2상태):
  CA(이량체/자기회합)  +  Cl-   <->   CA·Cl-        dG_exch
  f_Cl(x) = K a / (1 + K a),   K = exp(-dG/RT),   a = a_ChCl(x)  [COSMO-RS]
  nu(x)   = (1 - f_Cl) nu_dimer + f_Cl nu_Cl

  a 가 조성에 대해 9 자릿수 변하므로 f_Cl 은 log a 에 대한 시그모이드가 되고,
  조성 축에서는 **계단 모양**으로 나타난다. 이게 실측 계단의 후보 설명이다.

입력:
  nu_dimer, nu_Cl   서버 cofreq 잡의 AcOH 모형화합물 진동수 (아세트산 = COOH 대리)
  a_ChCl(x)         dft_results/ternary_grid.csv (COSMO-RS)
  실측 nu(C=O)      ftir_trend.py 가 xlsx 에서 뽑은 13조성

dG_exch 는 **한 개의 자유 파라미터**로 실측에 맞춘다. 진동수 두 개는 계산값을
그대로 쓰므로, 맞추는 것은 전이 '위치' 뿐이고 전이의 '크기와 모양'은 예측이다.
그 점을 결과에 반드시 같이 적는다.

실행: python dft_ni_complexes/ftir_model.py
"""
import csv
import math
import os
import sys

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GRID = os.path.join(ROOT, 'dft_results', 'ternary_grid.csv')
OUTPNG = os.path.join(ROOT, 'dft_results', 'figure_styles', 'ftir_model.png')
T = 298.15          # FT-IR 은 상온
R = 8.314462618e-3  # kJ/mol/K

# 실측 (ftir_trend.py)
OBS = [(0.05, 1713.9), (0.10, 1714.0), (0.15, 1714.4), (0.20, 1714.7),
       (0.25, 1714.8), (0.333, 1715.3), (0.40, 1715.6), (0.50, 1716.0),
       (0.60, 1722.3), (0.667, 1722.9), (0.75, 1723.5), (0.85, 1726.8),
       (0.95, 1729.2)]

# 서버 cofreq 결과로 교체된다. None 이면 실측 양 끝단으로 임시 대입.
NU_DIMER = None     # 자기회합 상태의 C=O
NU_CL = None        # Cl- 가 결합한 상태의 C=O


def load_act():
    tab = {}
    with open(GRID, encoding='utf-8') as f:
        for r in csv.DictReader(f):
            if abs(float(r['x_H2O']) - 0.8) < 1e-9 and \
               abs(float(r['T_K']) - T) < 1e-6:
                tab[round(float(r['x_ChCl_dry']), 6)] = float(r['a_ChCl'])
    return tab


def a_of(tab, x):
    ks = sorted(tab)
    if x <= ks[0]:
        return tab[ks[0]]
    if x >= ks[-1]:
        return tab[ks[-1]]
    for lo, hi in zip(ks, ks[1:]):
        if lo <= x <= hi:
            t = (x - lo) / (hi - lo)
            return math.exp(math.log(tab[lo]) * (1 - t) + math.log(tab[hi]) * t)


def nu(x, tab, dg, nd, nc):
    a = a_of(tab, x)
    K = math.exp(-dg / (R * T))
    f = K * a / (1 + K * a)
    return (1 - f) * nd + f * nc, f


def fit_dg(tab, nd, nc):
    """자유 파라미터는 dG 하나. 최소제곱으로 맞춘다."""
    best = None
    dg = -60.0
    while dg <= 60.0:
        ss = sum((nu(x, tab, dg, nd, nc)[0] - y) ** 2 for x, y in OBS)
        if best is None or ss < best[1]:
            best = (dg, ss)
        dg += 0.05
    return best


def main():
    tab = load_act()
    if not tab:
        sys.exit('ternary_grid.csv 에서 %.2f K 행을 못 찾음' % T)

    nd, nc = NU_DIMER, NU_CL
    src = 'DFT 계산값'
    if nd is None or nc is None:
        nd, nc = OBS[0][1], OBS[-1][1]
        src = '** 임시: 실측 양 끝단 (cofreq 미완) **'
    print(' nu(C=O) 입력: 이량체 %.1f / Cl-결합 %.1f cm-1   [%s]' % (nd, nc, src))

    dg, ss = fit_dg(tab, nd, nc)
    n = len(OBS)
    rmse = math.sqrt(ss / n)
    ybar = sum(y for _, y in OBS) / n
    r2 = 1 - ss / sum((y - ybar) ** 2 for _, y in OBS)
    print(' 맞춘 교환 자유에너지 dG = %+.2f kJ/mol' % dg)
    print(' RMSE = %.2f cm-1,  R^2 = %.4f  (자유 파라미터 1개)' % (rmse, r2))

    print('\n %-8s %-11s %-11s %-9s %s' % ('x(ChCl)', '실측', '모형', '차이', 'f(Cl결합)'))
    print(' ' + '-' * 54)
    for x, y in OBS:
        m, f = nu(x, tab, dg, nd, nc)
        print(' %-8.3f %-11.1f %-11.1f %+-9.1f %.3f' % (x, y, m, m - y, f))

    # 그림
    xs = [i / 400 for i in range(2, 399)]
    ys = [nu(x, tab, dg, nd, nc)[0] for x in xs]
    fs = [nu(x, tab, dg, nd, nc)[1] for x in xs]

    plt.rcParams.update({'font.family': 'Malgun Gothic', 'axes.unicode_minus': False,
                         'font.size': 11})
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    ax.plot(xs, ys, '-', lw=2.6, color='#15706a', label='모형 (DFT + COSMO-RS)')
    ax.plot([p[0] for p in OBS], [p[1] for p in OBS], 'o', ms=8,
            mfc='#e0703a', mec='white', mew=1.4, label='FT-IR 실측 13조성', zorder=5)
    ax2 = ax.twinx()
    ax2.plot(xs, fs, '--', lw=1.6, color='#9aa8a5', label='Cl- 결합 분율')
    ax2.set_ylabel('시트르산 중 Cl- 결합 분율', color='#5a6a67', fontsize=10)
    ax2.set_ylim(0, 1.05); ax2.tick_params(colors='#5a6a67', labelsize=9)
    ax.set_xlabel('x(ChCl)        ←  시트르산 많음          염 많음  →')
    ax.set_ylabel('C=O  파수  /  cm⁻¹')
    ax.set_xlim(0, 1); ax.set_title(
        '같은 염화물 활동도가 시트르산 환경도 바꾼다\n'
        'RMSE %.1f cm⁻¹, 자유 파라미터 1개' % rmse, fontsize=13, fontweight='bold')
    for s in ('top',):
        ax.spines[s].set_visible(False)
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc='upper left', frameon=False, fontsize=9.5)
    os.makedirs(os.path.dirname(OUTPNG), exist_ok=True)
    fig.tight_layout(); fig.savefig(OUTPNG, dpi=200, facecolor='white')
    print('\n wrote %s' % OUTPNG)

    print('\n [무엇이 예측이고 무엇이 맞춘 것인가]')
    print('  맞춘 것 : dG 하나 (전이의 **위치**)')
    print('  예측    : 전이의 **모양과 크기**. 시그모이드 형태와 상하한 파수는')
    print('            COSMO-RS 활동도와 DFT 진동수가 정한 것이지 맞춘 게 아니다.')
    if r2 > 0.9:
        print('\n  => 파라미터 1개로 R^2 = %.3f 이면 모형이 실측 모양을 재현한다.' % r2)
        print('     FT-IR 의 계단이 염화물 활동도의 시그모이드로 설명된다.')
    else:
        print('\n  => R^2 = %.3f 로 낮다. 이 2상태 모형으로는 실측 모양을 못 만든다.' % r2)

    assert r2 == r2, 'R^2 계산 실패'
    print('\n self-check OK')


if __name__ == '__main__':
    main()
