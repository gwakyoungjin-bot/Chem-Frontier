"""Ni · Co · Mn 화학종을 한 지도에 올려 '분리 영역도'를 낸다.

왜 여러 금속인가:
  Nature Commun 2021 (10.1038/s41467-021-26814-7) — 염화물 농도로 Co 는 음이온성
  CoCl4(2-) 가 되고 Ni 는 양이온으로 남는다. 이 갈라짐이 Co/Ni 분리의 근거다.
  금속마다 넘어가는 조건이 다르므로, 여러 금속을 겹치면 '어느 조합을 가를 수 있는가'
  가 영역으로 나타난다. 한 금속만으로는 단방향 "넘지 마라" 밖에 안 나온다.

  NCM622 는 Ni:Co:Mn = 6:2:2 다. Mn 을 빼면 금속의 1/3 을 무시하는 것이므로 넣는다.

Mn 에 대한 주의:
  Mn(II) 는 d5 고스핀이라 d-d 전이가 전부 스핀금지다 -> 분광학적으로는 거의 무색.
  따라서 Mn 은 dG(열역학)만 계산하고 dd(분광) 는 하지 않는다.
  ** 그리고 Ni < Mn < Co 같은 순서를 문헌으로 확정하지 못했다 (검색 실패).
     즉 이 계산은 '알려진 순서를 재현하는지' 를 보는 **검증 기회**이지,
     이미 아는 답을 확인하는 게 아니다. 결과 해석 시 이 점을 분명히 할 것. **

보정값 전용:
  Ni 에서 앵커 2점으로 구한 (shift, alpha) 를 Co·Mn 에 그대로 적용한다.
  shift 는 Cl- 한 개당 용매화 오차, alpha 는 이온쌍 대리변수 -> 자유 Cl- 환산 지수.
  둘 다 **용매 모델의 성질**이지 금속의 성질이 아니고, 세 사다리 모두 전하가
  +2 -> -2 로 똑같이 변한다. 따라서 금속 간 차이는 순수하게 DFT dG 차이에서 나온다.
  (가정이므로 한계 절에 반드시 기재)

실행:
    python cobalt_window.py --water     # 권장 (실험 조건의 물 포함)
    python cobalt_window.py             # 2성분 활동도
"""
import csv, math, os, sys

import speciation as S

HERE = os.path.dirname(os.path.abspath(__file__))


def ladder(sym):
    """(이름, Cl 개수, 물 개수, 기하) — 금속 기호만 바꾼 같은 사다리."""
    return [('%s_aq6' % sym,    0, 6, 'Oh'),
            ('%sCl1_aq5' % sym, 1, 5, 'Oh'),
            ('%sCl2_aq4' % sym, 2, 4, 'Oh'),
            ('%sCl3_aq1' % sym, 3, 1, 'Td'),
            ('%sCl4' % sym,     4, 0, 'Td')]


# (기호, 결과폴더, 색)
METALS = [('Ni', HERE, 'tab:blue'),
          ('Co', os.path.join(HERE, 'cobalt'), 'tab:red'),
          ('Mn', os.path.join(HERE, 'manganese'), 'tab:green')]

AMU = {'Cl': 34.969}


def harvest_metal(sym, d):
    """그 금속 폴더의 freq -> {이름: (dH, dS)}. 참조종(H2O, Cl-)은 Ni 폴더 것 재사용."""
    lad = ladder(sym)
    imag = []
    for name, n, m, geom in lad:
        c, mx = S.count_imaginary(os.path.join(d, '%s_freq.out' % name))
        if c:
            imag.append('%s: 허수 %d개 (최대 %.0f cm-1)' % (name, c, mx))

    def dG_at(T):
        G = {}
        for name, n, m, geom in lad:
            g = S.read_gibbs_eff(os.path.join(d, '%s_freq.out' % name), T)
            if g is not None:
                G[name] = g * S.HARTREE
        gw = S.read_gibbs_eff(os.path.join(HERE, 'H2O_freq.out'), T)
        ecl = S.read_sp_energy(os.path.join(HERE, 'Cl_ion_sp.out'))
        ref = '%s_aq6' % sym
        if ref not in G or gw is None or ecl is None:
            return None
        gw *= S.HARTREE
        gcl = S.atom_gibbs(ecl, AMU['Cl'], T)
        return {name: G[name] + (6 - m) * gw - G[ref] - n * gcl
                for name, n, m, geom in lad if name in G}

    lo, hi = dG_at(S.T_LO), dG_at(S.T_HI)
    if lo is None or hi is None:
        return None, imag, [n for n, _, _, _ in lad]
    th = {k: ((lo[k] + S.T_LO * (-(hi[k] - lo[k]) / (S.T_HI - S.T_LO))),
              -(hi[k] - lo[k]) / (S.T_HI - S.T_LO))
          for k in lo if k in hi}
    return th, imag, [n for n, _, _, _ in lad if n not in th]


def f_td(th, lad, aCl, aW, T, shift, alpha):
    num, tot = {}, 0.0
    ln_a = alpha * math.log(max(aCl, 1e-300))
    ln_w = math.log(max(aW, 1e-300))
    for name, n, m, geom in lad:
        if name not in th:
            continue
        dH, dS = th[name]
        v = math.exp(max(-700.0, min(700.0,
            -((dH - T * dS) + n * shift) / (S.R_KJ * T) + n * ln_a - (6 - m) * ln_w)))
        num[name] = v
        tot += v
    if tot <= 0:
        return 0.0
    td = {n for n, _, _, g in lad if g == 'Td'}
    return sum(v for k, v in num.items() if k in td) / tot


def main():
    use_water = '--water' in sys.argv
    aW = 1.0

    # --- Ni 로 보정값을 구한다 (앵커가 Ni 에만 있다) -----------------------
    th_ni, miss_ni, imag_ni = S.harvest()
    if th_ni is None:
        sys.exit('Ni 결과 부족: %s' % ', '.join(miss_ni))
    grid = S.load_activity(use_water)
    extra, _ = S.harvest_extra()
    shift, alpha = S.calibrate(th_ni, grid, aW, S.FTD_MAX, extra)
    if shift is None:
        sys.exit('Ni 앵커 보정 실패 — 전용할 보정값이 없습니다.')
    if alpha is None:
        alpha = 1.0
        print('⚠ alpha 결정 실패 -> 1.0. 지도를 정량으로 읽지 마세요.')

    print('활동도 모델: %s'
          % ('3성분 ChCl:CA:H2O (x_H2O=0.80, 실험조건)' if use_water else '2성분 (물 없음)'))
    print('Ni 보정값을 전 금속에 전용: shift=%+.1f kJ/mol per Cl-, alpha=%.4f'
          % (shift, alpha))

    # --- 금속별 열역학 수확 -------------------------------------------------
    data = []
    for sym, d, color in METALS:
        if sym == 'Ni':
            th, imag, missing = th_ni, imag_ni, miss_ni
        else:
            th, imag, missing = harvest_metal(sym, d)
        if th is None:
            print('  %s: 결과 없음 -> 건너뜀 (%s)' % (sym, d))
            continue
        for wmsg in imag:
            print('  ⚠ %s' % wmsg)
        if missing:
            print('  %s: 빠진 종 %s (그 종 없이 계산)' % (sym, ', '.join(missing)))
        data.append((sym, ladder(sym), th, color))
    if not data:
        sys.exit('쓸 수 있는 금속이 없습니다.')
    print('사용 금속: %s' % ', '.join(s for s, _, _, _ in data))

    print('\n반응 열역학 (각 금속의 아쿠아 착물 기준, 보정 전)')
    print('  %-12s %-4s %11s %11s %9s' % ('종', '기하', 'dH kJ/mol', 'dS J/mol/K', 'dG(60C)'))
    for sym, lad, th, _ in data:
        for name, n, m, geom in lad:
            if name in th:
                dH, dS = th[name]
                print('  %-12s %-4s %+11.1f %+11.1f %+9.1f'
                      % (name, geom, dH, dS * 1000, dH - 333.15 * dS))

    # --- 전환선 ------------------------------------------------------------
    xs = [round(0.05 * i, 2) for i in range(1, 20)]
    Ts = [298.15 + 2.5 * i for i in range(0, 41)]
    # 전환선. 교차가 없을 때 '범위 아래' 인지 '범위 위' 인지 구분해야 한다.
    #  None  = 격자 최저온도에서 이미 사면체 (선이 범위 아래에 있음)
    #  'high'= 격자 최고온도까지도 팔면체 (선이 범위 위에 있음)
    lines = {}
    for sym, lad, th, _ in data:
        ys = []
        for x in xs:
            f0 = f_td(th, lad, S.a_cl(grid, x, Ts[0]), aW, Ts[0], shift, alpha)
            if f0 >= 0.5:
                ys.append(None)            # 이미 넘어가 있음
                continue
            prev, t = None, 'high'
            for T in Ts:
                f = f_td(th, lad, S.a_cl(grid, x, T), aW, T, shift, alpha)
                if prev and prev[1] < 0.5 <= f:
                    T0, f0b = prev
                    t = T0 + (T - T0) * (0.5 - f0b) / (f - f0b) - 273.15
                    break
                prev = (T, f)
            ys.append(t)
        lines[sym] = ys

    syms = [s for s, _, _, _ in data]
    print('\n전환선 (f_Td = 0.5 가 되는 온도, °C)')
    print('  %-8s %s' % ('x(ChCl)', ' '.join('%-10s' % s for s in syms)))
    rows = []
    for i, x in enumerate(xs):
        cells = []
        r = dict(x_ChCl=x)
        for s in syms:
            v = lines[s][i]
            cells.append('%-10s' % ('<25' if v is None else ('>125' if v == 'high' else '%.0f' % v)))
            r['T_%s_C' % s] = ('below_25' if v is None else ('above_125' if v == 'high' else round(v, 1)))
        print('  %-8.2f %s' % (x, ' '.join(cells)))
        rows.append(r)
    out = os.path.join(HERE, 'window.csv')
    with open(out, 'w', newline='', encoding='utf-8') as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0]))
        wr.writeheader(); wr.writerows(rows)
    print('\nwrote', out)

    # --- 분리 영역도 -------------------------------------------------------
    # 색 = '음이온이 된 금속의 조합'. 금속이 늘어도 축을 더 쓰지 않는다.
    try:
        import numpy as np, matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.colors import ListedColormap, BoundaryNorm
        from matplotlib.patches import Patch
    except ImportError:
        print('matplotlib 없음 -> 그림 생략')
        return

    Tgrid = [298.15 + 2.5 * i for i in range(0, 41)]
    cats, Z = {}, []
    for T in Tgrid:
        row = []
        for x in xs:
            a = S.a_cl(grid, x, T)
            anionic = tuple(s for s, lad, th, _ in data
                            if f_td(th, lad, a, aW, T, shift, alpha) >= 0.5)
            if anionic not in cats:
                cats[anionic] = len(cats)
            row.append(cats[anionic])
        Z.append(row)
    Z = np.array(Z)

    def label(t):
        # 그림 안 글자는 영문으로. 서버 matplotlib 에 한글 폰트가 없어 네모로 깨진다.
        if not t:
            return 'all cationic\n(no separation)'
        rest = [s for s in syms if s not in t]
        if not rest:
            return 'all anionic\n(no separation)'
        return '%s anionic / %s cationic\n-> separable' % ('+'.join(t), '+'.join(rest))

    order = sorted(cats, key=lambda t: len(t))
    remap = {cats[t]: i for i, t in enumerate(order)}
    Z = np.vectorize(remap.get)(Z)
    palette = ['#e8e8e8', '#9ecae1', '#6baed6', '#3182bd', '#08519c'][:len(order)]

    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    ax.pcolormesh(xs, [t - 273.15 for t in Tgrid], Z, shading='auto',
                  cmap=ListedColormap(palette),
                  norm=BoundaryNorm(range(len(order) + 1), len(order)))
    for sym, lad, th, color in data:
        xv = [x for x, y in zip(xs, lines[sym]) if isinstance(y, float)]
        yv = [y for y in lines[sym] if isinstance(y, float)]
        if xv:
            ax.plot(xv, yv, '-', lw=2.5, color=color, label='%s  Oh → Td' % sym)
    for x, lab in S.MEAS:
        ax.plot(x, 60, 'o', ms=6, mfc='none', mec='k', mew=1.5, zorder=5)
        ax.annotate(lab, (x, 60), textcoords='offset points', xytext=(0, 7),
                    ha='center', fontsize=7)
    ax.axhline(60, color='k', ls=':', lw=1)
    ax.axvline(0.333, color='orange', ls='--', lw=1.2)
    ax.annotate('conventional 1:2', (0.325, 124), rotation=90, va='top',
                ha='right', color='orange', fontsize=8)
    ax.set_xlim(0.05, 0.95); ax.set_ylim(25, 125)
    ax.set_xlabel('x(ChCl)   [ChCl : citric acid, water-free basis]')
    ax.set_ylabel('temperature / °C')
    ax.set_title('separation regime map%s' % ('  (with water, x_H2O=0.80)' if use_water else ''))
    h = [Patch(facecolor=palette[i], edgecolor='0.6', label=label(t))
         for i, t in enumerate(order)]
    ax.legend(handles=h + ax.get_legend_handles_labels()[0],
              loc='center left', bbox_to_anchor=(1.01, 0.5), fontsize=8, frameon=False)
    fig.tight_layout()
    p = os.path.join(HERE, 'window_map.png')
    fig.savefig(p, dpi=160, bbox_inches='tight')
    print('wrote', p)
    print('\n한계: Ni 보정값을 Co·Mn 에 전용 (용매 오차가 금속 무관이라는 가정).')
    print('      a_Cl 은 이온쌍 대리변수. 절대값 인용 금지.')
    print('      Mn 은 d5 고스핀이라 d-d 가 전부 스핀금지 -> 분광 검증 불가, dG 만 유효.')


if __name__ == '__main__':
    main()
