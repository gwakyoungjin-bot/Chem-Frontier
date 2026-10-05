"""금속별 독립 앵커로 보정한 분리 영역도.

이전 방식의 문제 (2026-09-27):
  Ni 하나에만 앵커를 걸고 그 보정값을 Co·Mn 에 전용했다. 결과가 물리적으로 틀렸다
  (Co·Mn 이 25 C 에서도 100% 사면체). 그리고 Ni 앵커조차 **틀린 조건**에 걸려 있었다
  — Hartley 는 ChCl:EG + 물 몰분율 0.045 인데, 우리 계(ChCl:CA + 물 0.80)의 활동도에
  걸고 있었다. HBD 효과 48.7배 + 물 효과 18~168배가 그대로 오차로 들어갔다.

이 스크립트가 고치는 것:
  1. 앵커 활동도를 **그 논문 자신의 계**에서 계산해 쓴다 (아래 ANCHORS)
     -> 모델이 조성이 아니라 활동도에 걸리므로 HBD·물 불일치가 동시에 해소된다
  2. **금속마다 자기 앵커로 shift 를 따로 맞춘다** (전용하지 않는다)
  3. Busato 2022 로 Ni 를 독립 검증한다 (앵커로 쓰지 않고 통과 여부만 본다)

alpha(이온쌍 대리변수 지수)는 쓰지 않는다:
  금속당 앵커가 하나뿐이라 자유 파라미터도 하나여야 한다. 파라미터를 늘리면
  "맞춘 것"이지 "예측"이 아니게 된다.

실행:
    python anchored_map.py
"""
import csv, math, os, sys

import numpy as np
import speciation as S

HERE = os.path.dirname(os.path.abspath(__file__))

# --- 앵커: 활동도는 각 논문 자신의 계에서 COSMO-RS 로 계산한 값 ----------------
# (재현: dft_ni_complexes/hbd_effect.py 및 ternary_water.py 와 같은 방식)
ANCHORS = {
    'Ni': dict(a=1.8325e-01, T=368.15, f=0.50,
               src='Hartley 2025 (10.1021/acs.jpcc.5c05771): ChCl:EG 1:2, x_H2O=0.045, '
                   'Oh->Td 전환 90~100 C -> 95 C 에서 f_Td=0.50'),
    'Co': dict(a=1.3908e-02, T=298.15, f=0.40,
               src='Mannucci 2026 (10.1021/acs.inorgchem.6c01344): ChCl/CoCl2.6H2O 1:2, '
                   'x_H2O=0.80, XAS+MCR 로 Td 40% (실온). 활동도는 ChCl+H2O 이성분 근사'),
}
# Mn 은 앵커가 없다. Co 의 shift 를 전용하되 그림에서 점선 + '미보정' 으로 표기한다.
MN_BORROWS = 'Co'

# 독립 검증용 (앵커로 쓰지 않음)
# Busato 2022: NiCl2.6H2O:urea 1:3.5, 실온, 염화물 극도로 풍부한데도 Ni 완전 팔면체.
# 그 계의 a_ChCl 을 우리 모델로 계산할 수 없으므로(ChCl 이 없는 금속-DES),
# '우리 계에서 가장 염화물이 많은 조건' 으로 대신 확인한다 -> 거기서도 f_Td 가 작아야 한다.
BUSATO_NOTE = ('Busato 2022 (10.1021/acs.inorgchem.2c00864): 염화물 과잉 MDES 실온에서 '
               'Ni 완전 팔면체(f_Td=0). 모델이 실온에서 f_Td<0.05 면 통과.')

COLOR = {'Ni': 'tab:blue', 'Co': 'tab:red', 'Mn': 'tab:green'}


def ladder(sym):
    return [('%s_aq6' % sym,    0, 6, 'Oh'),
            ('%sCl1_aq5' % sym, 1, 5, 'Oh'),
            ('%sCl2_aq4' % sym, 2, 4, 'Oh'),
            ('%sCl3_aq1' % sym, 3, 1, 'Td'),
            ('%sCl4' % sym,     4, 0, 'Td')]


DIRS = {'Ni': HERE, 'Co': os.path.join(HERE, 'cobalt'),
        'Mn': os.path.join(HERE, 'manganese')}


def harvest(sym):
    """준조화 보정된 (dH, dS). 참조종 H2O/Cl- 는 Ni 폴더 것 공용."""
    d, lad = DIRS[sym], ladder(sym)
    imag = []
    for name, n, m, g in lad:
        c, mx = S.count_imaginary(os.path.join(d, '%s_freq.out' % name))
        if c:
            imag.append('%s: 허수 %d개 (최대 %.0f cm-1)' % (name, c, mx))

    def dG(T):
        G = {}
        for name, n, m, g in lad:
            v = S.read_gibbs_eff(os.path.join(d, '%s_freq.out' % name), T)
            if v is not None:
                G[name] = v * S.HARTREE
        gw = S.read_gibbs_eff(os.path.join(HERE, 'H2O_freq.out'), T)
        ecl = S.read_sp_energy(os.path.join(HERE, 'Cl_ion_sp.out'))
        ref = '%s_aq6' % sym
        if ref not in G or gw is None or ecl is None:
            return None
        gw *= S.HARTREE
        gcl = S.atom_gibbs(ecl, 34.969, T)
        return {k: G[k] + (6 - m) * gw - G[ref] - n * gcl
                for k, n, m, g in lad if k in G}

    lo, hi = dG(S.T_LO), dG(S.T_HI)
    if lo is None or hi is None:
        return None, imag
    out = {}
    for k in lo:
        if k in hi:
            dS = -(hi[k] - lo[k]) / (S.T_HI - S.T_LO)
            out[k] = (lo[k] + S.T_LO * dS, dS)
    return out, imag


def f_td(th, lad, a, T, shift):
    num, tot = {}, 0.0
    ln_a = math.log(max(a, 1e-300))
    for name, n, m, g in lad:
        if name not in th:
            continue
        dH, dS = th[name]
        # 물 활동도는 3성분 계산에 이미 반영돼 있으므로 여기서 따로 곱하지 않는다
        v = math.exp(max(-700.0, min(700.0,
            -((dH - T * dS) + n * shift) / (S.R_KJ * T) + n * ln_a)))
        num[name] = v
        tot += v
    if tot <= 0:
        return 0.0
    td = {k for k, _, _, g in lad if g == 'Td'}
    return sum(v for k, v in num.items() if k in td) / tot


def fit_shift(th, lad, a, T, target):
    """그 금속의 앵커 한 점을 정확히 재현하는 shift."""
    g = lambda s: f_td(th, lad, a, T, s) - target
    lo, hi = -600.0, 600.0
    if g(lo) * g(hi) > 0:
        return None
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if g(lo) * g(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def main():
    grid = S.load_activity(True)        # 3성분 ChCl:CA:H2O, x_H2O=0.80 (우리 실험 조건)

    th, shift, imag_all = {}, {}, []
    for sym in ('Ni', 'Co', 'Mn'):
        t, im = harvest(sym)
        if t is None:
            print('%s: 결과 없음 -> 건너뜀' % sym)
            continue
        th[sym] = t
        imag_all += im
    if 'Ni' not in th:
        sys.exit('Ni 결과가 없습니다.')

    print('=' * 74)
    print('금속별 독립 앵커 보정')
    print('=' * 74)
    for sym in list(th):
        src = ANCHORS.get(sym)
        if src is None:
            continue
        s = fit_shift(th[sym], ladder(sym), src['a'], src['T'], src['f'])
        if s is None:
            print('  %s: 앵커를 맞추는 shift 없음 -> 제외' % sym)
            del th[sym]
            continue
        shift[sym] = s
        print('  %s  shift = %+8.1f kJ/mol per Cl-' % (sym, s))
        print('      앵커: %s' % src['src'])
        print('      확인: a=%.3e, T=%.0f C -> f_Td=%.3f (목표 %.2f)'
              % (src['a'], src['T'] - 273.15,
                 f_td(th[sym], ladder(sym), src['a'], src['T'], s), src['f']))
    if 'Mn' in th and MN_BORROWS in shift:
        shift['Mn'] = shift[MN_BORROWS]
        print('  Mn  shift = %+8.1f kJ/mol per Cl-  ** %s 것을 빌림 — 앵커 없음(미보정) **'
              % (shift['Mn'], MN_BORROWS))

    if imag_all:
        print('\n⚠ 허수진동수 (준조화 보정 적용됨):')
        for w in imag_all:
            print('   ', w)

    # --- 독립 검증: Busato ------------------------------------------------
    print('\n독립 검증 (앵커로 쓰지 않음)')
    print('  %s' % BUSATO_NOTE)
    a_rt = S.a_cl(grid, 0.95, 298.15)
    v = f_td(th['Ni'], ladder('Ni'), a_rt, 298.15, shift['Ni'])
    print('  -> 우리 계 최대 염화물(x=0.95), 25 C 에서 f_Td(Ni) = %.4f  [%s]'
          % (v, '통과' if v < 0.05 else '불통과'))

    # --- 전환선 ------------------------------------------------------------
    xs = [round(0.05 * i, 2) for i in range(1, 20)]
    Ts = [298.15 + 2.5 * i for i in range(0, 41)]
    lines = {}
    for sym in th:
        ys = []
        for x in xs:
            f0 = f_td(th[sym], ladder(sym), S.a_cl(grid, x, Ts[0]), Ts[0], shift[sym])
            if f0 >= 0.5:
                ys.append(None); continue          # 범위 아래에서 이미 전환
            prev, t = None, 'high'
            for T in Ts:
                f = f_td(th[sym], ladder(sym), S.a_cl(grid, x, T), T, shift[sym])
                if prev and prev[1] < 0.5 <= f:
                    T0, f0b = prev
                    t = T0 + (T - T0) * (0.5 - f0b) / (f - f0b) - 273.15
                    break
                prev = (T, f)
            ys.append(t)
        lines[sym] = ys

    syms = [s for s in ('Ni', 'Co', 'Mn') if s in th]
    print('\n전환선 (f_Td = 0.5 가 되는 온도, °C).  <25 = 범위 아래, >125 = 범위 위')
    print('  %-8s %s' % ('x(ChCl)', ' '.join('%-10s' % s for s in syms)))
    rows = []
    for i, x in enumerate(xs):
        cells, r = [], dict(x_ChCl=x)
        for s in syms:
            v = lines[s][i]
            txt = '<25' if v is None else ('>125' if v == 'high' else '%.0f' % v)
            cells.append('%-10s' % txt)
            r['T_%s_C' % s] = txt
        print('  %-8.2f %s' % (x, ' '.join(cells)))
        rows.append(r)
    out = os.path.join(HERE, 'window.csv')
    with open(out, 'w', newline='', encoding='utf-8') as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0])); wr.writeheader(); wr.writerows(rows)
    print('\nwrote', out)

    # --- 60 C 우리 시료 ----------------------------------------------------
    print('\n60 °C 측정 5조성의 사면체 분율')
    print('  %-6s %-6s %-11s %s' % ('시료', 'x', 'a_ChCl', ' '.join('%-9s' % s for s in syms)))
    frows = []
    for x, lab in S.MEAS:
        a = S.a_cl(grid, x, 333.15)
        v = [f_td(th[s], ladder(s), a, 333.15, shift[s]) for s in syms]
        print('  %-6s %-6.2f %-11.2e %s' % (lab, x, a, ' '.join('%-9.4f' % q for q in v)))
        r = dict(sample=lab, x_ChCl='%.2f' % x, a_ChCl='%.4e' % a)
        r.update({'f_Td_%s' % s: '%.4f' % q for s, q in zip(syms, v)})
        frows.append(r)
    # 화면에만 찍고 끝나면 원고가 이 수치를 인용할 근거가 남지 않는다 -> 파일로 고정
    fout = os.path.join(HERE, 'ftd_60C.csv')
    with open(fout, 'w', newline='', encoding='utf-8') as f:
        wr = csv.DictWriter(f, fieldnames=list(frows[0])); wr.writeheader(); wr.writerows(frows)
    print('wrote', fout)

    # 보정값도 같이 남긴다 (원고 Q5 가 인용하는 값)
    sout = os.path.join(HERE, 'shifts.csv')
    with open(sout, 'w', newline='', encoding='utf-8') as f:
        wr = csv.writer(f); wr.writerow(['metal', 'shift_kJ_per_Cl', 'anchor_src'])
        for s in syms:
            src = ANCHORS.get(s, {}).get('src', 'borrowed from %s' % MN_BORROWS)
            wr.writerow([s, '%.1f' % shift[s], src])
    print('wrote', sout)

    # --- 그림 -----------------------------------------------------------
    import plot_regime as PR, plot_styles as PS
    Tc = [t - 273.15 for t in Ts]
    sub = ('DFT (r2SCAN-3c, quasi-harmonic) + COSMO-RS activities at the experimental' + chr(10) +
           'water content (x H2O = 0.80).  Each metal anchored to its own measurement:' + chr(10) +
           'Ni = Hartley 2025,  Co = Mannucci 2026.  Mn has no anchor (dashed).')
    d = os.path.join(HERE, "figure_styles")
    os.makedirs(d, exist_ok=True)

    # C: 검증된 범주형 (기본 산출물)
    PR.draw(xs, Tc, lines, S.MEAS, os.path.join(HERE, 'window_map.png'), sub)
    PR.draw(xs, Tc, lines, S.MEAS, os.path.join(d, 'C_categorical.png'), sub)

    # B: 같은 구조 + viridis 이산 색
    PS.regime_viridis(xs, Tc, lines, os.path.join(d, 'B_regime_viridis.png'), sub)

    # A: 금속별 f_Td 연속장 3패널 (viridis 본래 용도)
    fields = {}
    for sym in th:
        fields[sym] = [[f_td(th[sym], ladder(sym), S.a_cl(grid, x, T), T, shift[sym])
                        for x in xs] for T in Ts]
    PS.small_multiples(xs, Tc, fields, os.path.join(d, 'A_viridis_panels.png'))


if __name__ == '__main__':
    main()
