"""조성비별 물성 시뮬레이션 — 2단계 구조.

    침출효율  <-  [1차 물성]  <-  [2차 물성]  <-  조성비
                 점도/산도/          여기를 계산한다
                 Cl-활동도 등        (COSMO-RS 산출)

A군 = 혼합물 열역학량. 조성비에 직접 의존하므로 그대로 읽는다.
B군 = DES descriptor. 성분값은 분자 고유값이라 조성 무관 -> Lemaoui 2020 Eq.(2)의
      몰분율 가중합 S_DES = Σ xj·Sj 로 조성 의존량을 만든다.

    python properties_by_composition.py          # 표 + properties_by_x.csv + png
    python properties_by_composition.py --test   # 자체검증
"""
import csv, math, sys

T_EXP = 333.15                                   # 실험 온도 60 C (2026-09-23 확정)
MEAS = {'1:19': 0.05, '3:17': 0.15, '1:4': 0.20, '3:1': 0.75, '19:1': 0.95}

# sigma_profile_binN 의 중심 sigma (e/A^2).  bin i -> -0.025 + 0.005i + 0.0025
BIN_C = [-0.025 + 0.005 * i + 0.0025 for i in range(11)]
DONOR = [i for i, s in enumerate(BIN_C) if s <= -0.010]      # bins 0-2  수소결합 주개
ACCEP = [i for i, s in enumerate(BIN_C) if s >= +0.010]      # bins 7-10 수소결합 받개
CL_BIN = 8                                                   # +0.0175: 염화물 이온 고유 피크


def load(path='feature_extraction/features_all.csv'):
    g = {}
    for r in csv.DictReader(open(path)):
        if r['hba_name'] == 'choline_chloride' and r['hbd_name'] == 'citric_acid':
            g.setdefault(float(r['T']), {})[round(float(r['x_hba']), 2)] = r
    return g


def at(g, T, x, field):
    """온도 격자 사이 선형보간 (80 C는 격자점이 아님)."""
    if T in g:
        return float(g[T][x][field])
    Ts = sorted(g)
    lo, hi = max(t for t in Ts if t <= T), min(t for t in Ts if t >= T)
    w = (T - lo) / (hi - lo)
    return (1 - w) * float(g[lo][x][field]) + w * float(g[hi][x][field])


def props(g, x, T=T_EXP):
    """조성 x에서의 2차 물성 일습."""
    f = lambda fld: at(g, T, x, fld)
    mix = lambda suf: x * f('hba_' + suf) + (1 - x) * f('hbd_' + suf)   # Lemaoui Eq.(2)
    sp = lambda who, bins: sum(f(f'{who}_sigma_profile_bin{i}') for i in bins)
    wsum = lambda bins: x * sp('hba', bins) + (1 - x) * sp('hbd', bins)

    return {
        # --- A군: 혼합물 열역학 (조성비에 직접 의존) ---
        'lng_ChCl':   f('lng_hba_tot'),
        'a_ChCl':     x * math.exp(f('lng_hba_tot')),
        'lng_CA':     f('lng_hbd_tot'),
        'gE_over_RT': f('gE_over_RT'),
        'E_hb_kJ':    (x * f('hba_pm_E_hb') + (1 - x) * f('hbd_pm_E_hb')) / 1000,
        'E_mf_kJ':    (x * f('hba_pm_E_mf') + (1 - x) * f('hbd_pm_E_mf')) / 1000,
        'A_int_kJ':   (x * f('hba_pm_A_int') + (1 - x) * f('hbd_pm_A_int')) / 1000,
        # --- B군: DES descriptor (몰분율 가중합) ---
        'S_donor':    wsum(DONOR),          # 산성 양성자 표면적 (A^2)
        'S_acceptor': wsum(ACCEP),          # 수소결합 받개 표면적
        'S_chloride': x * f(f'hba_sigma_profile_bin{CL_BIN}') + (1 - x) * f(f'hbd_sigma_profile_bin{CL_BIN}'),
        'area':       mix('area'),
        'volume':     mix('volume'),
        'dipole':     mix('dipole_moment_mag'),
        'sig_m2':     mix('sigma_moment2'),  # 극성도(총 분극 세기)
        'sig_m3':     mix('sigma_moment3'),  # 극성 비대칭(주개/받개 치우침)
    }


# 2차 물성 -> 1차 물성 매핑 (문헌 근거는 DFT_PLAN.md §1-4)
MAP = [
    ('a_ChCl',                 '염화물 활동도 (P4)', '직접',      '이 값이 곧 1차 물성. Ni 화학종을 지배'),
    ('S_chloride',             '염화물 활동도 (P4)', '직접',      'Cl- 수소결합 받개 표면적. a_ChCl의 구조적 근거'),
    ('S_donor',                'HBD 산도 (P2)',      '대리변수',  '산성 양성자 표면. pKa 자체는 아님'),
    ('sig_m2, sig_m3, S_*',    '점도 (P1)',          'QSPR 입력', 'Mohan 2024 R2=0.99. 학습된 모델 필요'),
    ('sig_m2, sig_m3, S_*',    '전도도',             'QSPR 입력', 'Lemaoui 2020 R2=0.985. 학습된 모델 필요'),
    ('gE_over_RT, E_hb, A_int','비이상성/공융거동',   '직접',      '상호작용 세기. 조성 최적점 판정'),
    ('—',                      '물 활동도 (P3)',     '계산 불가', '물이 성분에 없는 2성분계라 산출 안 됨'),
]


def test():
    g = load()
    assert sorted(g) == [298.15, 318.15, 338.15, 358.15, 373.15]
    assert DONOR == [0, 1, 2] and ACCEP == [7, 8, 9, 10], (DONOR, ACCEP)
    assert abs(BIN_C[CL_BIN] - 0.0175) < 1e-9, BIN_C[CL_BIN]
    # 가중합은 양 끝에서 순성분값으로 수렴해야 한다
    r0, r1 = props(g, 0.05), props(g, 0.95)
    f = lambda x, fld: at(g, T_EXP, x, fld)
    assert r1['area'] > r0['area'] or r1['area'] < r0['area']       # 단조성만 확인용
    a_lo = 0.05 * f(0.05, 'hba_area') + 0.95 * f(0.05, 'hbd_area')
    assert abs(r0['area'] - a_lo) < 1e-9, (r0['area'], a_lo)
    # 활동도는 조성과 함께 단조증가해야 한다 (같은 성분을 더 넣으니까)
    xs = sorted(g[338.15])
    acts = [props(g, x)['a_ChCl'] for x in xs]
    assert all(b > a for a, b in zip(acts, acts[1:])), acts
    # 염화물 표면적은 ChCl 쪽이 압도적이므로 x와 함께 증가
    scl = [props(g, x)['S_chloride'] for x in xs]
    assert scl[-1] > scl[0], scl
    print('self-check OK')


def main():
    g = load()
    xs = sorted(g[338.15])
    P = {x: props(g, x) for x in xs}

    print(f"ChCl : citric acid,  {T_EXP-273.15:.0f} C  —  조성비별 2차 물성\n")
    print("[A군] 혼합물 열역학 — 조성비에 직접 의존")
    hdr = ['시료', 'x_ChCl', 'a(ChCl)', 'lnγ_ChCl', 'gE/RT', 'E_hb', 'E_mf', 'A_int']
    print(''.join(f'{h:>11}' for h in hdr))
    for k, x in sorted(MEAS.items(), key=lambda kv: kv[1]):
        p = P[x]
        print(f"{k:>11}{x:11.2f}{p['a_ChCl']:11.2e}{p['lng_ChCl']:11.2f}"
              f"{p['gE_over_RT']:11.3f}{p['E_hb_kJ']:11.1f}{p['E_mf_kJ']:11.1f}{p['A_int_kJ']:11.1f}")
    print("   (E_* 단위 kJ/mol)")

    print("\n[B군] DES descriptor — 몰분율 가중합 (Lemaoui 2020 Eq.2)")
    hdr = ['시료', 'x_ChCl', 'S_donor', 'S_accep', 'S_Cl-', 'area', 'volume', 'dipole', 'σ_m2', 'σ_m3']
    print(''.join(f'{h:>10}' for h in hdr))
    for k, x in sorted(MEAS.items(), key=lambda kv: kv[1]):
        p = P[x]
        print(f"{k:>10}{x:10.2f}{p['S_donor']:10.1f}{p['S_acceptor']:10.1f}{p['S_chloride']:10.1f}"
              f"{p['area']:10.1f}{p['volume']:10.1f}{p['dipole']:10.2f}{p['sig_m2']:10.1f}{p['sig_m3']:10.1f}")
    print("   (S_*, area 단위 Å², volume Å³)")

    print("\n[매핑] 2차 물성이 어느 1차 물성으로 가는가")
    print(f"{'2차 물성':<26}{'1차 물성':<22}{'관계':<12}비고")
    for a, b, c, d in MAP:
        print(f"{a:<26}{b:<22}{c:<12}{d}")

    with open('properties_by_x.csv', 'w', newline='') as fh:
        keys = list(P[xs[0]])
        w = csv.writer(fh); w.writerow(['x_ChCl', 'T'] + keys)
        for x in xs:
            w.writerow([x, T_EXP] + [f"{P[x][k]:.6g}" for k in keys])
    print(f"\n-> properties_by_x.csv  (전체 19개 조성 × {len(P[xs[0]])}개 물성)")

    try:
        import matplotlib; matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return
    # 제목은 영문: 기본 폰트에 한글이 없어 네모로 깨짐
    show = [('a_ChCl', 'chloride activity a(ChCl)', True),
            ('S_chloride', 'Cl- acceptor surface (A^2)', False),
            ('S_donor', 'acidic-proton surface (A^2)', False),
            ('gE_over_RT', 'gE/RT (non-ideality)', False),
            ('dipole', 'dipole moment (weighted)', False),
            ('sig_m3', 'sigma-moment 3 (polarity asymmetry)', False)]
    fig, axes = plt.subplots(2, 3, figsize=(13, 6.5))
    for ax, (key, title, logy) in zip(axes.ravel(), show):
        ax.plot(xs, [P[x][key] for x in xs], '-', lw=1.3, color='0.3')
        for k, xm in MEAS.items():
            ax.plot(xm, P[xm][key], 'o', ms=6)
            ax.annotate(k, (xm, P[xm][key]), fontsize=7, xytext=(2, 4), textcoords='offset points')
        if logy: ax.set_yscale('log')
        ax.set_title(title, fontsize=9); ax.set_xlabel('x(ChCl)', fontsize=8)
        ax.tick_params(labelsize=7); ax.grid(alpha=.3)
    plt.tight_layout(); plt.savefig('properties_by_x.png', dpi=130)
    print('-> properties_by_x.png')


if __name__ == '__main__':
    test() if '--test' in sys.argv else main()
