"""DFT 자유에너지 + COSMO-RS 활동도 -> 조성-온도 가이드라인 지도.

    python speciation.py            # ORCA 출력 수확 -> speciation.csv, guideline_map.png
    python speciation.py --test     # ORCA 없이 자체검증 (합성 dG)
    python speciation.py --awater 0.3 --ftdmax 0.02    # 가정값 민감도 확인

논리:
    Ni(H2O)6(2+) + n Cl- -> NiCl_n(H2O)_m + (6-m) H2O
    dG_n = G(NiCl_n) + (6-m)G(H2O) - G(Ni_aq6) - n G(Cl-)        <- DFT(Freq)
    f_n  ∝ exp(-(dG_n + n*shift)/RT) * a_Cl^(n*alpha) / a_W^(6-m) <- COSMO-RS 활동도
    f_Td = (f_NiCl3_aq1 + f_NiCl4) / sum(f)

왜 보정(shift, alpha)이 필요한가:
  * shift — 전하가 +2에서 -2까지 바뀌는 이온 반응의 절대 dG 는 암시적 용매로
    수십 kJ/mol 틀린다. 오차 대부분이 Cl- 한 개당 용매화 에너지이므로
    "Cl- 개수에 비례하는 상수" 하나로 흡수된다.
  * alpha — a_Cl 을 최대 4제곱하므로 이온쌍 근사 대리변수의 오차가 증폭된다.
    실측: a(x=0.95,60C)/a(x=0.333,95C) = 240배 -> 4제곱하면 3.2e9배.
    이러면 dG 가 무엇이든 염 과잉 조성이 f_Td=1 로 포화되어 지도가 무의미해진다.
    대리변수와 자유 Cl- 활동도를  a_Cl- = C * (a_ChCl)^alpha  로 두고 alpha 를 고정한다
    (C 는 shift 에 흡수).

두 파라미터를 고정하는 실측 앵커 2개:
  1. Hartley et al., J Phys Chem C 2025, 10.1021/acs.jpcc.5c05771
     ChCl계 DES에서 Ni(II) 팔면체->사면체 전환이 90~100 C (EXAFS 실측)
     -> 95 C, x=0.333 에서 f_Td = 0.50
  2. 우리 UV-Vis 미검출 (조성비_논리흐름.md)
     3sigma = 0.075 A, eps(Td) ~ 180 -> Td < 0.42 mM
     -> 60 C, x=0.95 에서 f_Td <= 0.05 (명목값)

  ⚠ 이건 보정이지 검증이 아니다. 미검출을 모델에 먹였으므로
    "모델이 미검출을 맞혔다"고 말하면 순환논증이다.
    말할 수 있는 것: 두 실측 조건을 동시에 만족하는 모델로 '다른 조건'을 예측한다.

한계 (보고서에 반드시 기재):
  1. a_Cl 은 자유 Cl- 활동도가 아니라 ChCl 중성 이온쌍 근사의 대리변수.
     절대값 인용 금지, 순서와 자릿수만 사용 (조성비_논리흐름.md §2).
  2. 물은 ChCl:CA 이성분 COSMO-RS 모델의 성분이 아니다 -> a_W 는 입력 가정값(기본 1.0).
  3. 앵커1 문헌은 ChCl:EG / ChCl:urea 계. 우리 ChCl:CA 와 HBD 가 다르다.
  4. 앵커2 의 f_Td 상한은 총 Ni 농도를 모르므로 명목값이다.
  5. 암시적 용매만 사용 -> DES 2차 배위권(콜린 양이온, 시트르산) 없음. 시트르산 착물 미포함.
"""
import csv, math, os, re, sys

R_KJ = 8.314462618e-3          # kJ/mol/K
HARTREE = 2625.4996            # kJ/mol
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# 사다리. (이름, Cl 개수 n, 물 개수 m, 기하)
LADDER = [('Ni_aq6',    0, 6, 'Oh'),
          ('NiCl1_aq5', 1, 5, 'Oh'),
          ('NiCl2_aq4', 2, 4, 'Oh'),
          ('NiCl3_aq1', 3, 1, 'Td'),
          ('NiCl4',     4, 0, 'Td')]

ANCHOR1_T, ANCHOR1_X, ANCHOR1_F = 368.15, 0.333, 0.50   # Hartley 2025, 95 C
ANCHOR2_T, ANCHOR2_X = 333.15, 0.95                     # 우리 UV-Vis, 60 C, 19:1
FTD_MAX = 0.05                                          # 미검출 -> 분율 명목 상한

MEAS = [(0.05, '1:19'), (0.15, '3:17'), (0.20, '1:4'), (0.75, '3:1'), (0.95, '19:1')]

# 시트르산 착물 (선택). citrate/ 결과가 없으면 통째로 무시되고 나머지는 그대로 동작한다.
#   [Ni(H2O)6]2+ + H3Cit -> [Ni(H3Cit)(H2O)4]2+ + 2 H2O
# ** 이 반응은 양변의 전하가 같다(+2 -> +2). **
# Cl 사다리(+2 -> -2)와 달리 이온 용매화 오차가 거의 상쇄되므로
# shift 보정을 적용하지 않는다. 즉 Cl 사다리보다 신뢰도가 높은 숫자다.
# (이름, 물 개수 m, 기하, 시트르산 개수)
EXTRA = [('NiCit_aq4', 4, 'Oh', 1)]
CIT_DIR = 'citrate'


# ---------------------------------------------------------------- ORCA 파싱
def read_gibbs(path, T):
    """ORCA freq 출력에서 지정 온도의 Final Gibbs free energy (Hartree)."""
    if not os.path.exists(path):
        return None
    txt = open(path, encoding='utf-8', errors='ignore').read()
    blocks = re.split(r'THERMOCHEMISTRY AT\s+([\d.]+)\s*K', txt)
    for i in range(1, len(blocks) - 1, 2):
        if abs(float(blocks[i]) - T) < 0.5:
            m = re.search(r'Final Gibbs free energy\s*\.*\s*(-?\d+\.\d+)', blocks[i + 1])
            if m:
                return float(m.group(1))
    return None


QH = True      # 준조화 보정 사용 (--raw 로 끔). thermo_qh.py 참조.


def read_gibbs_eff(path, T):
    """준조화 보정을 적용한 Gibbs (Hartree).

    ORCA 는 허수모드를 분배함수에서 빼버려, 허수가 있는 종의 엔트로피가 과소평가된다.
    그게 기준물질이면 오차가 모든 dG 에 실린다 (Co·Mn 이 25 C 에서도 100% 사면체로
    튀었던 원인). thermo_qh 가 |nu|<100 cm-1 모드를 전부 100 cm-1 로 치환해
    모든 종을 같은 규칙으로 다룬다.
    """
    g = read_gibbs(path, T)
    if g is None or not QH:
        return g
    import thermo_qh as Q
    c, _, _ = Q.correction(path, T)
    return g if c is None else g + c / HARTREE


def read_sp_energy(path):
    """단일점 전자에너지 (Hartree)."""
    if not os.path.exists(path):
        return None
    hits = re.findall(r'FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)',
                      open(path, encoding='utf-8', errors='ignore').read())
    return float(hits[-1]) if hits else None


def count_imaginary(path):
    """허수진동수 개수와 가장 큰 크기(cm-1). Freq 결과 건전성 점검용."""
    if not os.path.exists(path):
        return None, 0.0
    txt = open(path, encoding='utf-8', errors='ignore').read()
    m = re.search(r'VIBRATIONAL FREQUENCIES(.*?)(NORMAL MODES|$)', txt, re.S)
    if not m:
        return None, 0.0
    v = [float(x) for x in re.findall(r':\s+(-?\d+\.\d+)\s*cm\*\*-1', m.group(1))]
    neg = [x for x in v if x < -1.0]
    return len(neg), (abs(min(neg)) if neg else 0.0)


def atom_gibbs(E_hartree, mass_amu, T, mult=1):
    """단원자 이온의 G (kJ/mol). 진동 없음 -> 병진(Sackur-Tetrode) + 전자 축퇴만."""
    kB, h, NA = 1.380649e-23, 6.62607015e-34, 6.02214076e23
    m = mass_amu / 1000.0 / NA
    V = kB * T / 101325.0                                   # 1 atm 기준 분자당 부피
    q = (2 * math.pi * m * kB * T / h ** 2) ** 1.5 * V
    S = R_KJ * (math.log(q) + 2.5) + R_KJ * math.log(mult)  # kJ/mol/K
    H = 2.5 * R_KJ * T                                      # 3/2RT(병진) + RT(pV)
    return E_hartree * HARTREE + H - T * S


# ------------------------------------------------------- COSMO-RS 활동도
def load_activity(water=False):
    """활동도 표 -> {T: {x_ChCl(건조기준): (a_ChCl, a_CA)}}

    water=False : features_all.csv (ChCl:CA 2성분, 물 없음)
    water=True  : ternary_grid.csv  (ChCl:CA:H2O 3성분, x_H2O=0.80 = 실험 조건)

    왜 물 버전이 필요한가 (2026-09-27):
      레시피가 DES 3000 mg + 물 1270 uL 라 **물 몰분율이 0.77~0.82** 다.
      2성분 모델에는 물이 없어서 이 실험계를 기술하지 못한다.
      물을 넣으면 조성 간 활동도 대비가 5e7배 -> 9.6e4배로 약 500배 압축된다.
      즉 조성비의 지렛대가 2성분 모델이 말하는 것보다 훨씬 작다.
    """
    g = {}
    if water:
        cands = [os.path.join(d, 'ternary_grid.csv') for d in (HERE, ROOT)]
        p = next((c for c in cands if os.path.exists(c)), None)
        if p is None:
            raise SystemExit('ternary_grid.csv 가 없습니다. 먼저:'
                             '  python ternary_water.py --grid')
        for r in csv.DictReader(open(p, encoding='utf-8')):
            if abs(float(r['x_H2O']) - 0.80) > 1e-9:
                continue
            g.setdefault(float(r['T_K']), {})[round(float(r['x_ChCl_dry']), 2)] = (
                float(r['a_ChCl']), float(r['a_CA']))
        if not g:
            raise SystemExit('ternary_grid.csv 에 x_H2O=0.80 행이 없습니다.')
        return g

    cands = [os.path.join(d, 'feature_extraction', 'features_all.csv') for d in (ROOT, HERE)]
    p = next((c for c in cands if os.path.exists(c)), None)
    if p is None:
        raise SystemExit('features_all.csv 를 찾을 수 없습니다. 찾아본 곳:\n  ' +
                         '\n  '.join(cands))
    for r in csv.DictReader(open(p, encoding='utf-8')):
        if r['hba_name'] == 'choline_chloride' and r['hbd_name'] == 'citric_acid':
            x = round(float(r['x_hba']), 2)
            g.setdefault(float(r['T']), {})[x] = (
                x * math.exp(float(r['lng_hba_tot'])),
                (1.0 - x) * math.exp(float(r['lng_hbd_tot'])))
    return g


def _lerp(d, k):
    """1차원 딕셔너리 {키: 값} 선형보간. 격자 밖은 끝값 고정."""
    ks = sorted(d)
    k = min(max(k, ks[0]), ks[-1])
    lo = max(t for t in ks if t <= k)
    hi = min(t for t in ks if t >= k)
    if hi == lo:
        return d[lo]
    return d[lo] + (d[hi] - d[lo]) * (k - lo) / (hi - lo)


def _a(grid, x, T, idx):
    """활동도 보간. 자릿수가 크게 벌어지므로 **로그 공간**에서 선형보간한다.
    (선형공간에서 보간하면 작은 값 쪽이 뭉개진다)"""
    ln_at_T = {t: _lerp({k: math.log(max(v[idx], 1e-300)) for k, v in row.items()}, x)
               for t, row in grid.items()}
    return math.exp(_lerp(ln_at_T, T))


def a_cl(grid, x, T):
    """ChCl 활동도 (자유 Cl- 의 대리변수 — 절대값 인용 금지)."""
    return _a(grid, x, T, 0)


def a_ca(grid, x, T):
    """시트르산 활동도."""
    return _a(grid, x, T, 1)


# ------------------------------------------------------------- 화학종 분포
def fractions(th, aCl, aW, T, shift, alpha=1.0, extra=None, aCA=None):
    """반환: ({이름: 분율}, f_Td).

    th = {이름: (dH kJ/mol, dS kJ/mol/K)} — Ni_aq6 기준 반응의 엔탈피·엔트로피.
    dG(T) = dH - T*dS 로 온도를 다룬다. ORCA freq 를 두 온도(298/333 K)에서
    받아 van't Hoff 로 분리한 값이다.

    온도축이 물리적이려면 dS 가 반드시 필요하다. dG 를 한 온도에서만 읽고
    1/RT 로만 온도를 넣으면, shift 부호에 따라 "온도를 올리면 Td 가 준다"는
    거꾸로 된 지도가 나온다 (Hartley 전환과 정반대).
    Oh->Td 는 물을 여러 개 방출하므로 dS 가 크게 양수이고, 그게 전환의 실제 동력이다.

    ln a_Cl- = alpha * ln a_ChCl (상수는 shift 에 흡수).
    """
    num, tot = {}, 0.0
    ln_a = alpha * math.log(max(aCl, 1e-300))
    ln_w = math.log(max(aW, 1e-300))
    for name, n, m, geom in LADDER:
        if name not in th:
            continue
        dH, dS = th[name]
        lnK = -((dH - T * dS) + n * shift) / (R_KJ * T)
        v = math.exp(max(-700.0, min(700.0, lnK + n * ln_a - (6 - m) * ln_w)))
        num[name] = v
        tot += v
    if tot <= 0:
        return {k: 0.0 for k in num}, 0.0
    # 시트르산 착물: 전하중립 반응이라 shift 를 적용하지 않는다 (EXTRA 주석 참조)
    if extra and aCA is not None:
        ln_c = math.log(max(aCA, 1e-300))
        for name, m, geom, ncit in EXTRA:
            if name not in extra:
                continue
            dH, dS = extra[name]
            lnK = -(dH - T * dS) / (R_KJ * T)
            v = math.exp(max(-700.0, min(700.0, lnK + ncit * ln_c - (6 - m) * ln_w)))
            num[name] = v
            tot += v
        if tot <= 0:
            return {k: 0.0 for k in num}, 0.0

    f = {k: v / tot for k, v in num.items()}
    td = set(n for n, _, _, g in LADDER if g == 'Td') | set(
        n for n, _, g, _ in EXTRA if g == 'Td')
    return f, sum(v for k, v in f.items() if k in td)


def _bisect(fn, lo, hi, iters=200):
    """부호가 바뀌는 구간에서 근을 찾는다. 부호가 안 바뀌면 None."""
    flo = fn(lo)
    if flo * fn(hi) > 0:
        return None
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        fm = fn(mid)
        if flo * fm <= 0:
            hi = mid
        else:
            lo, flo = mid, fm
    return 0.5 * (lo + hi)


def solve_shift(th, grid, aW, alpha=1.0, extra=None):
    """앵커1(95 C, x=0.333 에서 f_Td=0.5)을 만족하는 shift."""
    a = a_cl(grid, ANCHOR1_X, ANCHOR1_T)
    c = a_ca(grid, ANCHOR1_X, ANCHOR1_T)
    return _bisect(lambda s: fractions(th, a, aW, ANCHOR1_T, s, alpha, extra, c)[1] - ANCHOR1_F,
                   -600.0, 600.0)


def calibrate(th, grid, aW, ftd_max=FTD_MAX, extra=None):
    """두 앵커로 (shift, alpha) 동시 결정.

    바깥: alpha 를 움직여 앵커2(60 C, x=0.95 에서 f_Td = ftd_max)를 맞춘다.
    안쪽: alpha 마다 앵커1 이 성립하도록 shift 를 푼다.
    반환 (shift, alpha). alpha 가 None 이면 앵커2 를 못 맞춘 것.
    """
    a2 = a_cl(grid, ANCHOR2_X, ANCHOR2_T)
    c2 = a_ca(grid, ANCHOR2_X, ANCHOR2_T)

    def resid(alpha):
        s = solve_shift(th, grid, aW, alpha, extra)
        if s is None:
            return None
        return fractions(th, a2, aW, ANCHOR2_T, s, alpha, extra, c2)[1] - ftd_max

    # 로그 간격으로 훑는다. 근이 아주 작은 alpha(1e-3 급)에 있을 수 있어서
    # 선형 격자(0.02 간격)로는 통째로 놓친다 — 실제로 놓쳤던 버그.
    scan = [(10 ** (-4 + 4.5 * i / 240.0), None) for i in range(241)]   # 1e-4 ~ 3.16
    scan = [(a, resid(a)) for a, _ in scan]
    scan = [(a, r) for a, r in scan if r is not None]
    if not scan:
        return None, None
    for (a1, r1), (a2_, r2) in zip(scan, scan[1:]):
        if r1 * r2 <= 0:
            f = lambda A: (resid(A) if resid(A) is not None else 1e9)
            alpha = _bisect(f, a1, a2_)
            if alpha is not None:
                return solve_shift(th, grid, aW, alpha, extra), alpha
    best = min(scan, key=lambda p: abs(p[1]))[0]
    return solve_shift(th, grid, aW, best, extra), None


# ------------------------------------------------------------------ 수확
T_LO, T_HI = 298.15, 333.15      # ORCA 인풋의 %freq Temp 두 점


def _dG_at(T):
    """온도 T 에서 dG (Ni_aq6 기준, kJ/mol) + 빠진 종. 없으면 (None, ...)."""
    G = {}
    for name, n, m, geom in LADDER:
        g = read_gibbs_eff(os.path.join(HERE, '%s_freq.out' % name), T)
        if g is not None:
            G[name] = g * HARTREE
    gw = read_gibbs_eff(os.path.join(HERE, 'H2O_freq.out'), T)
    ecl = read_sp_energy(os.path.join(HERE, 'Cl_ion_sp.out'))
    missing = [n for n, _, _, _ in LADDER if n not in G]
    if gw is None:
        missing.append('H2O')
    if ecl is None:
        missing.append('Cl_ion')
    if 'Ni_aq6' not in G or gw is None or ecl is None:
        return None, missing
    gw *= HARTREE
    gcl = atom_gibbs(ecl, 34.969, T)                 # 35Cl
    return ({name: G[name] + (6 - m) * gw - G['Ni_aq6'] - n * gcl
             for name, n, m, geom in LADDER if name in G}, missing)


def harvest_extra():
    """citrate/ 의 freq 결과 -> {이름: (dH, dS)}. 없으면 빈 dict (무시하고 진행)."""
    d = os.path.join(HERE, CIT_DIR)
    out, warn = {}, []
    gw = {T: read_gibbs_eff(os.path.join(HERE, 'H2O_freq.out'), T) for T in (T_LO, T_HI)}
    gref = {T: read_gibbs_eff(os.path.join(HERE, 'Ni_aq6_freq.out'), T) for T in (T_LO, T_HI)}
    glig = {T: read_gibbs_eff(os.path.join(d, 'H3Cit_freq.out'), T) for T in (T_LO, T_HI)}
    for name, m, geom, ncit in EXTRA:
        gc = {T: read_gibbs_eff(os.path.join(d, '%s_freq.out' % name), T) for T in (T_LO, T_HI)}
        if any(v is None for v in list(gc.values()) + list(glig.values())
               + list(gw.values()) + list(gref.values())):
            continue
        ni, mx = count_imaginary(os.path.join(d, '%s_freq.out' % name))
        if ni:
            warn.append('%s: 허수 %d개 (최대 %.0f cm-1)' % (name, ni, mx))
        dG = {}
        for T in (T_LO, T_HI):
            dG[T] = HARTREE * (gc[T] + (6 - m) * gw[T] - gref[T] - ncit * glig[T])
        dS = -(dG[T_HI] - dG[T_LO]) / (T_HI - T_LO)
        out[name] = (dG[T_LO] + T_LO * dS, dS)
    return out, warn


def harvest():
    """ORCA 출력 -> {이름: (dH, dS)}. 반환 (th, 빠진종, 허수진동수경고).

    두 온도의 dG 를 van't Hoff 로 분리한다:  dS = -(dG_hi - dG_lo)/(T_hi - T_lo),
    dH = dG_lo + T_lo*dS.  dH, dS 가 이 구간에서 일정하다고 가정 (표준 근사).
    지도를 125 C 까지 그리므로 외삽이 들어간다 — 한계로 기재할 것.
    """
    imag = []
    for name, n, m, geom in LADDER:
        p = os.path.join(HERE, '%s_freq.out' % name)
        ni, mx = count_imaginary(p)
        if ni:
            imag.append('%s: 허수 %d개 (최대 %.0f cm-1)' % (name, ni, mx))
    lo, missing = _dG_at(T_LO)
    hi, _ = _dG_at(T_HI)
    if lo is None or hi is None:
        return None, missing, imag
    th = {}
    for name in lo:
        if name in hi:
            dS = -(hi[name] - lo[name]) / (T_HI - T_LO)
            th[name] = (lo[name] + T_LO * dS, dS)
    return th, missing, imag


# ------------------------------------------------------------------ 출력
def write_outputs(th, shift, alpha, grid, aW, extra=None):
    xs = [round(0.05 * i, 2) for i in range(1, 20)]
    Ts = [298.15 + 5 * i for i in range(0, 21)]      # 25 ~ 125 C
    rows = []
    for T in Ts:
        for x in xs:
            a = a_cl(grid, x, T)
            f, ftd = fractions(th, a, aW, T, shift, alpha, extra, a_ca(grid, x, T))
            rows.append(dict(x_ChCl=x, T_C=round(T - 273.15, 1), a_ChCl='%.3e' % a,
                             f_Td=round(ftd, 4), **{k: round(v, 4) for k, v in f.items()}))
    out = os.path.join(HERE, 'speciation.csv')
    with open(out, 'w', newline='', encoding='utf-8') as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
        wr.writeheader(); wr.writerows(rows)
    print('wrote %s  (%d행)' % (out, len(rows)))

    try:
        import numpy as np, matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib 없음 -> 그림 생략 (csv 는 생성됨)')
        return

    Tc = [t - 273.15 for t in Ts]
    Z = np.array([[fractions(th, a_cl(grid, x, T), aW, T, shift, alpha,
                             extra, a_ca(grid, x, T))[1]
                   for x in xs] for T in Ts])
    fig, ax = plt.subplots(figsize=(8, 5.5))
    im = ax.pcolormesh(xs, Tc, Z, cmap='viridis', vmin=0, vmax=1, shading='auto')
    if Z.min() < 0.5 < Z.max():
        cs = ax.contour(xs, Tc, Z, levels=[0.5], colors='white', linewidths=2)
        ax.clabel(cs, fmt={0.5: 'Td 50%'}, fontsize=9)
    for x, lab in MEAS:
        ax.plot(x, 60, 'o', ms=7, mfc='none', mec='red', mew=1.8)
        ax.annotate(lab, (x, 60), textcoords='offset points', xytext=(0, 8),
                    ha='center', color='red', fontsize=8)
    ax.axvline(0.333, color='orange', ls='--', lw=1.5)
    ax.annotate('conventional 1:2', (0.325, 122), color='orange', fontsize=8,
                rotation=90, va='top', ha='right')
    ax.axhline(60, color='red', ls=':', lw=1.2)
    ax.annotate('our 60 C', (0.02, 61), color='red', fontsize=8)
    ax.axhspan(90, 100, color='w', alpha=0.18)
    ax.annotate('Hartley 2025: Oh->Td 90-100 C', (0.5, 95), color='w',
                fontsize=8, ha='center', va='center')
    ax.set_xlabel('x(ChCl)    [ChCl : citric acid]')
    ax.set_ylabel('temperature / C')
    ax.set_title('predicted tetrahedral Ni(II) fraction  (2-anchor calibrated)')
    fig.colorbar(im, ax=ax, label='f(Td)')
    fig.tight_layout()
    p = os.path.join(HERE, 'guideline_map.png')
    fig.savefig(p, dpi=160)
    print('wrote', p)


# ------------------------------------------------------------------ 검증
def test():
    """ORCA 없이 로직만 검증. 합성 (dH, dS): Cl 가 붙으면 엔탈피 불리,
    Oh->Td 에서 물을 여러 개 방출하므로 엔트로피가 크게 유리 (전환의 실제 동력)."""
    # Td 종은 물을 5~6개 방출하므로 dS 가 커야 한다. 200 J/mol/K 급이 정상.
    # (dS 가 작으면 60->95 C 사이에서 전환이 충분히 가팔라지지 않아 앵커2 를 못 맞춘다)
    th = {'Ni_aq6':    (0.0,   0.000),
          'NiCl1_aq5': (20.0,  0.010),
          'NiCl2_aq4': (45.0,  0.025),
          'NiCl3_aq1': (85.0,  0.225),
          'NiCl4':     (105.0, 0.265)}
    grid = load_activity()

    # 활동도는 x 에 대해 단조 증가, 격자점은 보간해도 그대로
    a = [a_cl(grid, x, 333.15) for x in (0.05, 0.2, 0.5, 0.75, 0.95)]
    assert all(a[i] < a[i + 1] for i in range(len(a) - 1)), a
    assert abs(a_cl(grid, 0.20, 338.15) - grid[338.15][0.20][0]) < 1e-9

    # 분율 합 = 1
    f, ftd = fractions(th, 1e-3, 1.0, 333.15, 0.0)
    assert abs(sum(f.values()) - 1.0) < 1e-12, sum(f.values())
    assert 0.0 <= ftd <= 1.0

    # Cl 활동도가 커지면 Td 분율 단조 증가 (물리적 방향성)
    prev = -1.0
    for lg in range(-8, 1):
        _, t = fractions(th, 10.0 ** lg, 1.0, 333.15, 0.0)
        assert t >= prev - 1e-12, (lg, t, prev)
        prev = t

    # 온도를 올리면 Td 증가 — dS 를 쓰기 때문에 성립 (dG 한 점만 쓰면 뒤집힐 수 있음)
    seq = [fractions(th, 1e-2, 1.0, T, 0.0)[1] for T in (298.15, 333.15, 368.15, 398.15)]
    assert all(seq[i] <= seq[i + 1] + 1e-12 for i in range(len(seq) - 1)), seq

    # alpha=1 단일앵커로는 앵커2 가 깨진다 = 두 번째 파라미터가 필요한 이유
    s1 = solve_shift(th, grid, 1.0, 1.0)
    _, bad = fractions(th, a_cl(grid, ANCHOR2_X, ANCHOR2_T), 1.0, ANCHOR2_T, s1, 1.0)
    assert bad > 0.9, bad

    # 두 앵커 동시 보정이 실제로 둘 다 만족하는가
    shift, alpha = calibrate(th, grid, 1.0)
    assert alpha is not None, '앵커2 를 맞추는 alpha 없음'
    _, t1 = fractions(th, a_cl(grid, ANCHOR1_X, ANCHOR1_T), 1.0, ANCHOR1_T, shift, alpha)
    _, t2 = fractions(th, a_cl(grid, ANCHOR2_X, ANCHOR2_T), 1.0, ANCHOR2_T, shift, alpha)
    assert abs(t1 - ANCHOR1_F) < 1e-4, t1
    assert abs(t2 - FTD_MAX) < 1e-3, t2

    # 보정 후 60 C 의 측정 5조성은 전부 Td 소수 (전환온도보다 한참 아래)
    for x, lab in MEAS:
        _, t60 = fractions(th, a_cl(grid, x, 333.15), 1.0, 333.15, shift, alpha)
        assert t60 <= FTD_MAX + 1e-6, (lab, t60)

    # 보정 후에도 온도 방향은 유지되어야 한다
    seq2 = [fractions(th, a_cl(grid, 0.95, T), 1.0, T, shift, alpha)[1]
            for T in (313.15, 333.15, 353.15, 373.15)]
    assert all(seq2[i] <= seq2[i + 1] + 1e-9 for i in range(len(seq2) - 1)), seq2

    # Cl- 병진 엔트로피가 상식 범위인가 (문헌 ~153 J/mol/K)
    S = (2.5 * R_KJ * 298.15 - atom_gibbs(0.0, 34.969, 298.15)) / 298.15 * 1000
    assert 145 < S < 165, S

    # 시트르산 종을 넣어도 분율 합이 1 이고, 시트르산 활동도가 커지면 그 분율이 는다
    ex = {'NiCit_aq4': (-15.0, 0.010)}
    f0, _ = fractions(th, 1e-3, 1.0, 333.15, 0.0, 1.0, ex, 1e-4)
    f1, _ = fractions(th, 1e-3, 1.0, 333.15, 0.0, 1.0, ex, 1e-1)
    assert abs(sum(f0.values()) - 1.0) < 1e-12 and 'NiCit_aq4' in f0
    assert f1['NiCit_aq4'] > f0['NiCit_aq4'], (f0['NiCit_aq4'], f1['NiCit_aq4'])

    # 시트르산 결과가 없으면(extra=None) 기존 동작과 완전히 같아야 한다
    fa, ta = fractions(th, 1e-3, 1.0, 333.15, 0.0)
    fb, tb = fractions(th, 1e-3, 1.0, 333.15, 0.0, 1.0, None, 1e-2)
    assert fa == fb and ta == tb

    # CA 활동도: x(ChCl) 이 커지면 시트르산 활동도는 줄어야 한다
    ca = [a_ca(grid, x, 333.15) for x in (0.05, 0.3, 0.6, 0.9)]
    assert all(ca[i] > ca[i + 1] for i in range(len(ca) - 1)), ca

    # 종이 일부 빠져도 죽지 않는가 (freq 실패 대비)
    part = {k: v for k, v in th.items() if k != 'NiCl3_aq1'}
    fp, _ = fractions(part, 1e-2, 1.0, 333.15, 0.0)
    assert abs(sum(fp.values()) - 1.0) < 1e-12 and len(fp) == 4

    # dS 가 얕으면 온도전환이 완만해 두 앵커를 동시에 못 맞춘다 -> alpha=None 로 보고해야 함
    flat = {k: (h, sv * 0.4) for k, (h, sv) in th.items()}
    _, a_flat = calibrate(flat, grid, 1.0)
    assert a_flat is None, '얕은 dS 인데 보정이 성공해버림: %s' % a_flat

    print('self-check OK   shift=%+.1f kJ/mol per Cl-,  alpha=%.3f' % (shift, alpha))
    print('  alpha<1 = 이온쌍 대리변수가 자유 Cl- 활동도의 조성의존성을 과대평가함을 보정')


def main():
    if '--test' in sys.argv:
        test(); return
    arg = lambda k, d: (float(sys.argv[sys.argv.index(k) + 1]) if k in sys.argv else d)
    aW, ftd_max = arg('--awater', 1.0), arg('--ftdmax', FTD_MAX)

    th, missing, imag = harvest()
    extra, imag2 = harvest_extra()      # citrate/ 없으면 {} 라서 이후 경로가 그대로 동작
    imag = imag + imag2
    if imag:
        print('⚠ 허수진동수 발견 — dG 신뢰도 저하:')
        for line in imag:
            print('   ', line)
        print('   50 cm-1 이하면 무시 가능(물 회전). 크면 그 종을 재최적화할 것.\n')
    if th is None:
        print('ORCA 결과 부족: %s' % ', '.join(missing))
        print('freq 가 끝난 뒤 다시 실행하세요. 로직 점검: python speciation.py --test')
        return
    if missing:
        print('빠진 종 (그 종 없이 계산): %s\n' % ', '.join(missing))

    print('반응 열역학 (Ni_aq6 기준, 보정 전)')
    for name, n, m, geom in LADDER:
        if name in th:
            dH, dS = th[name]
            print('  %-11s n(Cl)=%d  %-3s  dH=%+8.1f kJ/mol  dS=%+7.1f J/mol/K  dG(60C)=%+8.1f'
                  % (name, n, geom, dH, dS * 1000, dH - 333.15 * dS))

    globals()['QH'] = '--raw' not in sys.argv
    use_water = '--water' in sys.argv
    grid = load_activity(use_water)
    print('\n활동도 모델: %s'
          % ('3성분 ChCl:CA:H2O (x_H2O=0.80, 실험 레시피 그대로)' if use_water
             else '2성분 ChCl:CA (물 없음) — 실험계와 다름. --water 권장'))
    shift, alpha = calibrate(th, grid, aW, ftd_max, extra)
    if shift is None:
        print('\n앵커1 실패: 이 dG 로는 95 C 에서 f_Td=0.5 를 만들 수 없습니다.')
        print('  dG 부호/단위, 허수진동수 오염을 확인하세요.')
        return
    if alpha is None:
        alpha = 1.0
        print('\n⚠ 앵커2(UV-Vis 상한)를 맞추는 alpha 가 없습니다. alpha=1 로 진행합니다.')
        print('  염 과잉 조성이 포화될 수 있으니 지도를 정량으로 읽지 마세요.')

    print('\n보정 (실측 앵커 2점)')
    print('  앵커1 Hartley 2025 : %.0f C, x=%.3f 에서 f_Td=%.2f'
          % (ANCHOR1_T - 273.15, ANCHOR1_X, ANCHOR1_F))
    print('  앵커2 우리 UV-Vis  : %.0f C, x=%.2f 에서 f_Td=%.2f (미검출 환산, a_W=%.2f)'
          % (ANCHOR2_T - 273.15, ANCHOR2_X, ftd_max, aW))
    print('  -> shift = %+.1f kJ/mol per Cl-,  alpha = %.4f' % (shift, alpha))
    if alpha < 0.05:
        print('  ** 경고: alpha 가 사실상 0 입니다. 두 앵커를 동시에 맞추려면')
        print('     조성 의존성이 없어야 한다는 뜻 -> 지도가 온도축으로만 변합니다.')
        print('     이 경우 가이드라인은 "조성을 고르라"가 아니라 "온도를 넘겨라"가 됩니다.')
        print('     조성비를 주장하려면 자유 Cl- 활동도를 제대로 내는 모델이 필요합니다.')
    print('     |shift| > 50 이면 dG 자체를 의심할 것 (이온 용매화 오차 통상 범위)')

    if extra:
        print()
        print('시트르산 착물 (전하중립 반응 -> shift 보정 없음, Cl 사다리보다 신뢰도 높음)')
        for name, m, geom, ncit in EXTRA:
            if name in extra:
                dH, dS = extra[name]
                print('  %-11s %-3s  dH=%+8.1f kJ/mol  dS=%+7.1f J/mol/K  dG(60C)=%+8.1f'
                      % (name, geom, dH, dS * 1000, dH - 333.15 * dS))
    else:
        print('\n시트르산 결과 없음 -> Cl/H2O 만으로 계산 (citrate/ 확인)')

    print('\n60 C 측정 5조성의 예상 화학종 분율')
    hdr = [n for n, _, _, _ in LADDER] + [n for n, _, _, _ in EXTRA if n in (extra or {})]
    print('  %-5s %-6s %-10s %s | f_Td'
          % ('시료', 'x', 'a_ChCl', ' '.join('%-10s' % h for h in hdr)))
    for x, lab in MEAS:
        a = a_cl(grid, x, 333.15)
        f, t = fractions(th, a, aW, 333.15, shift, alpha, extra, a_ca(grid, x, 333.15))
        print('  %-5s %-6.2f %-10.2e %s | %.4f'
              % (lab, x, a, ' '.join('%-10.4f' % f.get(h, 0.0) for h in hdr), t))

    write_outputs(th, shift, alpha, grid, aW, extra)
    print('\n주의: a_Cl 은 이온쌍 근사 대리변수, a_W=%.2f 는 가정값. 절대값 인용 금지.' % aW)
    print('      이 지도는 앵커 2점으로 보정된 유효모델이다 (검증 아님, 보정임).')


if __name__ == '__main__':
    main()
