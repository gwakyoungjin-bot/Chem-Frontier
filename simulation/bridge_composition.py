"""조성비 -> COSMO-RS 활동도 -> 예상 Ni 화학종 -> 예상 d-d 밴드.
실측 UV-Vis 5조성과 DFT 계산을 잇는 다리.

    python bridge_composition.py          # 실험온도(60 C) 기준 표
    python bridge_composition.py --test    # 보간 자체검증

주의: features_all.csv의 온도 격자는 25/45/65/85/100 C 뿐이라 60 C는 격자점이 없다.
      45 C와 65 C를 선형보간해서 쓴다 (w = 0.75).
"""
import csv, math, sys

T_EXP = 333.15                                  # 실험 온도 60 C (2026-09-23 확정)
MEAS = {'1:19': 0.05, '3:17': 0.15, '1:4': 0.20, '3:1': 0.75, '19:1': 0.95}
# 활동도 구간 -> 우세 화학종 (DFT로 계산하는 5종 중)
GUESS = [(1e-3, 'Ni_aq6  [Ni(H2O)6]2+   팔면체·약한 밴드 (395/725/1176nm)'),
         (1e-1, 'NiCl2_aq4 부근         팔면체 유지, 약간 적색이동'),
         (9e9,  'NiCl4 / NiCl3_aq1      사면체·강한 밴드 (가시영역)')]


def load(path='feature_extraction/features_all.csv'):
    rows = [r for r in csv.DictReader(open(path))
            if r['hba_name'] == 'choline_chloride' and r['hbd_name'] == 'citric_acid']
    g = {}
    for r in rows:
        g.setdefault(float(r['T']), {})[round(float(r['x_hba']), 2)] = r
    return g


def interp(g, T, x, field):
    """온도 격자 사이를 선형보간. 격자점이면 그대로."""
    Ts = sorted(g)
    if T in g:
        return float(g[T][x][field])
    lo = max(t for t in Ts if t <= T)
    hi = min(t for t in Ts if t >= T)
    w = (T - lo) / (hi - lo)
    return (1 - w) * float(g[lo][x][field]) + w * float(g[hi][x][field])


def test():
    g = load()
    assert sorted(g) == [298.15, 318.15, 338.15, 358.15, 373.15], sorted(g)
    x = 0.20
    lo, hi = interp(g, 338.15, x, 'lng_hba_tot'), interp(g, 358.15, x, 'lng_hba_tot')
    mid = interp(g, 348.15, x, 'lng_hba_tot')
    assert min(lo, hi) <= mid <= max(lo, hi), (lo, mid, hi)      # 보간값은 사이에
    assert abs(mid - (lo + hi) / 2) < 1e-9, mid                  # 중점이면 평균
    assert abs(interp(g, 338.15, x, 'lng_hba_tot') - lo) < 1e-12  # 격자점은 그대로
    assert all(x in g[338.15] for x in MEAS.values()), '측정 조성이 격자에 없음'
    print('self-check OK')


def main():
    g = load()
    f = lambda x, fld: interp(g, T_EXP, x, fld)
    act = lambda x: x * math.exp(f(x, 'lng_hba_tot'))
    xs = sorted(g[338.15])
    xopt = min(xs, key=lambda x: f(x, 'gE_over_RT'))

    print(f"ChCl : citric acid,  {T_EXP - 273.15:.0f} C (실험온도, 45/65 C 보간)\n")
    print(f"{'시료':>6}{'x_ChCl':>8}{'gE/RT':>9}{'ChCl 활동도':>14}{'최적점과 거리':>14}")
    for k, x in sorted(MEAS.items(), key=lambda kv: kv[1]):
        print(f"{k:>6}{x:8.2f}{f(x,'gE_over_RT'):9.3f}{act(x):14.2e}{abs(x-xopt):14.2f}")

    print(f"\n계산상 비이상성 최대 조성  x_ChCl = {xopt:.2f}  (gE/RT = {f(xopt,'gE_over_RT'):.3f})")
    print(f"  = ChCl:CA 몰비 약 {xopt/(1-xopt):.2f} : 1   <- 측정 5점 중 아무도 여기 없음")
    print(f"  문헌 관행 1:2 (x=0.333) 근처 격자점 x=0.35 의 gE/RT = {f(0.35,'gE_over_RT'):.3f}")

    lo, hi = min(act(x) for x in MEAS.values()), max(act(x) for x in MEAS.values())
    print(f"\n5개 시료의 ChCl 활동도 범위: {lo:.1e} ~ {hi:.1e}  ({math.log10(hi/lo):.0f}자릿수 차이)")
    print("-> 염화물 활동도가 자릿수로 갈리므로 Ni 화학종이 바뀌어야 정상.\n")

    print(f"{'시료':>6}{'활동도':>11}  예상 우세 화학종 / 예상 d-d 거동")
    for k, x in sorted(MEAS.items(), key=lambda kv: kv[1]):
        a = act(x)
        print(f"{k:>6}{a:11.1e}  {next(s for t, s in GUESS if a < t)}")

    print("\n[한계] lng는 ChCl을 '중성 이온쌍 1성분'으로 다룬 값(HANDOFF C-4)이라")
    print("       자유 Cl- 활동도 자체가 아니라 그 대리변수다. 순서와 자릿수 차이는 유효.")
    print("[한계] 60 C는 온도격자(25/45/65/85/100 C)에 없어 선형보간. 45/65 C 각각에서도")
    print("       gE 최소점은 x=0.40~0.45 구간으로 동일하므로 결론은 보간에 의존하지 않는다.")


if __name__ == '__main__':
    test() if '--test' in sys.argv else main()
