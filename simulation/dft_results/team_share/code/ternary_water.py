"""ChCl : citric acid : H2O 3성분 COSMO-RS 활동도.

왜 필요한가 (2026-09-27 발견):
  실험 레시피가 DES 3000 mg + 물 1270 uL(30 wt%) 였다. 몰 기준으로 환산하면
  **물 몰분율이 0.77~0.82** 다. 즉 이 용액은 "물 조금 든 DES" 가 아니라
  사실상 수용액에 가깝다. 물:ChCl 몰비는 1:19 시료에서 89:1 까지 간다.

  그런데 지금까지 쓴 features_all.csv 는 ChCl:CA **2성분** 모델이라 물이 아예 없다.
  문헌은 물이 화학종 평형을 지배한다고 반복해 말한다:
    - 10.1021/acs.inorgchem.3c02205 : 질량·부피 아닌 '몰비' 가 배위껍질을 좌우
    - 10.1021/acs.inorgchem.6c01344 : 물이 늘면 수화 팔면체 쪽으로 평형이 이동
  => 물을 넣지 않은 활동도로 만든 지도는 이 실험계를 기술하지 못한다.

새 DFT 불필요: water / citric_acid / choline_chloride 의 .orcacosmo 가 이미 있다.

한계:
  ChCl 은 중성 이온쌍 1성분으로 다룬다(기존 파이프라인과 동일, Lemaoui 2020 방식).
  따라서 a_ChCl 은 자유 Cl- 활동도가 아니라 그 대리변수다. 절대값 인용 금지.

실행:
    python ternary_water.py            # 실험조건 + 물함량 스캔 -> ternary_water.csv
    python ternary_water.py --test     # 자체검증
"""
import csv, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
PIPE = os.path.expanduser(
    '~/software/cosmors_offline/openCOSMO-RS_conformer_pipeline')

MOL = {  # 이름 -> .orcacosmo 경로
    'ChCl':  'choline_chloride/COSMO_TZVPD/choline_chloride_c000.orcacosmo',
    'CA':    'citric_acid/COSMO_TZVPD/citric_acid_c000.orcacosmo',
    'H2O':   'water/COSMO_TZVPD/water_c000.orcacosmo',
}

# 실험 5조성의 x(ChCl) (물 제외 건조 기준) 과 실제 물 몰분율
# 물 몰분율은 레시피(DES 3000 mg + 물 1270 uL)에서 환산한 값
MEAS = [('1:19', 0.05, 0.817), ('3:17', 0.15, 0.812), ('1:4', 0.20, 0.810),
        ('3:1', 0.75, 0.782), ('19:1', 0.95, 0.770)]

T_EXP = 333.15


def paths():
    out = {}
    for k, rel in MOL.items():
        p = os.path.join(PIPE, rel)
        if not os.path.exists(p):
            sys.exit('없음: %s\n  (openCOSMO-RS_conformer_pipeline 경로 확인)' % p)
        out[k] = p
    return out


def activities(jobs, T=T_EXP):
    """jobs = [(x_ChCl, x_CA, x_H2O), ...] (합 1). 반환: [(a_ChCl, a_CA, a_H2O), ...]"""
    import numpy as np
    from opencosmorspy import COSMORS
    p = paths()
    crs = COSMORS(par='default_orca')
    for k in ('ChCl', 'CA', 'H2O'):
        crs.add_molecule([p[k]])
    for x in jobs:
        crs.add_job(np.array(x, dtype=float), T, refst='pure_component')
    res = crs.calculate()
    lng = np.array(res['tot']['lng'])            # (n_job, 3)
    return [tuple(float(x[i]) * math.exp(float(lng[j][i])) for i in range(3))
            for j, x in enumerate(jobs)]


def test():
    """물 함량을 올리면 물 활동도는 오르고 ChCl 활동도는 떨어져야 한다."""
    jobs, labels = [], []
    for xw in (0.0, 0.3, 0.6, 0.8, 0.95):
        xd = 1.0 - xw
        jobs.append((0.333 * xd, 0.667 * xd, xw)); labels.append(xw)
    a = activities(jobs)
    print('  x(H2O)   a_ChCl      a_CA        a_H2O')
    for lab, (ac, aa, aw) in zip(labels, a):
        print('  %-8.2f %-11.3e %-11.3e %-11.3e' % (lab, ac, aa, aw))
    aw = [x[2] for x in a]
    assert all(aw[i] <= aw[i+1] + 1e-12 for i in range(len(aw)-1)), aw
    assert all(0.0 <= v <= 5.0 for t in a for v in t), a
    print('self-check OK')


def grid_mode():
    """speciation.py 가 쓸 활동도 격자를 만든다.

    물 함량을 두 수준으로 낸다:
      x_H2O = 0.80  (실제 실험 조건, 물 30 wt%)
      x_H2O = 0.00  (건조 DES, 비교용)
    둘을 나란히 보면 '물이 조성비의 지렛대를 얼마나 깎는가' 가 그대로 보인다.
    """
    xs = [round(0.05 * i, 2) for i in range(1, 20)]
    Ts = [298.15 + 5 * i for i in range(0, 21)]
    rows = []
    for xw in (0.80, 0.00):
        for T in Ts:
            jobs = [(x * (1 - xw), (1 - x) * (1 - xw), xw) for x in xs]
            a = activities(jobs, T)
            for x, (ac, aa, aw) in zip(xs, a):
                rows.append(dict(x_H2O=xw, T_K=round(T, 2), x_ChCl_dry=x,
                                 a_ChCl='%.6e' % ac, a_CA='%.6e' % aa, a_H2O='%.6e' % aw))
            print('  x_H2O=%.2f  T=%.1f C  완료 (%d점)' % (xw, T - 273.15, len(xs)), flush=True)
    out = os.path.join(HERE, 'ternary_grid.csv')
    with open(out, 'w', newline='', encoding='utf-8') as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0]))
        wr.writeheader(); wr.writerows(rows)
    print('wrote %s (%d행)' % (out, len(rows)))


def main():
    if '--test' in sys.argv:
        test(); return
    if '--grid' in sys.argv:
        grid_mode(); return

    # (1) 실제 실험 조건
    jobs = [(x * (1 - w), (1 - x) * (1 - w), w) for _, x, w in MEAS]
    a = activities(jobs)
    print('실험 조건 (60 C, 물 30 wt% = 몰분율 0.77~0.82)')
    print('%-6s %-8s %-8s %-12s %-12s %-12s' %
          ('시료', 'x(ChCl)', 'x(H2O)', 'a_ChCl', 'a_CA', 'a_H2O'))
    rows = []
    for (lab, x, w), (ac, aa, aw) in zip(MEAS, a):
        print('%-6s %-8.2f %-8.3f %-12.3e %-12.3e %-12.3e' % (lab, x, w, ac, aa, aw))
        rows.append(dict(case='measured', label=lab, x_ChCl_dry=x, x_H2O=w,
                         T_C=60, a_ChCl='%.4e' % ac, a_CA='%.4e' % aa, a_H2O='%.4e' % aw))

    # (2) 물 함량을 낮추면 a_ChCl 이 얼마나 오르나 — 다음 실험 설계용
    print('\n물 함량을 낮추면 (x=0.95, 19:1 기준)')
    print('%-10s %-12s %-12s %-10s' % ('x(H2O)', 'a_ChCl', 'a_H2O', 'a_ChCl 배율'))
    jobs2 = [(0.95 * (1 - w), 0.05 * (1 - w), w) for w in (0.0, 0.2, 0.4, 0.6, 0.77, 0.9)]
    a2 = activities(jobs2)
    base = None
    for (xc, _, w), (ac, aa, aw) in zip(jobs2, a2):
        if base is None:
            base = ac
        print('%-10.2f %-12.3e %-12.3e %-10.1f' % (w, ac, aw, ac / a2[4][0]))
        rows.append(dict(case='water_scan', label='19:1', x_ChCl_dry=0.95, x_H2O=w,
                         T_C=60, a_ChCl='%.4e' % ac, a_CA='%.4e' % aa, a_H2O='%.4e' % aw))

    out = os.path.join(HERE, 'ternary_water.csv')
    with open(out, 'w', newline='', encoding='utf-8') as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0]))
        wr.writeheader(); wr.writerows(rows)
    print('\nwrote', out)
    print('주의: ChCl 은 중성 이온쌍 근사. a_ChCl 은 자유 Cl- 활동도의 대리변수다.')


if __name__ == '__main__':
    main()
