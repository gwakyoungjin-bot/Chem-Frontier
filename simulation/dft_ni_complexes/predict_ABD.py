"""인과 검증용 예측 3종 — 실험 전에 값을 먼저 박아두기 위한 계산.

배경:
  FT-IR 의 용매 구조 전이와 우리가 계산한 금속 화학종 전이가 같은 조성에서
  일어나지만, 인과를 설명하려던 두 시도가 모두 실패했다
  (수소결합 화학량론 -> 물이 지배, COSMO-RS enth 분해 -> 이상치 없음).
  남은 길은 **예측을 먼저 내고 실험으로 맞추는 것**이다.

전제 (이게 참이면 아래 예측이 성립한다):
  금속 화학종 전이는 **염화물 활동도가 임계값 a* 에 도달할 때** 일어난다.
  현재 계(x_H2O = 0.80)에서 전이는 x = 0.45~0.50 이므로 그 구간의 a_ChCl 이 a* 다.
  a* 가 물질·물함량에 무관한 상수라면:

  A  물 함량을 바꾸면 전이 조성이 예측대로 이동해야 한다.
  B  HBD 를 바꾸면(말론산·타르타르산) 전이 조성이 예측대로 달라져야 한다.

  둘 다 틀리면 '염화물 활동도가 전이를 구동한다'는 우리 전제 자체가 틀린 것이다.
  이건 사후 설명이 아니라 **반증 가능한 예측**이다.

실행: python dft_ni_complexes/predict_ABD.py
출력: dft_results/predict_ABD.csv
"""
import csv
import math
import os
import sys

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'literature', 'offline_bundle',
                                'openCOSMO-RS_py', 'src'))
import numpy as np                                    # noqa: E402
from opencosmorspy import COSMORS                     # noqa: E402

B = os.path.join(ROOT, 'literature', 'offline_bundle', 'openCOSMO-RS_py',
                 'tests', 'COSMO_ORCA')
FE = os.path.join(ROOT, 'feature_extraction')
CHCL = os.path.join(FE, 'orcacosmo', 'choline_chloride', 'choline_chloride_c000.orcacosmo')
WATER = os.path.join(B, 'H2O', 'COSMO_TZVPD', 'H2O_c000.orcacosmo')
HBD_DIR = os.path.join(FE, 'handover_20260908', 'orcacosmo')
HBDS = {'citric_acid': ('시트르산', 3), 'malonic_acid': ('말론산', 2),
        'tartaric_acid': ('타르타르산', 2)}
T = 333.15
X_REF, XW_REF = 0.475, 0.80        # 현재 계의 전이 조성(0.45~0.50 중앙)과 물함량


def hbd_path(name):
    return os.path.join(HBD_DIR, name, '%s_c000.orcacosmo' % name)


def activities(hbd, xw, xs):
    """조성 목록에 대한 a_ChCl. a = x_ChCl(전체 기준) * gamma."""
    crs = COSMORS(par='default_orca')
    crs.add_molecule([CHCL]); crs.add_molecule([hbd]); crs.add_molecule([WATER])
    for x in xs:
        crs.add_job(np.array([x * (1 - xw), (1 - x) * (1 - xw), xw]), T,
                    refst='pure_component')
    r = crs.calculate()
    out = []
    for i, x in enumerate(xs):
        xc = x * (1 - xw)
        out.append(xc * math.exp(float(r['tot']['lng'][i][0])))
    return out


def solve_x(hbd, xw, target, lo=0.02, hi=0.98):
    """a_ChCl = target 이 되는 조성을 이분법으로. 단조 증가를 가정한다."""
    for _ in range(18):
        mid = (lo + hi) / 2
        a = activities(hbd, xw, [mid])[0]
        if a < target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def main():
    cit = hbd_path('citric_acid')
    if not os.path.exists(cit):
        sys.exit('시트르산 프로파일 없음')

    a_star = activities(cit, XW_REF, [X_REF])[0]
    print(' 기준: 현재 계(시트르산, 물 %.2f)의 전이 조성 x = %.3f' % (XW_REF, X_REF))
    print(' 임계 염화물 활동도  a* = %.3e' % a_star)
    print(' -> 아래 예측은 전부 "a_ChCl 이 a* 에 닿는 조성" 이다.\n')

    rows = []

    print(' [A] 물 함량을 바꾸면 전이 조성이 어디로 가는가')
    print(' %-10s %-12s %s' % ('물 몰분율', '예측 전이 x', '현재 대비'))
    print(' ' + '-' * 42)
    for xw in (0.90, 0.80, 0.70, 0.60, 0.50, 0.30):
        xs = solve_x(cit, xw, a_star)
        d = xs - X_REF
        print(' %-10.2f %-12.3f %+.3f' % (xw, xs, d))
        rows.append(dict(test='A', hbd='citric_acid', x_H2O='%.2f' % xw,
                         x_transition='%.3f' % xs))

    print('\n [B] HBD 를 바꾸면 전이 조성이 어디로 가는가  (물 %.2f 고정)' % XW_REF)
    print(' %-14s %-8s %-12s %s' % ('HBD', 'COOH수', '예측 전이 x', '비고'))
    print(' ' + '-' * 52)
    for name, (ko, ncooh) in HBDS.items():
        p = hbd_path(name)
        if not os.path.exists(p):
            print(' %-14s (프로파일 없음)' % ko)
            continue
        try:
            xs = solve_x(p, XW_REF, a_star)
        except Exception as e:
            print(' %-14s 실패: %s' % (ko, e))
            continue
        note = '기준' if name == 'citric_acid' else ''
        print(' %-14s %-8d %-12.3f %s' % (ko, ncooh, xs, note))
        rows.append(dict(test='B', hbd=name, x_H2O='%.2f' % XW_REF,
                         x_transition='%.3f' % xs))

    out = os.path.join(ROOT, 'dft_results', 'predict_ABD.csv')
    with open(out, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['test', 'hbd', 'x_H2O', 'x_transition'])
        w.writeheader(); w.writerows(rows)
    print('\n wrote %s' % out)

    print('\n [해석]')
    print('  이 값들은 **실험 전에 계산된 예측**이다. 맞으면 염화물 활동도가')
    print('  전이를 구동한다는 전제가 지지되고, 틀리면 그 전제를 버려야 한다.')
    print('  어느 쪽이든 사후 설명이 아니라 검증이다.')

    print('\n ⚠ 한계')
    print('  1) a* 를 현재 계에서 뽑았으므로 A 의 x_H2O=0.80 행은 자명하게 맞는다.')
    print('     검증력은 다른 물함량과 다른 HBD 에만 있다.')
    print('  2) ChCl 을 중성 이온쌍으로 근사한 한계가 그대로 따라온다.')
    print('  3) 말론산·타르타르산 DES 가 그 조성에서 액체인지는 확인하지 않았다.')

    assert rows, '예측이 하나도 안 나왔다'
    print('\n self-check OK: 예측 %d건' % len(rows))


if __name__ == '__main__':
    main()
