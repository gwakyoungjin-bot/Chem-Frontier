"""염화물 활동도를 상호작용/조합 성분으로 분해해 FT-IR 전이와 대조한다.

왜 이걸 하나:
  a_ChCl 전체는 조성에 대해 매끄럽다(log 선형 R^2 = 0.978). 그래서 앞서
  "FT-IR 의 꺾임과 모양이 다르다 -> 다른 현상" 이라고 판정했다.
  그런데 lnγ 는 두 성분의 합이다:
     comb  조합 항 — 분자 크기·모양에서 오는 것. 구조상 **매끄러울 수밖에 없다**.
     enth  상호작용 항 — 수소결합·정전기. **FT-IR 이 보는 것이 바로 이쪽이다.**
  전체가 매끄러운 게 comb 이 지배해서라면, enth 만 떼어냈을 때 꺾일 수 있다.
  그러면 "다른 현상" 이 아니라 "같은 상호작용의 두 관측" 이 된다.

판정:
  enth 성분의 국소 기울기에 FT-IR 과 같은 위치(x ~ 0.5)에서 이상치가 나오는가.
  나오면 인과 연결의 근거가 되고, 안 나오면 앞선 판정(다른 현상)이 유지된다.

실행: python dft_ni_complexes/cosmo_decompose.py
"""
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
MOL = [
    os.path.join(FE, 'orcacosmo', 'choline_chloride', 'choline_chloride_c000.orcacosmo'),
    os.path.join(FE, 'handover_20260908', 'orcacosmo', 'citric_acid', 'citric_acid_c000.orcacosmo'),
    os.path.join(B, 'H2O', 'COSMO_TZVPD', 'H2O_c000.orcacosmo'),
]
X_H2O, T = 0.80, 333.15
XS = [0.05, 0.10, 0.15, 0.20, 0.25, 0.333, 0.40, 0.45, 0.50,
      0.55, 0.60, 0.667, 0.75, 0.85, 0.95]
# FT-IR 실측 C=O (ftir_trend.py). 비교용
CO = {0.05: 1713.9, 0.10: 1714.0, 0.15: 1714.4, 0.20: 1714.7, 0.25: 1714.8,
      0.333: 1715.3, 0.40: 1715.6, 0.50: 1716.0, 0.60: 1722.3, 0.667: 1722.9,
      0.75: 1723.5, 0.85: 1726.8, 0.95: 1729.2}


def main():
    for p in MOL:
        if not os.path.exists(p):
            sys.exit('프로파일 없음: %s' % p)

    crs = COSMORS(par='default_orca')
    for p in MOL:
        crs.add_molecule([p])
    for x in XS:
        crs.add_job(np.array([x * (1 - X_H2O), (1 - x) * (1 - X_H2O), X_H2O]),
                    T, refst='pure_component')
    r = crs.calculate()

    print(' ChCl 의 ln(활동도계수) 분해   (%.2f K, x_H2O = %.2f)' % (T, X_H2O))
    print(' %-8s %-11s %-11s %-11s %s'
          % ('x(ChCl)', 'ln γ 전체', 'enth(상호작용)', 'comb(조합)', 'C=O 실측'))
    print(' ' + '-' * 62)
    rows = []
    for i, x in enumerate(XS):
        tot = float(r['tot']['lng'][i][0])
        en = float(r['enth']['lng'][i][0])
        cb = float(r['comb']['lng'][i][0])
        rows.append((x, tot, en, cb))
        co = CO.get(x)
        print(' %-8.3f %-11.3f %-11.3f %-11.3f %s'
              % (x, tot, en, cb, ('%.1f' % co) if co else '—'))

    def slopes(vals):
        out = []
        for (x1, v1), (x2, v2) in zip(zip(XS, vals), list(zip(XS, vals))[1:]):
            out.append((x1, x2, (v2 - v1) / (x2 - x1)))
        return out

    print('\n [국소 기울기] 이상치가 어디에 있나  (중앙값 대비 배수)')
    for name, idx in (('전체 lnγ', 1), ('enth 상호작용', 2), ('comb 조합', 3)):
        sl = slopes([row[idx] for row in rows])
        mags = sorted(abs(s) for _, _, s in sl)
        med = mags[len(mags) // 2]
        worst = max(sl, key=lambda t: abs(t[2]))
        ratio = abs(worst[2]) / med if med else float('inf')
        print('  %-14s 최대 이상치 x=%.2f~%.2f  (중앙값의 %.1f배)'
              % (name, worst[0], worst[1], ratio))

    # FT-IR 과의 직접 대조: enth 기울기와 C=O 기울기의 이상치 위치가 같은가
    xs_co = sorted(CO)
    co_sl = [((a + b) / 2, (CO[b] - CO[a]) / (b - a))
             for a, b in zip(xs_co, xs_co[1:])]
    co_med = sorted(abs(s) for _, s in co_sl)[len(co_sl) // 2]
    co_out = [c for c, s in co_sl if abs(s) > 4 * co_med]

    en_sl = slopes([row[2] for row in rows])
    en_med = sorted(abs(s) for _, _, s in en_sl)[len(en_sl) // 2]
    en_out = [(a + b) / 2 for a, b, s in en_sl if abs(s) > 4 * en_med]

    print('\n [대조] 이상치가 나타난 조성')
    print('  FT-IR C=O      : %s' % (['%.2f' % c for c in co_out] or '없음'))
    print('  enth 상호작용   : %s' % (['%.2f' % c for c in en_out] or '없음'))

    print('\n [판정]')
    if not en_out:
        print('  enth 성분에도 이상치가 없다. 상호작용 항까지 매끄럽다.')
        print('  => COSMO-RS 는 이 조성에서 어떤 구조 전이도 내놓지 않는다.')
        print('     FT-IR 의 꺾임과 연결할 계산적 근거가 없다 — 앞선 판정 유지.')
        print('     (단, 평균장 이론이라 협동적 네트워크 붕괴를 원리상 못 잡는다는')
        print('      한계 때문일 수 있다. 계산의 부재가 현상의 부재는 아니다.)')
    else:
        near = [c for c in en_out if any(abs(c - d) < 0.12 for d in co_out)]
        if near:
            print('  enth 이상치가 FT-IR 이상치와 %s 부근에서 겹친다.'
                  % ', '.join('%.2f' % c for c in near))
            print('  => 두 관측이 같은 상호작용을 보고 있다는 계산적 근거가 된다.')
        else:
            print('  enth 에 이상치는 있으나 FT-IR 과 위치가 다르다.')
            print('  => 인과 연결의 근거로 쓸 수 없다.')

    assert len(rows) == len(XS)
    print('\n self-check OK: %d조성 분해' % len(rows))


if __name__ == '__main__':
    main()
