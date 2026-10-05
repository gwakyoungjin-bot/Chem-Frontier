"""HBD 를 바꾸면 염화물 활동도가 얼마나 달라지는가.

왜 필요한가 (2026-09-27):
  앵커 문헌은 대부분 릴라인(ChCl:urea) 이나 에타린(ChCl:EG) 로 되어 있는데
  우리 계는 ChCl:citric acid 다. 앵커를 쓰려면 "HBD 가 달라도 a_ChCl 은 비슷하다" 를
  **가정**해야 하는데, 그 가정을 숫자로 확인하지 않고 쓰는 건 위험하다.

  이 스크립트가 같은 조성·같은 물함량에서 HBD 만 바꿔 a_ChCl 을 비교한다.
    차이가 작으면 -> 앵커 전용이 정당화됨 (한계 절에 이 수치를 근거로 기재)
    차이가 크면   -> HBD 가 맞는 앵커를 더 찾아야 함

  새 DFT 불필요: 네 분자 모두 .orcacosmo 가 이미 있다.

비교 대상 물함량:
  0.045  Hartley 2025 (ChCl:EG, NiCl2.6H2O 0.1 M) — 현재 Ni 앵커1
  0.80   우리 실험 (DES 3000 mg + 물 1270 uL)
  0.00   건조 기준 (참고)

실행:
    python hbd_effect.py
    python hbd_effect.py --test
"""
import math, os, sys

PIPE = os.path.expanduser(
    '~/software/cosmors_offline/openCOSMO-RS_conformer_pipeline')
HERE = os.path.dirname(os.path.abspath(__file__))

MOL = {
    'ChCl': 'choline_chloride/COSMO_TZVPD/choline_chloride_c000.orcacosmo',
    'CA':   'citric_acid/COSMO_TZVPD/citric_acid_c000.orcacosmo',
    'EG':   'ethylene_glycol/COSMO_TZVPD/ethylene_glycol_c000.orcacosmo',
    'urea': 'urea/COSMO_TZVPD/urea_c000.orcacosmo',
    'H2O':  'water/COSMO_TZVPD/water_c000.orcacosmo',
}

# HBD -> (표시명, ChCl:HBD 관행 몰비에서의 x_ChCl(건조기준))
# 셋 다 1:2 로 맞춰 비교한다. 관행 몰비가 1:2 인 계(릴라인·에타린)와 동일 조건.
HBDS = ['CA', 'EG', 'urea']
X_CHCL_DRY = 1.0 / 3.0
WATERS = [0.00, 0.045, 0.40, 0.80]
T = 368.15          # Hartley 전환온도 구간 중앙 (95 C)


def path(k):
    p = os.path.join(PIPE, MOL[k])
    if not os.path.exists(p):
        sys.exit('없음: %s' % p)
    return p


def a_chcl(hbd, x_w, T=T):
    """ChCl:HBD(1:2) + 물 혼합물에서 ChCl 활동도."""
    import numpy as np
    from opencosmorspy import COSMORS
    crs = COSMORS(par='default_orca')
    for k in ('ChCl', hbd, 'H2O'):
        crs.add_molecule([path(k)])
    xd = 1.0 - x_w
    x = np.array([X_CHCL_DRY * xd, (1 - X_CHCL_DRY) * xd, x_w])
    crs.add_job(x, T, refst='pure_component')
    lng = np.array(crs.calculate()['tot']['lng'])[0]
    return float(x[0]) * math.exp(float(lng[0]))


def test():
    # 물이 늘면 어느 HBD 든 a_ChCl 이 줄어야 한다
    v = [a_chcl('EG', w) for w in (0.0, 0.5, 0.9)]
    assert all(v[i] > v[i + 1] for i in range(len(v) - 1)), v
    assert all(0 < x < 1 for x in v), v
    print('self-check OK')


def main():
    if '--test' in sys.argv:
        test(); return
    print('a_ChCl  (ChCl:HBD = 1:2 건조기준, T = %.0f C)' % (T - 273.15))
    print()
    print('%-8s %-12s %-12s %-12s %-14s' %
          ('x(H2O)', 'CA(우리)', 'EG(Hartley)', 'urea(릴라인)', 'HBD 최대/최소'))
    rows = []
    for w in WATERS:
        vals = {h: a_chcl(h, w) for h in HBDS}
        hi, lo = max(vals.values()), min(vals.values())
        print('%-8.3f %-12.3e %-12.3e %-12.3e %-14.1f배'
              % (w, vals['CA'], vals['EG'], vals['urea'], hi / lo))
        rows.append((w, vals, hi / lo))

    print()
    print('[비교] 같은 HBD 에서 물함량만 바꿀 때의 변화폭')
    for h in HBDS:
        v0 = a_chcl(h, 0.045)      # Hartley 조건
        v1 = a_chcl(h, 0.80)       # 우리 조건
        print('  %-5s : x_w 0.045 -> 0.80 에서 %.0f배 감소' % (h, v0 / v1))

    print()
    print('[판정 기준]')
    print('  물함량 효과 >> HBD 효과  이면 -> 물만 맞추면 HBD 는 달라도 된다')
    print('  두 효과가 비슷하면      -> HBD 도 맞는 앵커를 찾아야 한다')


if __name__ == '__main__':
    main()
