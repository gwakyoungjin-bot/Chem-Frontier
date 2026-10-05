"""sigma-프로파일로 '유효 수소결합 공여자 / Cl- 수용자' 비를 계산하고 임계를 예측한다.

왜 이걸 하나:
  앞선 hbond_threshold.py 는 'COOH 3개' 라는 화학량론만 센 사후 설명이었다.
  n(Cl- 당 수소결합 수)을 2~4 로 열어두고 실측 구간에 맞춘 것이라 예측이 아니었다.
  여기서는 **분자마다 실제 수소결합 능력을 sigma-프로파일에서 직접 적분**해서
  임계를 먼저 계산한다. 콜린의 -OH 와 물도 자동으로 포함된다.

원리:
  COSMO-RS 에서 표면 세그먼트의 screening charge density sigma 는 분자 전하와 부호가 반대다.
    수소결합 **공여자**(산성 H, delta+)  -> sigma < -sigma_hb
    수소결합 **수용자**(비공유전자쌍)     -> sigma > +sigma_hb
  표준 임계 sigma_hb = 0.0084 e/A^2 를 쓴다.
  각 분자의 공여 면적과 수용 면적을 적분하면, 조성마다
    (전체 공여 면적) / (Cl- 수용 면적)
  을 낼 수 있다. 이 값이 1 로 떨어지는 조성이 '공여자가 Cl- 를 다 못 덮는' 임계다.

이 계산은 sigma 프로파일에서 나오므로 **가정한 n 이 없다**. 그게 앞 버전과의 차이다.

실행: python dft_ni_complexes/hbond_sigma.py
"""
import math
import os
import sys

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = os.path.join(ROOT, 'feature_extraction', 'handover_20260908', 'orcacosmo')
ALT = os.path.join(ROOT, 'literature', 'offline_bundle',
                   'openCOSMO-RS_py', 'tests', 'COSMO_ORCA')
# ChCl 과 Cl- 는 다른 수확 폴더에 있다
ALT2 = os.path.join(ROOT, 'feature_extraction', 'orcacosmo')

SIGMA_HB = 0.0084          # e/A^2, COSMO-RS 표준 수소결합 임계
BOHR2_TO_A2 = 0.529177210903 ** 2
X_H2O = 0.80               # 실험 물 몰분율

OBS = (0.40, 0.60)         # 실측 전이 구간 (육안·FT-IR·계산 전체 포괄)
KINK = 0.50                # FT-IR 꺾임


def find(name):
    for p in (os.path.join(BASE, name, '%s_c000.orcacosmo' % name),
              os.path.join(ALT2, name, '%s_c000.orcacosmo' % name),
              os.path.join(ALT, '%s.orcacosmo' % name),
              os.path.join(ALT, name, 'COSMO_TZVPD', '%s_c000.orcacosmo' % name)):
        if os.path.exists(p):
            return p
    return None


def profile(path):
    """(공여 면적, 수용 면적, 총면적) in A^2.

    표면점 블록의 area 와 COSMO_corrected 전하를 짝지어 sigma 를 만든다.
    보정 전하를 쓰는 이유: 그게 opencosmorspy 가 실제로 쓰는 값이다.
    """
    lines = open(path, encoding='utf-8', errors='replace').read().splitlines()
    i0 = next(i for i, l in enumerate(lines) if l.startswith('# SURFACE POINTS'))
    hdr = i0 + 2                                   # 컬럼명 줄
    areas = []
    i = hdr + 1
    while i < len(lines):
        p = lines[i].split()
        if len(p) < 10:
            break
        areas.append(float(p[3]))                  # area (a.u. = bohr^2)
        i += 1
    j0 = next((k for k, l in enumerate(lines) if l.startswith('#COSMO_corrected')), None)
    charges = []
    if j0 is not None:
        k = j0 + 1
        while k < len(lines):
            s = lines[k].strip()
            if s.startswith('#'):
                break
            try:
                charges.append(float(s))
            except ValueError:
                pass
            k += 1
    if len(charges) < len(areas):                  # 보정전하가 없으면 원 전하 사용
        charges = []
        i = hdr + 1
        while i < len(lines):
            p = lines[i].split()
            if len(p) < 10:
                break
            charges.append(float(p[5]))
            i += 1

    don = acc = tot = 0.0
    for a_bohr, q in zip(areas, charges):
        a = a_bohr * BOHR2_TO_A2
        sig = q / a if a > 0 else 0.0
        tot += a
        if sig < -SIGMA_HB:
            don += a
        elif sig > SIGMA_HB:
            acc += a
    return don, acc, tot


def main():
    need = {'citric_acid': 'CA', 'choline_chloride': 'ChCl 이온쌍',
            'water': '물', 'chloride_anion': 'Cl- 단독'}
    prof = {}
    print(' sigma-프로파일 적분  (sigma_hb = %.4f e/A^2)' % SIGMA_HB)
    print(' %-18s %-10s %-10s %-10s %s' % ('분자', '공여 A^2', '수용 A^2', '총 A^2', '비고'))
    print(' ' + '-' * 62)
    for name, ko in need.items():
        p = find(name)
        if not p:
            print(' %-18s (파일 없음)' % name)
            continue
        d, a, t = profile(p)
        prof[name] = (d, a, t)
        print(' %-18s %-10.1f %-10.1f %-10.1f %s' % (name, d, a, t, ko))

    if not all(k in prof for k in ('citric_acid', 'choline_chloride', 'water')):
        sys.exit('필요한 프로파일이 부족하다')

    d_ca = prof['citric_acid'][0]
    d_ch, a_ch = prof['choline_chloride'][0], prof['choline_chloride'][1]
    d_w = prof['water'][0]

    print('\n 조성별  (전체 공여 면적) / (Cl- 수용 면적)')
    print('  분자수 기준: x ChCl + (1-x) CA,  물은 몰분율 %.2f' % X_H2O)
    print(' %-9s %-12s %-12s %-12s %s' % ('x(ChCl)', '공여 총합', 'Cl- 수용', '비', ''))
    print(' ' + '-' * 56)
    rows = []
    for x in (0.05, 0.10, 0.15, 0.20, 0.25, 0.333, 0.40, 0.45, 0.50,
              0.55, 0.60, 0.667, 0.75, 0.85, 0.95):
        # 건조 기준 1몰에 대해 물 몰수
        n_w = X_H2O / (1 - X_H2O)
        don = x * d_ch + (1 - x) * d_ca + n_w * d_w
        acc = x * a_ch
        r = don / acc if acc > 0 else float('inf')
        rows.append((x, don, acc, r))
        mark = ''
        if abs(x - KINK) < 1e-9:
            mark = '  <== FT-IR 꺾임'
        elif x == 0.45:
            mark = '  <== 계산 f_Td 상승'
        print(' %-9.3f %-12.1f %-12.1f %-12.2f%s' % (x, don, acc, r, mark))

    print('\n [예측] 비가 1 이 되는 조성 (공여자가 Cl- 를 다 못 덮기 시작)')
    xstar = None
    for (x1, _, _, r1), (x2, _, _, r2) in zip(rows, rows[1:]):
        if (r1 - 1) * (r2 - 1) <= 0 and r1 != r2:
            xstar = x1 + (x2 - x1) * (r1 - 1) / (r1 - r2)
            break
    if xstar is None:
        print('  범위 안에서 1 을 지나지 않는다. 물의 공여가 압도적이라')
        print('  Cl- 는 어느 조성에서도 수소결합 상대가 남아돈다.')
        print('  => 이 지표로는 임계가 안 나온다. 아래 [재검증] 참고.')
    else:
        print('  x* = %.3f   실측 구간 %.2f~%.2f 안? %s'
              % (xstar, OBS[0], OBS[1], 'O' if OBS[0] <= xstar <= OBS[1] else 'X'))

    print('\n [재검증] 물을 빼고 유기 공여자(CA + 콜린)만 세면')
    rows2 = []
    for x in (0.20, 0.333, 0.40, 0.45, 0.50, 0.55, 0.60, 0.75, 0.95):
        don = x * d_ch + (1 - x) * d_ca
        acc = x * a_ch
        rows2.append((x, don / acc if acc else float('inf')))
    for x, r in rows2:
        mark = '  <== FT-IR 꺾임' if abs(x - KINK) < 1e-9 else ''
        print('  x=%.3f   비 = %.2f%s' % (x, r, mark))
    xs2 = None
    for (x1, r1), (x2, r2) in zip(rows2, rows2[1:]):
        if (r1 - 1) * (r2 - 1) <= 0 and r1 != r2:
            xs2 = x1 + (x2 - x1) * (r1 - 1) / (r1 - r2)
            break
    if xs2:
        print('  x* = %.3f   실측 구간 안? %s'
              % (xs2, 'O' if OBS[0] <= xs2 <= OBS[1] else 'X'))
    else:
        print('  1 을 지나지 않음')

    print('\n ⚠ 한계')
    print('  1) ChCl 을 중성 이온쌍으로 다뤄 콜린 -OH 와 Cl- 를 분리하지 못했다.')
    print('     Cl- 수용 면적에 콜린 산소 기여가 섞여 있을 수 있다.')
    print('  2) sigma_hb = 0.0084 는 COSMO-RS 표준값이지 우리 계에 맞춘 값이 아니다.')
    print('  3) 면적 비는 경쟁의 세기(결합 에너지)를 반영하지 않는다. 물은 면적이 커도')
    print('     COOH 보다 약한 공여자다.')

    assert prof, '프로파일을 하나도 못 읽었다'
    print('\n self-check OK')


if __name__ == '__main__':
    main()
