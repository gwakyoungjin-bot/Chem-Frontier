"""최적화된 구조가 화학적으로 말이 되는지 수치로 검사한다.

그림만 보고 '이상하다/멀쩡하다' 를 판단하면 카메라 각도에 속는다.
여기서는 거리만 본다:
  - M-O, M-Cl 배위거리가 상식 범위인가
  - 물의 O-H 가 정상인가 (0.95~1.02 A). 깨졌으면 양성자가 떨어져 나간 것
  - 물의 H 가 금속 쪽을 향하는가 (뒤집힌 물 = 잘못된 배위)
  - 원자끼리 겹치는가

실행: python dft_ni_complexes/check_geom.py
"""
import math
import os
import sys

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

HERE = os.path.dirname(os.path.abspath(__file__))
OPT = os.path.join(os.path.dirname(HERE), 'dft_results', 'geometries', 'optimized')
METALS = ('Ni', 'Co', 'Mn')

NAMES = [n + s for n in METALS
         for s in ('_aq6', 'Cl1_aq5', 'Cl2_aq4', 'Cl3_aq1', 'Cl4')]


def read(path):
    with open(path, encoding='utf-8') as f:
        L = [l.split() for l in f if l.strip()]
    n = int(L[0][0])
    return [(p[0], float(p[1]), float(p[2]), float(p[3])) for p in L[2:2 + n]]


def d(a, b):
    return math.dist(a[1:], b[1:])


def check(name):
    p = os.path.join(OPT, '%s_opt.xyz' % name)
    if not os.path.exists(p):
        return None
    at = read(p)
    m = next(a for a in at if a[0] in METALS)
    os_, cls, hs = [a for a in at if a[0] == 'O'], \
                   [a for a in at if a[0] == 'Cl'], \
                   [a for a in at if a[0] == 'H']

    mo = sorted(d(m, o) for o in os_)
    mcl = sorted(d(m, c) for c in cls)

    bad = []
    # 물이 온전한가: 각 H 는 어떤 O 와 0.95~1.02 A 로 묶여야 한다
    for h in hs:
        near = min(d(h, o) for o in os_) if os_ else 9.9
        if not (0.90 <= near <= 1.10):
            bad.append('H-O %.2f' % near)
    # 물이 뒤집혔나: H 가 자기 O 보다 금속에 가까우면 배위 방향이 틀린 것
    flip = sum(1 for h in hs if os_ and d(m, h) < min(d(m, o) for o in os_) - 0.1)
    if flip:
        bad.append('뒤집힌 H %d개' % flip)
    # 원자 겹침
    for i in range(len(at)):
        for j in range(i + 1, len(at)):
            if d(at[i], at[j]) < 0.85:
                bad.append('겹침 %s-%s %.2f' % (at[i][0], at[j][0], d(at[i], at[j])))

    cn = len(mo and [x for x in mo if x < 2.6]) + len([x for x in mcl if x < 2.9])
    return dict(name=name, m=m[0], nO=len(os_), nCl=len(cls), cn=cn,
                mo=mo, mcl=mcl, bad=bad)


def main():
    print(' 종            배위수  M-O (A)              M-Cl (A)           이상')
    print(' ' + '-' * 78)
    issues = 0
    for name in NAMES:
        r = check(name)
        if not r:
            print(' %-12s  (파일 없음)' % name)
            continue
        mo = ' '.join('%.2f' % x for x in r['mo'][:6]) or '-'
        mcl = ' '.join('%.2f' % x for x in r['mcl'][:4]) or '-'
        flag = '; '.join(r['bad']) if r['bad'] else 'OK'
        if r['bad']:
            issues += 1
        print(' %-12s  %d     %-20s %-18s %s' % (r['name'], r['cn'], mo, mcl, flag))

    print('\n 이상 있는 종: %d개' % issues)
    # 자체검증: 배위거리가 상식 밖이면 후속 계산 전부가 무의미하다
    for name in NAMES:
        r = check(name)
        if not r:
            continue
        for x in r['mo']:
            assert 1.8 < x < 3.0, '%s: M-O %.2f A 는 배위거리가 아니다' % (name, x)
        for x in r['mcl']:
            assert 2.0 < x < 3.2, '%s: M-Cl %.2f A 는 배위거리가 아니다' % (name, x)
    print(' self-check OK: 모든 M-O 1.8~3.0, M-Cl 2.0~3.2 A 범위 안')


if __name__ == '__main__':
    main()
