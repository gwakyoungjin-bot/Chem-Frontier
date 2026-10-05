"""MCl3(H2O)3 6배위 vs MCl3(H2O) 4배위 — 어느 쪽이 실제로 안정한가.

이 계산이 답하는 질문 (민이형 지적):
  기존 사다리는 Cl 2개 -> 3개에서 배위수가 6 -> 4 로 한 번에 떨어진다. 그런데
  그 사이의 6배위 MCl3(H2O)3 를 후보에 넣은 적이 없었다. 그래서 "Cl 3개째에서
  사면체로 꺾인다" 가 계산 결과인지 후보 선택의 결과인지 구분이 안 됐다.

비교 반응 (원자·전하 balanced):
  [MCl3(H2O)3]-  ->  [MCl3(H2O)]- + 2 H2O
    dG > 0  6배위가 안정 -> 기존 사다리가 틀렸다. 전환점이 오른쪽으로 밀린다.
    dG < 0  4배위가 안정 -> 기존 결론이 후보 선택의 산물이 아니다.

fac / mer 중 낮은 쪽을 6배위 대표값으로 쓴다.
준조화 보정(|nu|<100 cm-1 치환)을 모든 종에 같은 규칙으로 적용한다.

서버에서 실행: python harvest_mcl3.py
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

import speciation as S                      # read_gibbs_eff (준조화 포함)

HARTREE = 2625.4996                         # kJ/mol
T = 333.15
HERE = os.path.dirname(os.path.abspath(__file__))
SUB = {'Ni': '', 'Co': 'cobalt', 'Mn': 'manganese'}


def g(path):
    p = os.path.join(HERE, path)
    if not os.path.exists(p):
        return None
    return S.read_gibbs_eff(p, T)


def main():
    gw = g('H2O_freq.out')
    if gw is None:
        sys.exit('H2O_freq.out 을 못 찾음 — 서버에서 실행해야 한다')

    print(' MCl3 화학종: 6배위 vs 4배위  (%.2f K, 준조화 보정)' % T)
    print(' %-4s %12s %12s %12s   %s'
          % ('금속', 'fac(6배위)', 'mer(6배위)', '4배위', 'dG(6->4) kJ/mol'))
    print(' ' + '-' * 70)

    rows = []
    for m, sub in SUB.items():
        fac = g(os.path.join('mcl3_oh', '%sCl3_aq3_fac_freq.out' % m))
        mer = g(os.path.join('mcl3_oh', '%sCl3_aq3_mer_freq.out' % m))
        td = g(os.path.join(sub, '%sCl3_aq1_freq.out' % m))
        if None in (fac, mer, td):
            print(' %-4s  (출력 누락: fac=%s mer=%s td=%s)'
                  % (m, fac is not None, mer is not None, td is not None))
            continue
        oh = min(fac, mer)                   # 낮은 이성질체가 대표
        which = 'fac' if fac <= mer else 'mer'
        # [MCl3(H2O)3]- -> [MCl3(H2O)]- + 2 H2O
        dg = (td + 2 * gw - oh) * HARTREE
        rows.append((m, which, dg, (mer - fac) * HARTREE))
        print(' %-4s %12.6f %12.6f %12.6f   %+8.1f'
              % (m, fac, mer, td, dg))

    # 표를 파일로 남긴다. 손으로 다시 계산하면 H2O 를 Gibbs 가 아닌 단일점
    # 에너지로 잘못 쓰기 쉽다 (실제로 한 번 틀렸다 — 4.9 kJ/mol 오차).
    out = os.path.join(os.path.dirname(HERE), 'dft_results', 'mcl3_6vs4.csv')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'w', newline='', encoding='utf-8-sig') as f:
        w = csv.writer(f)
        w.writerow(['metal', 'dG_6to4_kJmol', 'more_stable', 'lower_isomer',
                    'fac_minus_mer_kJmol', 'T_K', 'reaction'])
        for m, which, dg, iso in rows:
            w.writerow([m, '%+.1f' % dg,
                        '4-coordinate (Td)' if dg < 0 else '6-coordinate (Oh)',
                        which, '%+.1f' % iso, '%.2f' % T,
                        '[MCl3(H2O)3]- -> [MCl3(H2O)]- + 2 H2O, quasi-harmonic'])
    print('\n wrote %s' % out)

    print('\n 판정  (dG = [4배위 + 2H2O] - [6배위])')
    for m, which, dg, iso in rows:
        verdict = ('4배위가 안정 — 기존 사다리 유지' if dg < 0 else
                   '**6배위가 안정 — 사다리 수정 필요**')
        print('  %-4s dG = %+7.1f kJ/mol   %s' % (m, dg, verdict))
        print('       안정한 이성질체: %s (fac-mer 차 %+.1f kJ/mol)' % (which, iso))

    if rows and all(d < 0 for _, _, d, _ in rows):
        print('\n  => 세 금속 모두 4배위가 이긴다.')
        print('     "Cl 3개째에서 사면체로 꺾인다" 는 후보 선택의 결과가 아니라')
        print('     계산 결과다. 민이형 지적에 정면으로 답할 수 있다.')
    elif rows and any(d > 0 for _, _, d, _ in rows):
        print('\n  => 6배위가 이기는 금속이 있다. 사다리에 추가하고')
        print('     anchored_map.py 를 다시 돌려야 한다. 지도가 바뀐다.')

    assert rows, '수확된 결과가 없다'
    print('\n self-check OK: %d 금속' % len(rows))


if __name__ == '__main__':
    main()
