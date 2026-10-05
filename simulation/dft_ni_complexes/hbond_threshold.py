"""FT-IR 전이와 금속 화학종 전이를 하나의 원인으로 설명할 수 있는가.

가설 (공통 원인):
  Cl- 를 붙잡는 것은 COOH 의 수소결합(Cl-...H-O)이다. 그렇다면

    COOH 가 남아돌면  ->  Cl- 가 전부 묶임  ->  자유 Cl- 없음  ->  금속은 아쿠아
    COOH 가 모자라면  ->  자유 Cl- 출현      ->  금속이 사면체로

  그리고 **같은 임계점에서** COOH 자신의 환경도 바뀐다(묶인 COOH 비율이 꺾임)
  -> C=O 파수가 꺾인다.

  즉 FT-IR 이 보는 것과 우리가 계산한 것이 '다른 현상'이 아니라
  **같은 화학량론적 임계의 두 얼굴**일 수 있다.

임계 조성:
  시트르산 1분자에 COOH 3개. 무수 기준 COOH/Cl- = 3(1-x)/x.
  Cl- 하나가 받는 수소결합을 n 개라 하면 임계는  x* = 3/(3+n).
  n 을 가정값으로 두지 않고 2~4 범위를 전부 보고 실측 전이와 비교한다.

이건 DFT 결과가 아니라 **화학량론 논증**이다. 그렇게 표기해서 써야 한다.

실행: python dft_ni_complexes/hbond_threshold.py
"""
import math
import sys

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

# 실측 전이 구간 (앞선 분석에서)
OBS = {
    '육안 · 금속 화학종': (0.40, 0.60),
    'FT-IR · 용매 구조': (0.50, 0.60),
    '계산 · 금속 화학종': (0.45, 0.50),
}
# C=O 국소 기울기 이상치가 시작된 지점 (ftir_recheck.py)
KINK = 0.50

N_HB = [2, 2.5, 3, 3.5, 4]      # Cl- 하나가 받는 수소결합 수 후보
N_COOH = 3                      # 시트르산의 카르복실기 수


def x_star(n):
    """COOH/Cl- = n 이 되는 조성. 3(1-x)/x = n -> x = 3/(3+n)"""
    return N_COOH / (N_COOH + n)


def ratio(x):
    return N_COOH * (1 - x) / x if x > 0 else float('inf')


def main():
    print(' 실측 전이 구간')
    for k, (a, b) in OBS.items():
        print('   %-22s %.2f ↔ %.2f' % (k, a, b))
    lo = min(a for a, _ in OBS.values())
    hi = max(b for _, b in OBS.values())
    print('   %-22s %.2f ↔ %.2f' % ('-> 전체 포괄', lo, hi))

    print('\n 화학량론 임계  x* = 3/(3+n)')
    print(' %-8s %-10s %s' % ('n', 'x*', '실측 구간 안?'))
    print(' ' + '-' * 38)
    hits = []
    for n in N_HB:
        xs = x_star(n)
        ok = lo <= xs <= hi
        hits.append((n, xs, ok))
        print(' %-8.1f %-10.3f %s' % (n, xs, 'O' if ok else 'X'))

    print('\n 조성별 COOH / Cl- 비  (무수 기준)')
    print(' %-8s %-10s %s' % ('x(ChCl)', 'COOH/Cl-', ''))
    print(' ' + '-' * 36)
    for x in (0.05, 0.20, 0.333, 0.40, 0.45, 0.50, 0.60, 0.75, 0.95):
        r = ratio(x)
        mark = ''
        if abs(x - KINK) < 1e-9:
            mark = '  <== FT-IR 꺾임'
        elif x == 0.45:
            mark = '  <== 계산 f_Td 상승 시작'
        print(' %-8.2f %-10.2f%s' % (x, r, mark))

    print('\n [핵심]')
    r_kink = ratio(KINK)
    print('  FT-IR 이 꺾이는 x = %.2f 에서 COOH/Cl- = %.2f 이다.' % (KINK, r_kink))
    print('  즉 Cl- 하나당 카르복실기가 정확히 %.0f개인 지점이다.' % r_kink)
    print('  이 아래(COOH 과잉)에서는 Cl- 가 전부 묶이고,')
    print('  이 위(COOH 부족)에서는 묶이지 못한 Cl- 가 생긴다.')

    ok_n = [n for n, xs, ok in hits if ok]
    print('\n [판정]')
    if ok_n:
        print('  Cl- 당 수소결합 %s개를 가정하면 임계가 실측 전이 구간(%.2f~%.2f) 안에 든다.'
              % ('~'.join('%.1f' % n for n in (min(ok_n), max(ok_n))), lo, hi))
        print('  => 두 전이를 **하나의 화학량론적 임계**로 설명할 수 있다.')
        print('     FT-IR(용매)과 계산(금속)은 서로 인과가 아니라 **공통 원인의 두 관측**이다.')
    else:
        print('  어떤 n 을 가정해도 임계가 실측 구간 밖이다. 이 설명은 성립하지 않는다.')

    print('\n ⚠ 한계 — 이대로 주장하면 안 되는 이유')
    print('  1) 이건 화학량론 논증이지 DFT 결과가 아니다. 계산이 이 임계를 예측한 게 아니다.')
    print('  2) 콜린 양이온의 -OH 도 Cl- 에 수소결합한다. x 가 크면 콜린이 많아지므로')
    print('     실제 임계는 위 값보다 오른쪽으로 밀린다. 보정 안 했다.')
    print('  3) 물이 몰분율 0.80 이다. 물도 Cl- 에 수소결합하므로 무수 기준 비는 근사다.')
    print('  4) n 을 맞춰서 구간에 넣은 것이라 예측이 아니라 사후 설명이다.')

    assert any(ok for _, _, ok in hits), '임계가 실측 구간과 전혀 안 맞는다'
    print('\n self-check OK: n = %s 에서 부합' % ok_n)


if __name__ == '__main__':
    main()
