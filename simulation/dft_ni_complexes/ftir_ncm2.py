"""1;4 vs 1;4N 차이 스펙트럼 정밀 해석 — 무엇이 금속 때문인지 가른다.

1차 분석(ftir_ncm.py)의 오류를 고친다:
  원 스펙트럼에서 COO- 밴드를 찾으려 했는데, 원 스펙트럼은 자유 시트르산(약 3 M)이
  지배해서 금속이 만든 약한 밴드가 묻힌다. 실제로 탐색창 경계값(1660)이 잡혔고
  그걸로 낸 dv = 266 은 의미가 없다.
  **금속이 만든 밴드는 차이 스펙트럼에서 찾아야 한다.**

가려야 할 두 가설 (둘 다 C=O 감소 + COO- 증가를 낸다):
  (A) 금속-카르복실레이트 직접 배위
  (B) NCM 용해가 양성자를 소모 -> pH 상승 -> 시트르산이 그냥 탈양성자화
      (LiMO2 + H+ -> M(2+) + H2O. 금속 34 mM 이면 H+ 약 100 mM 소모)

  가르는 기준 둘:
    1) dv = COO- 비대칭 - 대칭
         160 미만  이좌/가교 배위       -> (A)
         160~200   자유 카르복실레이트   -> (B)
         200 초과  단좌 배위            -> (A)
    2) 크기. 직접 배위면 생성된 COO- 가 금속 배위권(금속당 최대 6)을 넘을 수 없다.

실행: python dft_ni_complexes/ftir_ncm2.py
"""
import csv
import math
import os
import sys

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

DESK = os.path.join(os.path.expanduser('~'), 'Desktop')
QUIET = (1900.0, 2100.0)     # 흡수 거의 없는 창 -> 베이스라인 0 기준
C_METAL = 0.034              # M. ANCHOR_VALUES.md 의 총 금속 농도
C_COOH = 3 * 4.22            # M. 1:4 의 시트르산 4.22 M x COOH 3개 (ph_estimate.py)


def load(fn):
    """%T CSV -> (파수, 베이스라인 보정된 흡광도)."""
    xs, ts = [], []
    with open(os.path.join(DESK, fn), encoding='utf-8-sig') as f:
        for row in csv.reader(f):
            try:
                xs.append(float(row[0])); ts.append(float(row[1]))
            except (ValueError, IndexError):
                continue                       # 헤더 두 줄
    a = [(-math.log10(t / 100.0) if t > 0 else float('nan')) for t in ts]
    q = [v for x, v in zip(xs, a) if QUIET[0] <= x <= QUIET[1] and v == v]
    base = sum(q) / len(q)
    return xs, [v - base for v in a]


def band(xs, ys, lo, hi, want='max'):
    """구간 내 극값과 '창 경계에 걸렸는지' 를 같이 돌려준다.
    경계값이면 진짜 봉우리가 아니라 옆 밴드의 어깨다 — 1차 분석이 이걸로 틀렸다."""
    w = [(x, y) for x, y in zip(xs, ys) if lo <= x <= hi and y == y]
    p = (max if want == 'max' else min)(w, key=lambda t: t[1])
    return p, (abs(p[0] - lo) < 1.5 or abs(p[0] - hi) < 1.5)


def main():
    x1, a1 = load('1;4.csv')             # NCM 있음
    x0, a0 = load('1;4 N.csv')           # NCM 없음
    assert x1 == x0, '파수 격자가 다르다'

    print(' 스케일 보정 민감도 — 배율을 바꿔도 결론이 버티는가')
    print(' %-9s %-24s %s' % ('배율', '최대 증가', '최대 감소'))
    print(' ' + '-' * 60)
    ups = []
    for k in (1.000, 0.980, 0.967, 0.950):
        d = [(x, p * k - q) for x, p, q in zip(x1, a1, a0) if 1300 <= x <= 1800]
        u = max(d, key=lambda t: t[1]); v = min(d, key=lambda t: t[1])
        ups.append(u[0])
        print(' %-9.3f %+.4f @ %4.0f cm-1      %+.4f @ %4.0f cm-1'
              % (k, u[1], u[0], v[1], v[0]))

    k = 0.967                            # C-O 영역으로 맞춘 ATR 접촉 보정
    dx = x1
    dy = [p * k - q for p, q in zip(a1, a0)]

    print('\n 금속이 만든 밴드 (차이 스펙트럼에서 탐색)')
    asym, e_a = band(dx, dy, 1500, 1650)
    sym,  e_s = band(dx, dy, 1330, 1470)
    co,   _   = band(dx, dy, 1660, 1780, 'min')
    print('  COO- 비대칭 (증가)  %4.0f cm-1  dA = %+.4f %s'
          % (asym[0], asym[1], '<- 창 경계, 신뢰 불가' if e_a else ''))
    print('  COO- 대칭   (증가)  %4.0f cm-1  dA = %+.4f %s'
          % (sym[0], sym[1], '<- 창 경계, 신뢰 불가' if e_s else ''))
    print('  C=O (COOH)  (감소)  %4.0f cm-1  dA = %+.4f' % (co[0], co[1]))

    dv = asym[0] - sym[0]
    verdict = ('이좌/가교 배위 (직접 배위)' if dv < 160 else
               '자유 카르복실레이트 (직접 배위 아님)' if dv <= 200 else
               '단좌 배위 (직접 배위)')
    print('\n  dv = %.0f - %.0f = %.0f cm-1  ->  %s' % (asym[0], sym[0], dv, verdict))

    a_co = max(v for x, v in zip(x0, a0) if 1700 <= x <= 1730)
    frac = abs(co[1]) / a_co
    formed = frac * C_COOH
    print('\n 크기 검증')
    print('  C=O 밴드가 %.2f %% 감소 (dA %.4f / A %.3f)' % (100 * frac, abs(co[1]), a_co))
    print('  COOH 총 %.1f M 의 %.2f %% = 약 %.3f M 이 COO- 로 전환' % (C_COOH, 100 * frac, formed))
    print('  총 금속 %.3f M  ->  금속 1개당 카르복실 %.1f 개' % (C_METAL, formed / C_METAL))
    print('  (직접 배위라면 배위권 최대 6개를 넘을 수 없다)')

    print('\n 판정')
    if dv > 200 or dv < 160:
        print('  dv 가 직접 배위 범위다.')
    else:
        print('  dv 가 자유 카르복실레이트 범위다 — 금속에 직접 붙은 신호가 아니다.')
    if formed / C_METAL > 6:
        print('  생성량이 배위권을 %.0f 배 넘는다. 직접 배위로는 설명이 안 된다.'
              % (formed / C_METAL / 6))
        print('  NCM 용해가 H+ 를 소모해 시트르산이 전반적으로 탈양성자화된 것으로 본다.')

    assert max(ups) - min(ups) < 60, '배율에 따라 밴드 위치가 흔들린다 = 아티팩트'
    print('\n self-check OK: 배율 0.95~1.00 에서 증가 밴드가 %.0f~%.0f cm-1 로 안정'
          % (min(ups), max(ups)))


if __name__ == '__main__':
    main()
