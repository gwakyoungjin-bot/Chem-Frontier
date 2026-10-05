"""앞선 판정(ftir_mechanism.py)을 다시 검증한다 — R^2 비교가 정당했는지부터.

앞 분석의 의심스러운 점:
  M1 (C=O ~ x) 과 M2 (C=O ~ log a_ChCl) 를 R^2 로 비교했는데,
  log a 자체가 x 에 대해 거의 선형(R^2 = 0.978)이다. 즉 두 예측변수가 **공선적**이다.
  공선적인 두 모형을 R^2 차이로 가르는 건 불안정하다. 다시 해야 한다.

대신 쓰는 검정 (모형 가정이 없거나 약한 것들):
  T1  국소 기울기 이상치 — 인접 구간의 d(C=O)/d(log a) 분포에서 튀는 구간이 있는가.
      a_ChCl 이 원인이라면 이 기울기가 대체로 일정해야 한다.
  T2  분할회귀 — 꺾임을 허용하면 적합이 유의하게 좋아지는가 (F 검정).
  T3  한 점 빼기(leave-one-out) — 결론이 특정 시료 하나에 의존하는가.
      이게 제일 중요하다. 계단이 한 쌍(0.50 vs 0.60)에만 기대고 있으면 약한 근거다.

실행: python dft_ni_complexes/ftir_recheck.py
"""
import math
import sys

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

# ftir_mechanism.py 출력에서 그대로 (실측 C=O, 계산 log10 a_ChCl)
DATA = [
    ("1;19", 0.050, 1713.9, -7.91), ("1;9", 0.100, 1714.0, -7.27),
    ("3;17", 0.150, 1714.4, -6.76), ("1;4", 0.200, 1714.7, -6.31),
    ("1;3", 0.250, 1714.8, -5.90), ("1;2", 0.333, 1715.3, -5.26),
    ("2;3", 0.400, 1715.6, -4.79), ("1;1", 0.500, 1716.0, -4.15),
    ("3;2", 0.600, 1722.3, -3.56), ("2;1", 0.667, 1722.9, -3.21),
    ("3;1", 0.750, 1723.5, -2.81), ("17;3", 0.850, 1726.8, -2.39),
    ("19;1", 0.950, 1729.2, -2.02),
]


def fit(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((a - mx) ** 2 for a in xs)
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    m = sxy / sxx if sxx else 0.0
    c = my - m * mx
    ss_res = sum((b - (m * a + c)) ** 2 for a, b in zip(xs, ys))
    ss_tot = sum((b - my) ** 2 for b in ys)
    return m, c, ss_res, (1 - ss_res / ss_tot if ss_tot else 0.0)


def seg_fit(xs, ys, brk):
    """꺾은선: 기울기가 brk 전후로 달라지는 모형. 파라미터 3개."""
    A = [(1.0, x, max(0.0, x - brk)) for x in xs]
    # 정규방정식 3x3 을 가우스 소거로
    n = 3
    M = [[sum(A[k][i] * A[k][j] for k in range(len(xs))) for j in range(n)] +
         [sum(A[k][i] * ys[k] for k in range(len(xs)))] for i in range(n)]
    for i in range(n):
        p = max(range(i, n), key=lambda r: abs(M[r][i]))
        M[i], M[p] = M[p], M[i]
        if abs(M[i][i]) < 1e-12:
            return None, None
        for r in range(n):
            if r != i:
                f = M[r][i] / M[i][i]
                for c in range(i, n + 1):
                    M[r][c] -= f * M[i][c]
    b = [M[i][n] / M[i][i] for i in range(n)]
    pred = [b[0] + b[1] * x + b[2] * max(0.0, x - brk) for x in xs]
    ss = sum((p - y) ** 2 for p, y in zip(pred, ys))
    return b, ss


def main():
    labs = [d[0] for d in DATA]
    xs = [d[1] for d in DATA]
    cos = [d[2] for d in DATA]
    las = [d[3] for d in DATA]

    print(' [T0] 예측변수가 공선적인가 — 그렇다면 R^2 비교는 무효')
    _, _, _, r2 = fit(xs, las)
    print('  log a ~ x  R^2 = %.4f' % r2)
    print('  -> %s\n' % ('공선적. 앞선 M1 vs M2 비교는 신뢰할 수 없다.'
                         if r2 > 0.95 else '충분히 독립적'))

    print(' [T1] 국소 기울기  d(C=O)/d(log a)  — a_ChCl 이 원인이면 일정해야 한다')
    sl = []
    for i in range(len(DATA) - 1):
        dla = las[i + 1] - las[i]
        s = (cos[i + 1] - cos[i]) / dla if dla else float('nan')
        sl.append((s, labs[i], labs[i + 1]))
    med = sorted(s for s, _, _ in sl)[len(sl) // 2]
    for s, a, b in sl:
        flag = '  <== 이상치 (중앙값의 %.0f배)' % (s / med) if s > 4 * med else ''
        print('  %-6s -> %-6s  %6.2f cm-1 per log a%s' % (a, b, s, flag))
    print('  중앙값 %.2f' % med)

    print('\n [T2] 꺾임을 허용하면 유의하게 좋아지는가 (F 검정)')
    _, _, ss1, _ = fit(las, cos)
    best = None
    for brk_i in range(2, len(xs) - 2):
        b, ss2 = seg_fit(las, cos, las[brk_i])
        if ss2 is not None and (best is None or ss2 < best[1]):
            best = (labs[brk_i], ss2, xs[brk_i])
    n = len(xs)
    F = ((ss1 - best[1]) / 1) / (best[1] / (n - 3))
    print('  직선 SS = %.2f   꺾은선 SS = %.2f   꺾임 위치 %s (x=%.2f)'
          % (ss1, best[1], best[0], best[2]))
    print('  F(1, %d) = %.1f   %s' % (n - 3, F,
          '-> 꺾임이 유의하다 (F > 10 이면 강함)' if F > 10 else '-> 유의하지 않다'))

    print('\n [T3] 한 점 빼기 — 결론이 시료 하나에 의존하는가')
    print('  %-8s %-12s %-12s %s' % ('제외', '직선 R^2', '꺾은선 SS감소', '결론 유지'))
    print('  ' + '-' * 52)
    flip = []
    for k in range(len(DATA)):
        sub = [d for i, d in enumerate(DATA) if i != k]
        sx = [d[1] for d in sub]; sc = [d[2] for d in sub]; sa = [d[3] for d in sub]
        _, _, s1, r2k = fit(sa, sc)
        bb = None
        for bi in range(2, len(sx) - 2):
            b, s2 = seg_fit(sa, sc, sa[bi])
            if s2 is not None and (bb is None or s2 < bb):
                bb = s2
        drop = (s1 - bb) / s1 if s1 else 0
        keep = drop > 0.5
        if not keep:
            flip.append(labs[k])
        print('  %-8s %-12.3f %-12.1f%% %s'
              % (labs[k], r2k, 100 * drop, 'O' if keep else 'X  <== 뒤집힘'))

    print('\n [최종 판정]')
    if not flip:
        print('  어느 시료를 빼도 꺾임이 유지된다. 한 점에 의존하지 않는다.')
        print('  => C=O 의 불연속은 실재한다. a_ChCl 은 매끄러우므로 둘은 다른 현상이다.')
        print('     FT-IR 은 염화물 활동도의 직접 증거가 **아니다** — 앞선 결론 유지.')
    else:
        print('  %s 을(를) 빼면 결론이 뒤집힌다.' % ', '.join(flip))
        print('  => 근거가 특정 시료에 의존한다. 단정하면 안 된다.')

    assert len(DATA) == 13
    print('\n self-check OK: 13조성, 한점빼기 %d회' % len(DATA))


if __name__ == '__main__':
    main()
