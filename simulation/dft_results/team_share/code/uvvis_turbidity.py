"""무색 조성끼리 UV-Vis 세기를 비교하기 위한 산란 보정.

왜 무색끼리만 비교하나:
  Beer-Lambert 로 A 를 농도 대리로 쓰려면 비교 대상의 **몰흡광계수가 같아야** 한다.
  계산상 무색 6조성(x <= 0.40)은 f_Td = 0 으로 전부 같은 팔면체 종이다.
  청록 4조성은 사면체라 eps 가 두 자릿수 다르므로 같은 축에서 비교할 수 없다.

왜 보정이 필요한가:
  미용해 NCM 입자가 빛을 산란시켜 전 파장에 걸쳐 A 를 올린다. 실제로 1;4 는
  A(800) = 0.322 로 장파장에서 오히려 **증가**한다. 전자전이는 장파장에서 0 으로
  수렴하므로, 이건 흡수가 아니라 탁도다. 빼지 않으면 '많이 녹았다' 로 오독한다.

보정 방법:
  700~800 nm 는 이 계에 전자전이가 없다고 보고, 그 구간에 A = k * lam^(-n) 을
  맞춰 전 파장으로 외삽한 뒤 뺀다. n 은 입자 크기에 따라 0(큰 입자)~4(Rayleigh)
  사이라 고정하지 않고 시료마다 맞춘다.

실행: python dft_ni_complexes/uvvis_turbidity.py
"""
import csv
import math
import os
import sys

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

PATH = os.path.join(os.path.expanduser('~'), 'Desktop', '최종.csv')
FIT = (700.0, 800.0)          # 전자전이가 없다고 보는 창
BAND = (330.0, 370.0)         # 관심 밴드
# 계산상 f_Td = 0 인 조성만 (같은 화학종 -> eps 동일 가정 가능)
PALE = ['1;19', '1;9', '3;17', '1;4', '1;3', '1;2', '2;3']
X = {'1;19': .05, '1;9': .10, '3;17': .15, '1;4': .20, '1;3': .25,
     '1;2': .333, '2;3': .40}
CA_M = {'1;19': 5.05, '1;9': 4.72, '3;17': 4.49, '1;4': 4.22, '1;3': 3.98,
        '1;2': 3.49, '2;3': 3.07}          # ph_estimate.py 기준 시트르산 농도


def load():
    with open(PATH, encoding='utf-8-sig') as f:
        rows = list(csv.reader(f))
    labs = [h.strip() for h in rows[0][1:] if h.strip()]
    wl, data = [], {l: [] for l in labs}
    for r in rows[1:]:
        if not r or not r[0].strip():
            continue
        wl.append(float(r[0]))
        for i, l in enumerate(labs, 1):
            data[l].append(float(r[i]) if r[i].strip() else float('nan'))
    return wl, data


def fit_scatter(wl, y):
    """log A = log k - n log lam 을 최소제곱으로. A<=0 인 점은 로그를 못 쓰니 뺀다."""
    pts = [(math.log(w), math.log(a)) for w, a in zip(wl, y)
           if FIT[0] <= w <= FIT[1] and a == a and a > 1e-4]
    if len(pts) < 20:
        return None, None
    n = len(pts)
    sx = sum(p[0] for p in pts); sy = sum(p[1] for p in pts)
    sxx = sum(p[0] ** 2 for p in pts); sxy = sum(p[0] * p[1] for p in pts)
    den = n * sxx - sx * sx
    if abs(den) < 1e-12:
        return None, None
    slope = (n * sxy - sx * sy) / den
    inter = (sy - slope * sx) / n
    return math.exp(inter), -slope          # k, n


def peak(wl, y, lo, hi):
    w = [(x, v) for x, v in zip(wl, y) if lo <= x <= hi and v == v]
    return max(w, key=lambda t: t[1]) if w else (None, None)


def main():
    wl, data = load()
    print(' 무색 조성만 — 산란 보정 전후 (같은 화학종이라 eps 비교 가능)')
    print(' %-6s %-7s %-9s %-7s %-10s %-10s %s'
          % ('시료', 'x(ChCl)', 'A(800)', 'n', 'A raw', 'A 보정', '보정 후 순위'))
    print(' ' + '-' * 74)

    out = []
    for l in PALE:
        if l not in data:
            continue
        y = data[l]
        k, n = fit_scatter(wl, y)
        raw_w, raw_a = peak(wl, y, *BAND)
        if k is None:
            # 장파장이 음수/0 -> 산란 없음으로 보고 보정 생략
            corr = raw_a
            n = float('nan')
        else:
            base = k * raw_w ** (-n)
            corr = raw_a - base
        a800 = [v for w, v in zip(wl, y) if abs(w - 800) < 0.6][0]
        out.append((l, raw_a, corr, a800, n, raw_w))

    order = sorted(out, key=lambda t: -t[2])
    rank = {t[0]: i + 1 for i, t in enumerate(order)}
    for l, raw, corr, a800, n, w in out:
        print(' %-6s %-7.2f %-9.3f %-7s %-10.3f %-10.3f %d위'
              % (l, X[l], a800, ('%.1f' % n) if n == n else '—', raw, corr, rank[l]))

    print('\n 보정 후 순서 (많이 녹은 쪽 -> 적게)')
    print('  ' + '  >  '.join('%s(%.3f)' % (l, c) for l, _, c, _, _, _ in order))

    # 시트르산 농도와의 상관 — 산 농도가 침출을 끌었다면 양의 상관이어야
    xs = [CA_M[l] for l, _, _, _, _, _ in out]
    ys = [c for _, _, c, _, _, _ in out]
    n_ = len(xs)
    mx, my = sum(xs) / n_, sum(ys) / n_
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = sum((a - mx) ** 2 for a in xs); syy = sum((b - my) ** 2 for b in ys)
    r = sxy / math.sqrt(sxx * syy) if sxx * syy > 0 else float('nan')
    print('\n 시트르산 농도 vs 보정 흡광도 상관계수  r = %+.2f' % r)
    if r > 0.5:
        print('  -> 산이 많을수록 많이 녹았다는 해석과 일치')
    elif r < -0.5:
        print('  -> 산 농도와 반대. 다른 요인이 지배한다')
    else:
        print('  -> 뚜렷한 상관 없음. 산 농도만으로는 설명 안 된다')

    print('\n ⚠ 남은 전제 두 가지')
    print('  1) NCM 없는 같은 조성 DES 의 UV-Vis 가 없어 DES 자체 흡수를 못 뺐다.')
    print('     345 nm 에 시트르산 기여가 있으면 이 순서가 바뀔 수 있다.')
    print('  2) 보정은 700~800 nm 외삽이다. 그 구간에도 약한 흡수가 있으면 과보정된다.')

    assert len(out) >= 5, '비교할 무색 조성이 부족하다'
    print('\n self-check OK: 무색 %d조성' % len(out))


if __name__ == '__main__':
    main()
