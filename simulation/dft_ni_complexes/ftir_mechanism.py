"""FT-IR 이 우리 메커니즘의 어느 고리를 지지하는지 검증한다.

우리가 주장하는 사슬:
    조성비  ->  염화물 활동도  ->  금속 화학종(Oh -> Td)
              (A)                (B)

FT-IR 은 시트르산의 C=O 와 O-H 를 본다. 금속 화학종을 직접 보지 않는다.
따라서 FT-IR 이 지지할 수 있는 건 잘해야 (A) 고리다. 그것도 조건이 있다.

**판정 기준**
  우리 COSMO-RS 는 a_ChCl 이 조성에 대해 매끄럽게(log 선형에 가깝게) 오른다고 한다.
  만약 C=O 이동이 (A) 때문이라면, C=O 는 log a_ChCl 에 대해 **매끄럽게** 움직여야 한다.
  반대로 C=O 가 특정 조성에서 **계단**을 보이면, 그건 염화물 활동도가 아니라
  DES 구조 자체의 전이를 보는 것이다 -> (A) 의 증거로 쓸 수 없다.

  구체적으로 세 모형을 같은 데이터에 맞춰 R^2 를 비교한다:
    M1  C=O ~ x(ChCl)            단순 희석/혼합
    M2  C=O ~ log10 a_ChCl       우리 메커니즘 (A)
    M3  C=O ~ 계단(x >= x0)      구조 전이
  어느 것도 압도적이지 않으면 FT-IR 은 (A) 의 증거가 되지 못한다.

O-H 영역도 본다. Cl- 는 IR 밴드가 없지만 Cl-...H-O 수소결합이 O-H 를 민다.
따라서 O-H 가 염화물의 유일한 간접 지표다.

실행: python dft_ni_complexes/ftir_mechanism.py
입력: 바탕화면 '새 스프레드시트 문서.xlsx', dft_results/ternary_grid.csv
"""
import csv
import math
import os
import re
import sys
import zipfile

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

XLSX = os.path.join(os.path.expanduser('~'), 'Desktop', '새 스프레드시트 문서.xlsx')
GRID = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    'dft_results', 'ternary_grid.csv')
CO_WIN, OH_WIN, QUIET = (1650.0, 1800.0), (2900.0, 3700.0), (1900.0, 2100.0)
T_GRID = 298.15


def parse_label(v):
    if v.strip().lower() in ('1.4n', '1,4n'):
        return None, None
    s = ('%.2f' % float(v)).rstrip('0')
    a, b = s.split('.')
    a, b = int(a), int(b)
    return '%d;%d' % (a, b), a / (a + b)


def read_sheet():
    z = zipfile.ZipFile(XLSX)
    ss = re.findall(r'<t[^>]*>(.*?)</t>',
                    z.read('xl/sharedStrings.xml').decode('utf-8'), re.S) \
        if 'xl/sharedStrings.xml' in z.namelist() else []
    x = z.read('xl/worksheets/sheet1.xml').decode('utf-8')
    grid = {}
    for rn, body in re.findall(r'<row[^>]*r="(\d+)"[^>]*>(.*?)</row>', x, re.S):
        for c in re.findall(r'<c\b[^>]*>.*?</c>|<c\b[^>]*/>', body, re.S):
            ref = re.search(r'r="([A-Z]+)\d+"', c)
            v = re.search(r'<v>(.*?)</v>', c, re.S)
            if ref and v:
                val = ss[int(v.group(1))] if re.search(r't="s"', c) else v.group(1)
                grid.setdefault(int(rn), {})[ref.group(1)] = val
    return grid


def load_activity():
    tab = {}
    with open(GRID, encoding='utf-8') as f:
        for r in csv.DictReader(f):
            if abs(float(r['x_H2O']) - 0.8) < 1e-9 and \
               abs(float(r['T_K']) - T_GRID) < 1e-6:
                tab[round(float(r['x_ChCl_dry']), 6)] = float(r['a_ChCl'])
    return tab


def interp_log(tab, x):
    ks = sorted(tab)
    if x <= ks[0]:
        return tab[ks[0]]
    if x >= ks[-1]:
        return tab[ks[-1]]
    for lo, hi in zip(ks, ks[1:]):
        if lo <= x <= hi:
            t = (x - lo) / (hi - lo)
            return math.exp(math.log(tab[lo]) * (1 - t) + math.log(tab[hi]) * t)


def peak(wl, a, lo, hi):
    idx = [i for i, w in enumerate(wl) if lo <= w <= hi]
    i = max(idx, key=lambda j: a[j])
    if i in (0, len(wl) - 1):
        return wl[i]
    y0, y1, y2 = a[i - 1], a[i], a[i + 1]
    d = y0 - 2 * y1 + y2
    return wl[i] + (0.5 * (y0 - y2) / d if d else 0.0) * (wl[i + 1] - wl[i])


def centroid(wl, a, lo, hi):
    """O-H 는 봉우리가 뭉툭해 위치가 불안정하다. 무게중심이 더 안정적이다."""
    pts = [(w, max(v, 0.0)) for w, v in zip(wl, a) if lo <= w <= hi]
    s = sum(v for _, v in pts)
    return sum(w * v for w, v in pts) / s if s > 0 else float('nan')


def fit(xs, ys):
    """최소제곱 직선과 R^2."""
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((a - mx) ** 2 for a in xs)
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    if sxx == 0:
        return 0.0, 0.0, 0.0
    m = sxy / sxx
    c = my - m * mx
    ss_res = sum((b - (m * a + c)) ** 2 for a, b in zip(xs, ys))
    ss_tot = sum((b - my) ** 2 for b in ys)
    return m, c, (1 - ss_res / ss_tot if ss_tot else 0.0)


def main():
    grid = read_sheet()
    hdr = grid[1]
    cols = {}
    for col, v in hdr.items():
        lab, xc = parse_label(v)
        if xc is not None:
            cols[col] = (lab, xc)

    rows = sorted(r for r in grid if r >= 2)
    wl = [float(grid[r]['A']) for r in rows if 'A' in grid[r]]

    act = load_activity()
    rec = []
    for col, (lab, xc) in cols.items():
        ts = [float(grid[r].get(col)) if grid[r].get(col) else float('nan')
              for r in rows]
        a = [(-math.log10(t / 100.0) if t == t and t > 0 else float('nan')) for t in ts]
        q = [v for w, v in zip(wl, a) if QUIET[0] <= w <= QUIET[1] and v == v]
        a = [v - sum(q) / len(q) for v in a]
        rec.append(dict(lab=lab, x=xc, co=peak(wl, a, *CO_WIN),
                        oh=centroid(wl, a, *OH_WIN),
                        la=math.log10(interp_log(act, xc))))
    rec.sort(key=lambda r: r['x'])

    print(' 조성별 실측값과 계산 활동도')
    print(' %-6s %-8s %-11s %-11s %s' % ('시료', 'x(ChCl)', 'C=O /cm-1',
                                          'O-H 무게중심', 'log10 a_ChCl'))
    print(' ' + '-' * 60)
    for r in rec:
        print(' %-6s %-8.3f %-11.1f %-11.1f %.2f'
              % (r['lab'], r['x'], r['co'], r['oh'], r['la']))

    xs = [r['x'] for r in rec]
    las = [r['la'] for r in rec]
    cos = [r['co'] for r in rec]
    ohs = [r['oh'] for r in rec]

    print('\n [검증 1] 우리 활동도는 조성에 대해 매끄러운가')
    m, c, r2 = fit(xs, las)
    print('  log10 a_ChCl ~ x   R^2 = %.4f  (1 에 가까우면 계단 없이 매끄럽다)' % r2)

    print('\n [검증 2] C=O 는 무엇을 따라가는가  (R^2 클수록 그 모형이 맞다)')
    for name, xv in (('M1  x(ChCl) 단순', xs), ('M2  log a_ChCl (우리 메커니즘)', las)):
        _, _, r2c = fit(xv, cos)
        print('  %-30s R^2 = %.4f' % (name, r2c))
    best = None
    for x0 in [v for v in sorted(set(xs)) if 0.1 < v < 0.95]:
        step = [1.0 if v >= x0 else 0.0 for v in xs]
        _, _, r2s = fit(step, cos)
        if best is None or r2s > best[1]:
            best = (x0, r2s)
    print('  %-30s R^2 = %.4f   (문턱 x0 = %.2f)'
          % ('M3  계단 모형', best[1], best[0]))

    print('\n [검증 3] O-H 는 염화물을 보는가  (Cl-...H-O 가 유일한 간접 지표)')
    _, _, r2o_x = fit(xs, ohs)
    _, _, r2o_a = fit(las, ohs)
    print('  O-H ~ x(ChCl)       R^2 = %.4f' % r2o_x)
    print('  O-H ~ log a_ChCl    R^2 = %.4f' % r2o_a)

    print('\n [판정]')
    _, _, r2_lin = fit(las, cos)
    if best[1] > r2_lin + 0.05:
        print('  C=O 는 계단 모형이 더 잘 맞는다 (R^2 %.3f vs %.3f).'
              % (best[1], r2_lin))
        print('  우리 활동도는 조성에 대해 매끄럽게(R^2 %.3f) 오르는데 C=O 는 꺾인다.' % r2)
        print('  => FT-IR 의 C=O 는 **염화물 활동도의 증거가 아니다**.')
        print('     DES 구조 자체가 x = %.2f 부근에서 재편된다는 별개의 증거다.' % best[0])
    else:
        print('  C=O 가 log a_ChCl 을 잘 따라간다 (R^2 %.3f).' % r2_lin)
        print('  => 조성 -> 염화물 활동도 고리의 실측 근거로 쓸 수 있다.')

    assert len(rec) >= 10, '조성이 부족하다'
    print('\n self-check OK: %d조성' % len(rec))


if __name__ == '__main__':
    main()
