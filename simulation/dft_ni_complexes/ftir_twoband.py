"""C=O 영역의 두 하위밴드 세기비로 시트르산 화학종 분포를 구한다.

왜 이게 앞선 분석들보다 나은가:
  지금까지는 C=O 의 **최대점 위치**를 추적했다. 그런데 검증에서 드러난 것은
  1650~1800 창에 **위치가 거의 고정된 두 밴드**(약 1716, 약 1722)가 공존하고
  **상대 세기만** 조성에 따라 바뀐다는 사실이다. 최대점이 6 cm-1 "점프" 한 것은
  두 밴드의 세기가 0.4 % 차이로 뒤바뀐 결과였다.

  진동분광에서 위치 고정 + 세기비 변화는 **두 화학종의 평형**의 표준 지문이다.
  그러면 올바른 관측량은 파수가 아니라 **세기비 = 화학종 분율**이다.

이 접근이 앞선 반론들을 피하는 이유:
  - 랩 간 절대 파수 비교를 하지 않는다 (한 스펙트럼 안의 비율만 쓴다)
  - 베이스라인 상수 오프셋이 비율에서 상쇄된다
  - "계단" 아티팩트가 애초에 생기지 않는다 (뒤바뀜을 일으킨 양 자체를 잰다)

정량 검정:
  두 상태 평형  A + n Cl-  <->  B  이면
      ln( f_B / (1 - f_B) )  =  ln K  +  n · ln a_Cl
  즉 **기울기 n 이 결합한 염화물 개수**다. n 이 1 근처로 나오면
  "시트르산 환경 변화가 염화물 1개 결합으로 일어난다" 는 정량적 결론이 된다.
  n 이 0 이나 음수면 염화물과 무관하다는 뜻이다. 반증 가능하다.

실행: python dft_ni_complexes/ftir_twoband.py
"""
import csv
import math
import os
import re
import sys
import zipfile

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
XLSX = os.path.join(os.path.expanduser('~'), 'Desktop', '새 스프레드시트 문서.xlsx')
GRID = os.path.join(ROOT, 'dft_results', 'ternary_grid.csv')
LO_WIN, HI_WIN = (1705.0, 1719.0), (1719.0, 1734.0)   # 두 하위밴드 창
QUIET = (1900.0, 2100.0)
T = 298.15
# 60 C 계산 사면체 분율 (ftd_60C.csv / ftd_scan)
FTD = {0.05: 0.0, 0.10: 0.0, 0.15: 0.0, 0.20: 0.0, 0.25: 0.0, 0.333: 0.0,
       0.40: 0.0, 0.50: 0.0002, 0.60: 0.0053, 0.667: 0.0347, 0.75: 0.2424,
       0.85: 0.7668, 0.95: 0.9614}


def parse_label(v):
    if v.strip().lower() in ('1.4n', '1,4n'):
        return None, None
    s = ('%.2f' % float(v)).rstrip('0')
    a, b = s.split('.')
    return '%d;%d' % (int(a), int(b)), int(a) / (int(a) + int(b))


def read_sheet():
    z = zipfile.ZipFile(XLSX)
    ss = re.findall(r'<t[^>]*>(.*?)</t>',
                    z.read('xl/sharedStrings.xml').decode('utf-8'), re.S) \
        if 'xl/sharedStrings.xml' in z.namelist() else []
    x = z.read('xl/worksheets/sheet1.xml').decode('utf-8')
    g = {}
    for rn, body in re.findall(r'<row[^>]*r="(\d+)"[^>]*>(.*?)</row>', x, re.S):
        for c in re.findall(r'<c\b[^>]*>.*?</c>|<c\b[^>]*/>', body, re.S):
            ref = re.search(r'r="([A-Z]+)\d+"', c)
            v = re.search(r'<v>(.*?)</v>', c, re.S)
            if ref and v:
                val = ss[int(v.group(1))] if re.search(r't="s"', c) else v.group(1)
                g.setdefault(int(rn), {})[ref.group(1)] = val
    return g


def load_act():
    tab = {}
    with open(GRID, encoding='utf-8') as f:
        for r in csv.DictReader(f):
            if abs(float(r['x_H2O']) - 0.8) < 1e-9 and abs(float(r['T_K']) - T) < 1e-6:
                tab[round(float(r['x_ChCl_dry']), 6)] = float(r['a_ChCl'])
    return tab


def a_of(tab, x):
    ks = sorted(tab)
    if x <= ks[0]:
        return tab[ks[0]]
    if x >= ks[-1]:
        return tab[ks[-1]]
    for lo, hi in zip(ks, ks[1:]):
        if lo <= x <= hi:
            t = (x - lo) / (hi - lo)
            return math.exp(math.log(tab[lo]) * (1 - t) + math.log(tab[hi]) * t)


def peak_in(wl, a, lo, hi):
    w = [(x, v) for x, v in zip(wl, a) if lo <= x <= hi and v == v]
    return max(w, key=lambda t: t[1]) if w else (None, None)


def fit(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((p - mx) ** 2 for p in xs)
    sxy = sum((p - mx) * (q - my) for p, q in zip(xs, ys))
    m = sxy / sxx
    c = my - m * mx
    ssr = sum((q - (m * p + c)) ** 2 for p, q in zip(xs, ys))
    sst = sum((q - my) ** 2 for q in ys)
    se = math.sqrt(ssr / (n - 2) / sxx) if n > 2 and sxx else float('nan')
    return m, c, (1 - ssr / sst if sst else 0), se


def main():
    g = read_sheet()
    cols = {}
    for k, v in g[1].items():
        lab, xc = parse_label(v)
        if xc is not None:
            cols[k] = (lab, xc)
    rows = sorted(r for r in g if r >= 2)
    wl = [float(g[r]['A']) for r in rows if 'A' in g[r]]
    tab = load_act()

    rec = []
    for k, (lab, xc) in cols.items():
        ts = [float(g[r].get(k)) if g[r].get(k) else float('nan') for r in rows]
        a = [(-math.log10(t / 100.0) if t == t and t > 0 else float('nan')) for t in ts]
        q = [v for w, v in zip(wl, a) if QUIET[0] <= w <= QUIET[1] and v == v]
        a = [v - sum(q) / len(q) for v in a]
        wA, vA = peak_in(wl, a, *LO_WIN)
        wB, vB = peak_in(wl, a, *HI_WIN)
        if None in (vA, vB) or vA <= 0 or vB <= 0:
            continue
        fB = vB / (vA + vB)
        rec.append(dict(lab=lab, x=xc, wA=wA, vA=vA, wB=wB, vB=vB, fB=fB,
                        la=math.log10(a_of(tab, xc))))
    rec.sort(key=lambda r: r['x'])

    print(' C=O 두 하위밴드 — 위치는 고정, 세기비만 변한다')
    print(' %-6s %-7s %-16s %-16s %-8s %s'
          % ('시료', 'x', '밴드 A (낮음)', '밴드 B (높음)', 'f_B', 'log a'))
    print(' ' + '-' * 70)
    for r in rec:
        print(' %-6s %-7.3f %6.0f (%.4f)   %6.0f (%.4f)   %-8.4f %.2f'
              % (r['lab'], r['x'], r['wA'], r['vA'], r['wB'], r['vB'], r['fB'], r['la']))

    wAs = [r['wA'] for r in rec]; wBs = [r['wB'] for r in rec]
    print('\n 밴드 위치가 실제로 고정인가')
    print('  A: %.0f ~ %.0f cm-1 (폭 %.0f)    B: %.0f ~ %.0f cm-1 (폭 %.0f)'
          % (min(wAs), max(wAs), max(wAs) - min(wAs),
             min(wBs), max(wBs), max(wBs) - min(wBs)))
    print('  비교: 최대점 이동 폭은 15.3 cm-1 였다. 밴드 자체는 그보다 훨씬 덜 움직인다.')

    # 두 상태 평형 검정:  ln(fB/(1-fB)) = lnK + n ln a
    xs_ln = [r['la'] * math.log(10) for r in rec]          # ln a
    ys_ln = [math.log(r['fB'] / (1 - r['fB'])) for r in rec]
    n_slope, c_ln, r2_ln, se = fit(xs_ln, ys_ln)
    print('\n [핵심 검정]  ln(f_B/(1-f_B)) = ln K + n · ln a_Cl')
    print('  기울기 n = %.3f ± %.3f      R^2 = %.4f' % (n_slope, se, r2_ln))
    print('  95%% 신뢰구간 n = [%.3f, %.3f]' % (n_slope - 1.96 * se, n_slope + 1.96 * se))
    lo, hi = n_slope - 1.96 * se, n_slope + 1.96 * se
    if lo <= 1.0 <= hi:
        print('  -> n = 1 이 신뢰구간 안에 있다. **염화물 1개 결합**과 부합한다.')
    elif lo > 0:
        print('  -> n = %.2f. 염화물이 관여하되 1개는 아니다.' % n_slope)
    else:
        print('  -> n 이 0 을 포함하거나 음수. 염화물과 무관하다.')

    # 조성 축 대조 (같은 형식으로)
    m2, c2, r2_x, se2 = fit([r['x'] for r in rec], ys_ln)
    print('\n  대조: 같은 좌변을 x 로 적합하면 R^2 = %.4f (활동도 %.4f)' % (r2_x, r2_ln))

    # 금속 분율과의 대조
    print('\n [금속과의 대조]  같은 a 에서 두 화학종 분율')
    print(' %-6s %-8s %-10s %s' % ('시료', 'log a', 'f_B (시트르산)', 'f_Td (금속, 60 C)'))
    print(' ' + '-' * 52)
    for r in rec:
        ft = FTD.get(round(r['x'], 3))
        print(' %-6s %-8.2f %-10.4f %s'
              % (r['lab'], r['la'], r['fB'], ('%.4f' % ft) if ft is not None else '-'))

    assert len(rec) >= 10
    print('\n self-check OK: %d조성' % len(rec))


if __name__ == '__main__':
    main()
