"""14조성 FT-IR 의 C=O 파수 추세 — 내 COSMO-RS 예측을 실측으로 검증한다.

내가 세운 예측 (ftir_predict.py):
  Cl- 가 COOH 의 O-H 에 수소결합하면 C=O 가 약해져 **낮은 파수로** 간다.
  문헌 근거: 순수 CA 1736 -> ChCl:CA DES 1722 (ChCl 을 넣으면 내려감).
  따라서 **x(ChCl) 이 커질수록 C=O 파수가 낮아져야** 한다.
  뒤집어 말하면 **CA 가 많아질수록 C=O 파수가 높아져야** 한다.

민이형 관찰:
  "CA 비율 증가하면서 15 cm-1 red shift" = CA 가 많아질수록 **낮아진다**.
  내 예측과 방향이 반대다. 둘 중 하나가 틀렸으므로 데이터로 가른다.

측정 방법:
  C=O 영역(1650~1800)에서 흡광도 최대점을 찾고, 그 주변 3점 포물선 맞춤으로
  격자(1 cm-1)보다 세밀한 위치를 낸다. 봉우리 꼭짓점만 보면 잡음에 1 cm-1 씩 흔들린다.

실행: python dft_ni_complexes/ftir_trend.py
입력: 바탕화면 '새 스프레드시트 문서.xlsx'
"""
import math
import os
import re
import sys
import zipfile

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

XLSX = os.path.join(os.path.expanduser('~'), 'Desktop', '새 스프레드시트 문서.xlsx')
CO_WIN = (1650.0, 1800.0)
QUIET = (1900.0, 2100.0)

# 헤더의 소수 표기 -> (라벨, x_ChCl).  '1.19' 는 1;19 을 뜻한다.
def parse_label(v):
    """헤더가 1;19 를 숫자 1.19 로 저장해 놓아서, 2;3 이 2.2999999999999998 로
    들어온다. 소수부를 원래 자릿수만큼 반올림해 복원한다."""
    if v.strip().lower() in ('1.4n', '1,4n'):
        return '1;4N', None
    # 엑셀이 '1;19' 를 숫자 1.19 로 저장해 2;3 이 2.2999999999999998 로 들어온다.
    # 소수 둘째자리로 반올림한 뒤 뒤 0 을 떼면 원래 표기가 복원된다.
    #   2.2999999999999998 -> '2.30' -> '2.3' -> 3
    #   1.19               -> '1.19'         -> 19
    s = ('%.2f' % float(v)).rstrip('0')
    whole, frac = s.split('.')
    a, b = int(whole), int(frac)
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
            if not ref or not v:
                continue
            val = v.group(1)
            if re.search(r't="s"', c):
                val = ss[int(val)]
            grid.setdefault(int(rn), {})[ref.group(1)] = val
    return grid


def peak(xs, ys, lo, hi):
    """흡광도 최대점 + 3점 포물선 보간. 격자보다 세밀한 위치를 낸다."""
    idx = [i for i, x in enumerate(xs) if lo <= x <= hi]
    i = max(idx, key=lambda j: ys[j])
    if i in (0, len(xs) - 1):
        return xs[i]
    y0, y1, y2 = ys[i - 1], ys[i], ys[i + 1]
    den = y0 - 2 * y1 + y2
    d = 0.5 * (y0 - y2) / den if den else 0.0
    return xs[i] + d * (xs[i + 1] - xs[i])


def main():
    grid = read_sheet()
    hdr = grid[1]
    cols = {}
    for col, v in hdr.items():
        try:
            cols[col] = parse_label(v)
        except Exception:
            pass

    rows = sorted(r for r in grid if r >= 2)
    xs = [float(grid[r]['A']) for r in rows if 'A' in grid[r]]

    res = []
    for col, (lab, xc) in cols.items():
        ts = []
        for r in rows:
            v = grid[r].get(col)
            ts.append(float(v) if v not in (None, '') else float('nan'))
        if len(ts) != len(xs):
            continue
        a = [(-math.log10(t / 100.0) if t == t and t > 0 else float('nan')) for t in ts]
        q = [v for x, v in zip(xs, a) if QUIET[0] <= x <= QUIET[1] and v == v]
        a = [v - sum(q) / len(q) for v in a]
        res.append((lab, xc, peak(xs, a, *CO_WIN)))

    known = sorted([r for r in res if r[1] is not None], key=lambda r: r[1])
    extra = [r for r in res if r[1] is None]

    print(' C=O 파수 실측  (x(ChCl) 오름차순 = CA 감소 방향)')
    print(' %-7s %-9s %-11s %s' % ('시료', 'x(ChCl)', 'C=O /cm-1', '앞 조성 대비'))
    print(' ' + '-' * 52)
    prev = None
    for lab, xc, p in known:
        d = '' if prev is None else '%+.1f' % (p - prev)
        print(' %-7s %-9.3f %-11.1f %s' % (lab, xc, p, d))
        prev = p
    for lab, _, p in extra:
        print(' %-7s %-9s %-11.1f  (NCM 없음, 1;4 와 비교)' % (lab, '-', p))

    lo, hi = known[0], known[-1]
    span = hi[2] - lo[2]
    print('\n 양 끝단: %s(CA 최다) %.1f  ->  %s(ChCl 최다) %.1f   총 %+.1f cm-1'
          % (lo[0], lo[2], hi[0], hi[2], span))

    # 단조성 검사: 몇 쌍이 순서를 지키는가
    ps = [p for _, _, p in known]
    dec = sum(1 for a, b in zip(ps, ps[1:]) if b < a)
    inc = sum(1 for a, b in zip(ps, ps[1:]) if b > a)
    print(' 인접쌍 %d개 중 감소 %d, 증가 %d' % (len(ps) - 1, dec, inc))

    print('\n 판정')
    if span < -3:
        print('  x(ChCl) 이 커질수록 C=O 가 내려간다 = CA 가 많을수록 높다.')
        print('  -> 내 예측과 **일치**. Cl- 수소결합 해석이 실측으로 지지된다.')
    elif span > 3:
        print('  x(ChCl) 이 커질수록 C=O 가 올라간다 = CA 가 많을수록 낮다.')
        print('  -> 내 예측과 **반대**. 민이형 관찰(CA 증가 -> red shift)이 맞다.')
        print('     Cl- 수소결합만으로는 설명이 안 된다. 다른 원인을 찾아야 한다.')
    else:
        print('  이동이 %.1f cm-1 로 작다. 추세를 주장하기 어렵다.' % span)

    assert len(known) >= 10, '조성이 10개 미만이면 추세를 말할 수 없다'
    print('\n self-check OK: %d 조성, 파수 격자 %.0f~%.0f cm-1'
          % (len(known), min(xs), max(xs)))


if __name__ == '__main__':
    main()
