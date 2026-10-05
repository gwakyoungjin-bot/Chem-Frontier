"""전체 UV-Vis (250~800 nm, 13조성) 분석 — 345 nm 밴드의 정체를 가른다.

왜 이걸 가려야 하나:
  345 nm 가 **진짜 전자전이**면 그건 LMCT 일 가능성이 높고, LMCT 는 리간드->금속
  전하이동이라 **시트르산이 배위해야 생긴다**. 그러면 mixed ligand 계산이 정당해지고
  TD-DFT 로 검증할 수 있다.
  반대로 그게 **산란/베이스라인 기울기**면 설명할 대상이 아니다.

가르는 방법:
  * 진짜 밴드 = 극대점이 있다. 2차 미분이 음수로 뚜렷하다.
  * 산란/오프셋 = 단조 감소(또는 평탄)하고 극대점이 없다. 장파장에서 0 으로 안 간다.
  * 그래서 (1) 원 스펙트럼의 극대, (2) 2차 미분, (3) 장파장 꼬리를 같이 본다.

추가 점검:
  바이알 사진에서는 1;4 가 **무색**, 3;1·19;1 이 청록이었다. 이 데이터가 그와
  맞는지 본다. 어긋나면 둘 중 하나가 다른 시료이거나 아티팩트다.

실행: python dft_ni_complexes/uvvis_full.py
입력: 바탕화면 '최종.csv'
"""
import csv
import os
import sys

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

PATH = os.path.join(os.path.expanduser('~'), 'Desktop', '최종.csv')
X = {'1;19': .05, '1;9': .10, '3;17': .15, '1;4': .20, '1;3': .25, '1;2': .333,
     '2;3': .40, '1;1': .50, '3;2': .60, '2;1': .667, '3;1': .75, '17;3': .85,
     '19;1': .95}
VIAL = {'1;19': '무색', '3;17': '무색', '1;4': '무색', '3;1': '청록', '19;1': '청록'}


def load():
    with open(PATH, encoding='utf-8-sig') as f:
        rows = list(csv.reader(f))
    hdr = [h.strip() for h in rows[0]]
    labs = [h for h in hdr[1:] if h]
    data = {l: [] for l in labs}
    wl = []
    for r in rows[1:]:
        if not r or not r[0].strip():
            continue
        wl.append(float(r[0]))
        for i, l in enumerate(labs, start=1):
            v = r[i].strip()
            data[l].append(float(v) if v else float('nan'))
    return wl, data


def smooth(y, w=9):
    """이동평균. 2차 미분을 내려면 먼저 잡음을 눌러야 한다."""
    out = []
    for i in range(len(y)):
        a, b = max(0, i - w // 2), min(len(y), i + w // 2 + 1)
        s = [v for v in y[a:b] if v == v]
        out.append(sum(s) / len(s) if s else float('nan'))
    return out


def maxima(wl, y, lo, hi, min_prom=0.004):
    """구간 내 국소 극대. 양옆보다 min_prom 이상 높아야 인정한다."""
    out = []
    for i in range(2, len(y) - 2):
        if not (lo <= wl[i] <= hi):
            continue
        if y[i] >= y[i - 1] and y[i] >= y[i + 1]:
            left = min(y[max(0, i - 25):i] or [y[i]])
            right = min(y[i + 1:i + 26] or [y[i]])
            prom = y[i] - max(left, right)
            if prom >= min_prom:
                out.append((wl[i], y[i], prom))
    # 가까운 것끼리 묶어 가장 높은 것만
    out.sort(key=lambda t: -t[1])
    keep = []
    for w_, v, p in out:
        if all(abs(w_ - k[0]) > 20 for k in keep):
            keep.append((w_, v, p))
    return sorted(keep)


def at(wl, y, target):
    i = min(range(len(wl)), key=lambda j: abs(wl[j] - target))
    return y[i]


def main():
    wl, data = load()
    print(' 파장 %.0f~%.0f nm, %d점, %d조성\n' % (min(wl), max(wl), len(wl), len(data)))

    labs = sorted(data, key=lambda l: X.get(l, 9))
    print(' 주요 파장 흡광도')
    print(' %-6s %-7s %7s %7s %7s %7s %7s %7s  %s'
          % ('시료', 'x(ChCl)', '250', '345', '400', '500', '700', '800', '바이알'))
    print(' ' + '-' * 78)
    for l in labs:
        y = data[l]
        print(' %-6s %-7.2f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f  %s'
              % (l, X.get(l, float('nan')), at(wl, y, 250), at(wl, y, 345),
                 at(wl, y, 400), at(wl, y, 500), at(wl, y, 700), at(wl, y, 800),
                 VIAL.get(l, '')))

    print('\n 실제 극대점 (돌출도 0.004 이상, 300~800 nm)')
    for l in labs:
        s = smooth(data[l])
        pk = maxima(wl, s, 300, 800)
        txt = ', '.join('%.0f nm (A=%.3f, 돌출 %.3f)' % p for p in pk) if pk else '없음 — 단조'
        print('  %-6s %s' % (l, txt))

    print('\n 345 nm 가 밴드인가 기울기인가 (평활 후 2차 미분; 음수 클수록 봉우리)')
    print(' %-6s %10s %10s %12s' % ('시료', 'A(345)', "d2A/dλ2", '판정'))
    print(' ' + '-' * 44)
    for l in labs:
        s = smooth(data[l], 15)
        i = min(range(len(wl)), key=lambda j: abs(wl[j] - 345))
        h = 8
        d2 = (s[i - h] - 2 * s[i] + s[i + h]) / (h * (wl[1] - wl[0])) ** 2
        v = '봉우리' if d2 < -2e-6 else ('평탄/기울기' if d2 < 2e-6 else '골')
        print(' %-6s %10.3f %10.2e %12s' % (l, at(wl, data[l], 345), d2, v))

    print('\n 장파장 꼬리 — 진짜 흡수는 800 nm 에서 0 으로 수렴한다')
    print(' %-6s %8s %8s %8s  %s' % ('시료', 'A(750)', 'A(800)', '기울기', '해석'))
    print(' ' + '-' * 52)
    for l in labs:
        y = data[l]
        a750, a800 = at(wl, y, 750), at(wl, y, 800)
        slope = a800 - a750
        note = ('꼬리 큼 — 산란/오프셋 의심' if a800 > 0.05 else
                '0 근처 — 정상')
        print(' %-6s %8.3f %8.3f %8.3f  %s' % (l, a750, a800, slope, note))

    # 사진과의 대조
    print('\n 바이알 색과 대조 (청록이면 600~700 nm 흡수가 있어야 한다)')
    for l in ('1;19', '3;17', '1;4', '3;1', '19;1'):
        if l not in data:
            continue
        a = at(wl, data[l], 650)
        print('  %-6s 사진 %-4s  A(650) = %+.3f  %s'
              % (l, VIAL[l], a,
                 '<- 어긋남' if (VIAL[l] == '청록') != (a > 0.05) else ''))

    assert len(data) >= 10, '조성이 너무 적다'
    print('\n self-check OK: %d조성 로드' % len(data))


if __name__ == '__main__':
    main()
