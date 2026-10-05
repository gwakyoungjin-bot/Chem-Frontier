"""1;4 (NCM 있음) vs 1;4 N (NCM 없음) FT-IR 대조 — 카르복실 배위 판정.

이 한 쌍이 왜 중요한가:
  지금까지 FT-IR 로는 "COO- 가 금속에 배위하는가" 를 판정할 수 없었다. 조성끼리만
  비교한 데이터였고, DES 안 시트르산(약 3 M)이 금속(약 34 mM)보다 100배 많아서
  카르복실 신호의 대부분이 배위하지 않은 자유 시트르산이기 때문이다.
  **같은 조성에서 금속만 넣고 뺀 쌍**이 있으면 그 문제가 사라진다. 차이 스펙트럼이
  곧 금속이 만든 변화다.

판정 기준 (문헌 통용):
  COO- 비대칭(~1550-1650)과 대칭(~1380-1450) 밴드의 간격 dv
      dv 200 이상  단좌 배위
      dv 160~200   자유 카르복실레이트
      dv 160 미만  이좌/가교 배위
  그리고 COOH 의 C=O(~1700-1740)가 줄고 COO- 가 늘어나면 탈양성자화 + 배위 신호.

데이터 다루는 원칙:
  * %T 가 100 을 넘는 구간이 있다(배경 문제). 그래서 흡광도 절대값을 쓰지 않고
    **같은 시료 안의 조용한 창(1900~2100)** 을 0 으로 잡아 상대비교만 한다.
  * 두 시료의 ATR 접촉 압력이 다르면 전체가 비례로 어긋난다. 그래서 차이를 보기 전에
    **공통 밴드 하나로 스케일을 맞춘다**. 이걸 안 하면 접촉 차이를 배위로 오독한다.

실행: python dft_ni_complexes/ftir_ncm.py
입력: 바탕화면의 '1;4.csv', '1;4 N.csv'
"""
import csv
import math
import os
import sys

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

DESK = os.path.join(os.path.expanduser('~'), 'Desktop')
FILES = {'NCM 있음': '1;4.csv', 'NCM 없음': '1;4 N.csv'}
QUIET = (1900.0, 2100.0)      # 흡수가 거의 없는 창 -> 베이스라인 0 기준
SCALE_BAND = (1150.0, 1250.0)  # C-O 영역. 시트르산이 지배하므로 두 시료에 공통
REGIONS = {
    'C=O (COOH)':      (1680.0, 1780.0),
    'COO- 비대칭':      (1530.0, 1660.0),
    'COO- 대칭':        (1360.0, 1460.0),
}


def load(path):
    xs, ys = [], []
    with open(path, encoding='utf-8-sig') as f:
        for row in csv.reader(f):
            if len(row) < 2:
                continue
            try:
                xs.append(float(row[0])); ys.append(float(row[1]))
            except ValueError:
                continue                      # 헤더 두 줄
    return xs, ys


def absorbance(xs, ts):
    """%T -> 흡광도. %T<=0 은 측정 불가 구간이라 버린다."""
    return [(-math.log10(t / 100.0) if t > 0 else float('nan')) for t in ts]


def window(xs, ys, lo, hi):
    return [(x, y) for x, y in zip(xs, ys) if lo <= x <= hi and y == y]


def mean_in(xs, ys, lo, hi):
    w = window(xs, ys, lo, hi)
    return sum(y for _, y in w) / len(w) if w else float('nan')


def peak_in(xs, ys, lo, hi):
    """구간 내 흡광도 최대점 = %T 최소점."""
    w = window(xs, ys, lo, hi)
    if not w:
        return None, None
    return max(w, key=lambda p: p[1])


def main():
    data = {}
    for lab, fn in FILES.items():
        p = os.path.join(DESK, fn)
        if not os.path.exists(p):
            sys.exit('파일 없음: %s' % p)
        xs, ts = load(p)
        a = absorbance(xs, ts)
        base = mean_in(xs, a, *QUIET)
        a = [v - base for v in a]             # 조용한 창을 0 으로
        data[lab] = (xs, ts, a)
        print(' %-9s %d점  %.0f~%.0f cm-1  베이스라인 보정 %+.4f'
              % (lab, len(xs), min(xs), max(xs), -base))

    (x1, t1, a1), (x0, t0, a0) = data['NCM 있음'], data['NCM 없음']
    assert x1 == x0, '두 파일의 파수 격자가 다르다 — 보간이 필요하다'

    # ATR 접촉 차이 보정: 공통 밴드 넓이로 스케일을 맞춘다
    s1 = mean_in(x1, a1, *SCALE_BAND)
    s0 = mean_in(x0, a0, *SCALE_BAND)
    k = s0 / s1 if s1 else 1.0
    print('\n 접촉 보정: C-O 영역 흡광도 NCM있음 %.4f / 없음 %.4f -> 배율 %.3f'
          % (s1, s0, k))
    a1s = [v * k for v in a1]

    print('\n 밴드 위치 (흡광도 최대점)')
    print(' %-14s %-22s %-22s %s' % ('영역', 'NCM 있음', 'NCM 없음', '이동'))
    print(' ' + '-' * 74)
    pos = {}
    for name, (lo, hi) in REGIONS.items():
        p1, v1 = peak_in(x1, a1s, lo, hi)
        p0, v0 = peak_in(x0, a0, lo, hi)
        pos[name] = (p1, p0)
        print(' %-14s %7.0f cm-1 (A=%.3f)  %7.0f cm-1 (A=%.3f)  %+.0f cm-1'
              % (name, p1, v1, p0, v0, p1 - p0))

    # dv = 비대칭 - 대칭
    for lab, i in (('NCM 있음', 0), ('NCM 없음', 1)):
        dv = pos['COO- 비대칭'][i] - pos['COO- 대칭'][i]
        verdict = ('이좌/가교 배위' if dv < 160 else
                   '자유 카르복실레이트' if dv <= 200 else '단좌 배위')
        print(' %-9s  dv = %.0f cm-1  ->  %s' % (lab, dv, verdict))

    # 차이 스펙트럼: 금속이 만든 변화만 남는다
    print('\n 차이 스펙트럼 (NCM있음 - NCM없음), 1300~1800 cm-1')
    d = [(x, p - q) for x, p, q in zip(x1, a1s, a0) if 1300 <= x <= 1800]
    dmax = max(d, key=lambda p: p[1])
    dmin = min(d, key=lambda p: p[1])
    rms = math.sqrt(sum(v * v for _, v in d) / len(d))
    noise = math.sqrt(sum((p - q) ** 2 for x, p, q in zip(x1, a1s, a0)
                          if QUIET[0] <= x <= QUIET[1])
                      / len(window(x1, a1s, *QUIET)))
    print('  최대 증가 %+.4f @ %.0f cm-1' % (dmax[1], dmax[0]))
    print('  최대 감소 %+.4f @ %.0f cm-1' % (dmin[1], dmin[0]))
    print('  RMS %.4f   |   조용한 창의 잡음 RMS %.4f   |   신호/잡음 %.1f'
          % (rms, noise, rms / noise if noise else float('inf')))

    print('\n 판정')
    if rms < 3 * noise:
        print('  차이가 잡음 수준이다 (RMS < 3x). 금속 투입으로 카르복실 환경이')
        print('  유의하게 바뀌지 않았다 -> COO- 배위는 검출되지 않음.')
        print('  => 물 배위를 주 화학종으로 본 계산 모델과 어긋나지 않는다.')
    else:
        print('  차이가 잡음보다 크다 (%.1f배). 금속이 카르복실 환경을 바꿨다.' % (rms / noise))
        print('  위 최대 증가 위치가 1550~1650 이면 COO- 배위 신호다.')
        print('  => 계산 모델에 금속-카르복실레이트 착물을 넣어야 한다.')

    # 자체검증: 스케일 보정이 공통 밴드를 실제로 맞췄는지
    chk = mean_in(x1, a1s, *SCALE_BAND)
    assert abs(chk - s0) < 1e-9, '스케일 보정이 적용되지 않았다'
    assert noise > 0, '잡음이 0 — 두 파일이 동일한 파일일 수 있다'
    print('\n self-check OK: 스케일 보정 적용됨, 두 파일이 서로 다름')


if __name__ == '__main__':
    main()
