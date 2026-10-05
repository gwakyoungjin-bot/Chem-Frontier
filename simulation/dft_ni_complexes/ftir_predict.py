"""FT-IR 14조성에 대한 예측표 — 계산과 실측을 붙여보기 위한 것.

목적:
  민이형 FT-IR 이 14조성을 다 찍었다. 우리 COSMO-RS 는 조성마다 염화물 활동도를
  내놓는다. 둘을 붙이면 **우리 모델이 틀렸는지 실측으로 판정할 수 있다.**

예측의 근거:
  ChCl 을 시트르산에 넣으면 Cl- 가 COOH 의 O-H 에 수소결합한다(Cl-...H-O).
  그러면 C=O 가 약해져 **더 낮은 파수로 이동**한다. 문헌값: 순수 CA 1736 -> DES 1722.
  ChCl 이 많아질수록 이 이동이 커져야 한다.

무엇이 걸려 있나:
  * 이동이 x(ChCl) 에 **단조**면 -> 우리 모델의 방향이 맞다
  * 단조가 아니거나 반대면 -> 우리 활동도 계산이 뭔가 놓치고 있다
  * 1;4 와 1;4N(NCM 첨가) 이 **같으면** -> 금속은 카르복실기에 거의 안 붙는다
    (= 우리 물 배위 모델이 맞다).  다르면 -> 우리 모델을 고쳐야 한다.

주의: a_ChCl 은 ChCl 을 중성 이온쌍으로 근사한 **대리변수**다. 절대값이 아니라
      조성 간 순서와 자릿수만 쓴다.

실행: python dft_ni_complexes/ftir_predict.py
입력: dft_results/ternary_grid.csv (서버 COSMO-RS 산출)
"""
import csv
import os
import sys

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

HERE = os.path.dirname(os.path.abspath(__file__))
GRID = os.path.join(os.path.dirname(HERE), 'dft_results', 'ternary_grid.csv')

# 라벨은 ChCl ; 시트르산 순서 (2026-09-30 확정)
SAMPLES = ['1;19', '3;17', '1;9', '1;4', '1;3', '1;2', '2;3', '1;1',
           '3;2', '2;1', '3;1', '17;3', '19;1']
UVVIS = {'1;19', '3;17', '1;4', '3;1', '19;1'}     # UV-Vis 도 찍은 5조성
T_PICK = 298.15                                    # FT-IR 은 상온 측정


def x_of(lab):
    a, b = (float(v) for v in lab.split(';'))
    return a / (a + b)


def load(x_h2o):
    """(x_ChCl_dry -> a_ChCl) 표. 해당 물함량·온도 행만."""
    out = {}
    with open(GRID, encoding='utf-8') as f:
        for r in csv.DictReader(f):
            if abs(float(r['x_H2O']) - x_h2o) < 1e-9 and \
               abs(float(r['T_K']) - T_PICK) < 1e-6:
                out[round(float(r['x_ChCl_dry']), 6)] = float(r['a_ChCl'])
    return out


def interp(tab, x):
    """격자 사이는 log 선형 보간. a 가 자릿수로 변하므로 선형보간은 틀린다."""
    import math
    ks = sorted(tab)
    if x <= ks[0]:
        return tab[ks[0]]
    if x >= ks[-1]:
        return tab[ks[-1]]
    for lo, hi in zip(ks, ks[1:]):
        if lo <= x <= hi:
            t = (x - lo) / (hi - lo)
            return math.exp(math.log(tab[lo]) * (1 - t) + math.log(tab[hi]) * t)


def main():
    dry, wet = load(0.0), load(0.8)
    if not dry and not wet:
        sys.exit('ternary_grid.csv 에서 %g K 행을 못 찾음' % T_PICK)

    rows = []
    for lab in SAMPLES:
        x = x_of(lab)
        rows.append((lab, x, interp(dry, x) if dry else None,
                     interp(wet, x) if wet else None))
    rows.sort(key=lambda r: r[1])

    print(' FT-IR 조성별 예측  (25 °C, ChCl;시트르산 순서)')
    print(' ' + '-' * 72)
    print(' %-7s %-8s %-13s %-13s %-8s %s'
          % ('시료', 'x(ChCl)', 'a_ChCl 무수', 'a_ChCl 물0.8', 'UV-Vis', 'C=O 예측'))
    print(' ' + '-' * 72)
    n = len(rows)
    for i, (lab, x, ad, aw) in enumerate(rows):
        # 파수 이동은 순위로만 말한다. 절대 파수를 예측할 근거가 우리에겐 없다.
        rank = '가장 높음' if i == 0 else ('가장 낮음' if i == n - 1 else '%d위' % (i + 1))
        print(' %-7s %-8.3f %-13.3e %-13.3e %-8s %s'
              % (lab, x, ad or 0, aw or 0, '○' if lab in UVVIS else '', rank))

    print()
    print(' 예측 1 — C=O 파수는 위에서 아래로 **단조 감소**해야 한다.')
    print('          (Cl- 가 COOH 에 수소결합할수록 C=O 가 약해진다)')
    print(' 예측 2 — 1;4 와 1;4N(NCM 첨가) 의 카르복실 영역이 **거의 같아야** 한다.')
    print('          금속 34 mM vs 시트르산 3 M 이라 붙어도 1%% 수준이기 때문.')
    print('          눈에 띄게 다르면 우리 물 배위 모델을 고쳐야 한다.')
    print()
    print(' 활동도가 %.0f 배 (무수 기준 양 끝단) 벌어지는데, 파수 이동은 그만큼'
          % (rows[-1][2] / rows[0][2]))
    print(' 벌어지지 않을 것이다 — 수소결합은 포화되기 때문. 순서만 본다.')

    # 자체검증: 활동도가 x 에 단조 증가해야 한다. 아니면 예측 자체가 성립 안 한다.
    a = [r[2] for r in rows]
    assert all(p < q for p, q in zip(a, a[1:])), 'a_ChCl 이 x 에 단조가 아니다'
    print('\n self-check OK: a_ChCl 이 x(ChCl) 에 단조 증가 (%d 조성)' % n)


if __name__ == '__main__':
    main()
