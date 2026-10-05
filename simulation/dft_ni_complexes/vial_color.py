"""바이알 사진에서 착색 여부를 수치화한다.

왜 필요한가: 원고가 인용하던 "0.004 / 0.221 / 0.238, 55배" 는 환산 코드도 출력도
저장소에 없어 근거를 댈 수 없었다(2026-09-30 감사). 그 수치를 재생산한다.

⚠ 한계 — 이 값은 흡광도가 아니다.
  휴대폰 카메라는 감마·화이트밸런스·자동노출이 걸려 Beer-Lambert 가 성립하지 않는다.
  게다가 진한 두 시료는 광학적으로 포화되어 있다(원고 Q16). 따라서 이 지표는
  **착색 여부의 순서 판정(ordinal)** 에만 쓰고, 농도·분율로 환산하지 않는다.
  정량은 희석 후 UV-Vis 가 필요하며 후속 과제다.

정의: A_red = -log10(R / G).  청록색 용액은 적색을 흡수하므로 R < G -> A_red > 0.
      무색이면 R ~ G -> A_red ~ 0.  G 를 기준으로 삼는 이유는 같은 사진 안의
      같은 조명을 쓰기 위해서다(외부 백색 기준이 없으므로).

실행: python dft_ni_complexes/vial_color.py
"""
import math
import os
import sys

# Windows 기본 콘솔이 cp949 라 한글/em-dash 출력에서 죽는다
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PHOTO = os.path.join(ROOT, 'photos', 'vials_260922.png')

# plot_final.py 와 동일한 crop (라벨, x, 사진상 중심 픽셀)
VIALS = [('1:19', 0.05, 1900), ('3:17', 0.15, 1050), ('1:4', 0.20, 1500),
         ('3:1', 0.75, 2350), ('19:1', 0.95, 2800)]
CROP_Y, CROP_W = (1900, 2120), 95


def a_red(px):
    """중앙값을 쓴다 — 유리 반사광(소수의 아주 밝은 픽셀)에 평균이 끌려간다."""
    r, g = float(np.median(px[:, :, 0])), float(np.median(px[:, :, 1]))
    return -math.log10(max(r, 1.0) / max(g, 1.0)), r, g


def main():
    src = Image.open(PHOTO).convert('RGB')
    print(' 시료    x      R     G     A_red')
    print(' ' + '-' * 38)
    out = {}
    for lab, x, cx in VIALS:
        px = np.asarray(src.crop((cx - CROP_W, CROP_Y[0], cx + CROP_W, CROP_Y[1])))
        a, r, g = a_red(px)
        out[lab] = a
        print(' %-6s %4.2f  %5.1f %5.1f  %+.4f' % (lab, x, r, g, a))

    pale = [out[k] for k in ('1:19', '3:17', '1:4')]
    tint = [out[k] for k in ('3:1', '19:1')]
    print('\n 무색군 3개 : 평균 %+.4f, 표준편차 %.4f' % (np.mean(pale), np.std(pale)))
    print(' 착색군 2개 : %+.4f, %+.4f' % (tint[0], tint[1]))
    print(' 두 군이 겹치는가 : %s' % ('예 — 이분 불가' if max(pale) >= min(tint) else '아니오'))
    print('\n 배수는 보고하지 않는다: 무색군이 0 근처라 분모가 불안정하다.')

    # 자체검증: 착색군이 무색군보다 확실히 크지 않으면 이 지표는 쓸 수 없다
    assert min(tint) > max(pale), 'ordinal separation failed'
    print(' self-check OK: 착색군 최소 > 무색군 최대')


if __name__ == '__main__':
    main()
