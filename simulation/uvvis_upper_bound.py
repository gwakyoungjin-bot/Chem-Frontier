"""UV-Vis 검출한계 -> 사면체 Ni 화학종 농도의 상한.

밴드가 안 보인 것도 정량 정보다. "안 보였다"는 곧 "이 농도 이상은 아니다"이다.
Beer-Lambert:  A = eps * c * l   ->   c < A_detect / (eps * l)

    python uvvis_upper_bound.py          # 표
    python uvvis_upper_bound.py --test   # 자체검증
"""
import glob, os, sys

L_CM = 1.0                  # 광로길이(cm). ※ 민이 형 확인 필요
SIGMA = 3.0                 # 검출한계 기준 (3 sigma)
# 문헌 전형값. DFT 끝나면 dd_bands.csv 의 eps_max 로 교체할 것
EPS = {'사면체 [NiCl4]2-': 150.0, '팔면체 [Ni(H2O)6]2+': 3.0}
FLAT = (420, 580)           # 밴드가 없는 평탄구간 — 여기서 노이즈를 잰다


def load(f):
    x, on = {}, False
    for ln in open(f, encoding='utf-8', errors='replace'):
        ln = ln.strip()
        if ln == 'XYDATA':
            on = True; continue
        if not on:
            continue
        p = ln.split(',')
        if len(p) != 2:
            break
        try:
            x[float(p[0])] = float(p[1])
        except ValueError:
            break
    return x


def noise_rms(spec, lo=FLAT[0], hi=FLAT[1]):
    """이웃 파장 간 차이의 RMS. 완만한 배경은 상쇄되고 순수 잡음만 남는다."""
    d = [spec[w] - spec[w + 1] for w in range(lo, hi) if w in spec and w + 1 in spec]
    return (sum(v * v for v in d) / len(d)) ** 0.5


def upper_bound(a_detect, eps, l=L_CM):
    """Beer-Lambert 역산: 이 흡광도 이하라면 농도는 얼마 미만인가 (M)."""
    return a_detect / (eps * l)


def test():
    # 잡음 없는 직선 -> RMS 0
    flat = {w: 0.5 for w in range(400, 600)}
    assert noise_rms(flat) < 1e-12
    # 기울기만 있는 배경 -> 이웃차가 일정하므로 RMS = 기울기 (배경이 아니라 잡음만 재는지 확인)
    ramp = {w: 0.001 * w for w in range(400, 600)}
    assert abs(noise_rms(ramp) - 0.001) < 1e-9, noise_rms(ramp)
    # Beer-Lambert 왕복
    assert abs(upper_bound(0.04, 150.0, 1.0) - 0.04 / 150) < 1e-15
    assert upper_bound(0.04, 3.0) > upper_bound(0.04, 150.0)   # 약한 흡수체일수록 상한 느슨
    print('self-check OK')


def main():
    files = sorted(glob.glob('1;4.csv_260920/*.csv'))
    if not files:
        sys.exit('1;4.csv_260920/*.csv 를 찾을 수 없습니다.')
    D = {os.path.basename(f)[:-4].replace('_', ':'): load(f) for f in files}

    rms = {k: noise_rms(s) for k, s in D.items()}
    allr = sum(([s[w] - s[w + 1] for w in range(*FLAT)] for s in D.values()), [])
    RMS = (sum(v * v for v in allr) / len(allr)) ** 0.5
    A_DET = SIGMA * RMS

    print(f"[1] 노이즈 ({FLAT[0]}-{FLAT[1]} nm 평탄구간, 이웃점 차의 RMS)")
    for k, v in sorted(rms.items(), key=lambda kv: kv[1]):
        print(f"    {k:>6}  RMS = {v:.4f} A")
    print(f"    전체 RMS = {RMS:.4f} A")
    print(f"    -> 검출한계 ({SIGMA:.0f}σ) = {A_DET:.3f} A\n")

    print(f"[2] 이 검출한계에서 허용되는 최대 농도  (광로 {L_CM:.0f} cm)")
    print(f"    {'화학종':<22}{'eps (M-1cm-1)':>15}{'농도 상한':>14}")
    ub = {}
    for name, e in EPS.items():
        ub[name] = upper_bound(A_DET, e)
        print(f"    {name:<22}{e:>15.0f}{ub[name]*1000:>11.2f} mM")
    print("    ※ eps는 문헌 전형값. DFT 완료 후 dd_bands.csv 의 eps_max 로 교체할 것")

    td = ub['사면체 [NiCl4]2-']
    print(f"\n[3] 핵심 결과")
    print(f"    밴드가 검출되지 않았다는 것은,")
    print(f"    사면체 Ni 화학종 농도가 {td*1000:.2f} mM 미만이라는 뜻이다.")
    print(f"    (팔면체는 흡수가 약해 {ub['팔면체 [Ni(H2O)6]2+']*1000:.1f} mM 까지 있어도 안 보인다)\n")

    print(f"[4] 총 Ni 농도를 알면 '사면체 비율의 상한'이 나온다  <- ICP 데이터 필요")
    print(f"    {'총 [Ni]':>12}{'사면체 최대비율':>16}   해석")
    for c in (0.1e-3, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1):
        frac = td / c
        if frac >= 1:
            note = "제약 없음 — 전부 사면체여도 안 보임"
        elif frac > 0.1:
            note = "약한 제약"
        else:
            note = "강한 제약 — 사면체 우세 가설과 충돌"
        print(f"    {c*1000:>9.1f} mM{min(frac,1.0)*100:>14.1f}%   {note}")

    print("""
[5] 이 결과를 어떻게 쓰나
    - 실험 데이터를 버리지 않는다. "검출 안 됨"을 정량적 상한으로 바꾼 것이다.
    - 총 [Ni]가 약 0.3 mM 이하면: 안 보이는 게 당연 -> 계산과 모순 없음
    - 총 [Ni]가 약 3 mM 이상이면: 사면체가 10% 넘게 있었다면 보였어야 함
      -> 염 과잉 조성에서도 사면체가 우세하지 않다는 뜻 -> 계산 가정을 수정해야 함
    => 어느 쪽이든 결론이 나온다. 필요한 건 ICP 한 줄.

[한계]
    - 광로길이 1 cm 가정. 다르면 상한이 그만큼 비례해서 바뀐다.
    - eps는 문헌 전형값. DFT 후 우리 계산값으로 교체할 것.
    - 받침대(offset)는 산란이 아니라 기기/블랭크 오프셋으로 보인다
      (짧은 파장에서 커지지 않고 평탄하거나 오히려 장파장에서 큼).
      따라서 시료 간 받침대 높이 차이는 화학 정보가 아니다.""")


if __name__ == '__main__':
    test() if '--test' in sys.argv else main()
