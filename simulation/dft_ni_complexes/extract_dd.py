"""*_dd.out (CASSCF/NEVPT2) -> 화학종별 d-d 밴드표 + 조성비별 예측 스펙트럼.

사용법 (계산 끝난 뒤 같은 폴더에서):
    python extract_dd.py            # 표 + dd_bands.csv + dd_spectra.png
    python extract_dd.py --test     # 파서 자체검증만 (ORCA 출력 없어도 실행됨)

파싱 원리: ORCA 흡수스펙트럼 표의 어느 버전이든 (에너지 cm-1, 파장 nm)이
나란히 오고 nm = 1e7/cm-1 을 만족한다. 이 관계로 컬럼을 찾으므로
ORCA 버전별 컬럼 변경에 영향을 받지 않는다.
"""
import glob, os, re, sys

# 소수점을 필수로 하고 앞뒤에 영숫자가 오지 못하게 막는다.
# 예전 정규식 `\d+\.?\d*(?:[eEdD]...)?` 은 ORCA 출력의 해시 문자열 `487d211` 을
# "487" + 지수 "d211" 로 읽어 float() 에서 죽었다. 우리가 뽑는 값(에너지·파장·진동자세기)은
# 전부 소수점을 가지므로 소수점 필수로 두면 그런 토큰을 애초에 건드리지 않는다.
FLOAT = re.compile(r'(?<![A-Za-z0-9_.])[-+]?\d+\.\d+(?:[eEdD][-+]?\d+)?(?![A-Za-z0-9_])')
VIS = (300., 1300.)          # 관심 영역(nm). Ni(II)의 세 밴드가 다 들어감
FWHM = 2000.                 # 가우시안 브로드닝 폭 (cm-1). d-d 밴드 전형값

# 진동자세기 f -> 몰흡광계수 eps_max (M-1 cm-1).  f = 4.32e-9 * integral(eps d_nu)
# 가우시안 밴드 가정: integral = 1.0645 * eps_max * FWHM
F2EPS = 1.0 / (4.32e-9 * 1.0645 * FWHM)


def transitions(text):
    """(파장nm, 진동자세기) 목록. 흡수스펙트럼 표에서만 뽑는다."""
    out = []
    for line in text.splitlines():
        f = [float(x.replace('D', 'E').replace('d', 'E')) for x in FLOAT.findall(line)]
        for i in range(len(f) - 2):
            e, lam = f[i], f[i + 1]
            if e > 1000 and lam > 0 and abs(1e7 / e - lam) < 1.0:
                out.append((lam, f[i + 2]))
                break
    # ORCA는 같은 표를 여러 번(CASSCF/NEVPT2) 찍으므로 중복 제거
    seen, uniq = set(), []
    for lam, fosc in out:
        k = round(lam, 1)
        if k not in seen:
            seen.add(k); uniq.append((lam, fosc))
    return uniq


def eps_max(fosc):
    """진동자세기를 몰흡광계수 최대값으로 환산 (가우시안 밴드 가정)."""
    return fosc * F2EPS


def lft_params(text):
    """AILFT가 찍은 10Dq/Racah 줄을 그대로 회수 (형식이 버전마다 달라 원문 보존)."""
    return [l.strip() for l in text.splitlines()
            if re.search(r'10\s*Dq|Racah', l, re.I)]


def spectrum(bands, grid):
    """가우시안 브로드닝. bands=[(nm, fosc)] -> grid(nm) 위의 상대흡광도."""
    y = [0.0] * len(grid)
    for lam, fosc in bands:
        if fosc <= 0:
            continue
        c0 = 1e7 / lam
        for i, g in enumerate(grid):
            y[i] += fosc * 2.71828 ** (-((1e7 / g - c0) / (FWHM / 1.6651)) ** 2)
    return y


def test():
    """ORCA 6 / ORCA 5 두 형식 + 잡음줄 섞어서 파서 검증."""
    sample = """
-----------------------------------------------------------------------------
         ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
-----------------------------------------------------------------------------
State   Energy  Wavelength   fosc         T2        TX        TY        TZ
        (cm-1)    (nm)                  (au**2)    (au)      (au)      (au)
-----------------------------------------------------------------------------
   1   8500.0   1176.5   0.000120   0.00123  0.01  0.02  0.03
   2  13793.1    725.0   0.000250   0.00251  0.01  0.02  0.03
  0-3A -> 5-3A  25316.5    395.0   0.000410   0.00412  0.01  0.02  0.03
Total run time: 1234.5 sec
  10Dq =   8500.0 cm-1
  Racah B =  900.0 cm-1
"""
    t = transitions(sample)
    assert len(t) == 3, t
    assert [round(l) for l, _ in t] == [1176, 725, 395], t
    assert abs(t[2][1] - 0.000410) < 1e-9, t
    # 'Total run time: 1234.5 sec' 같은 줄이 전이로 잡히면 안 됨
    assert all(300 < l < 1300 for l, _ in t), t
    assert len(lft_params(sample)) == 2
    # f=0.000410, FWHM=2000 -> eps ~ 45 M-1cm-1 규모여야 함 (자릿수 확인)
    e = eps_max(0.000410)
    assert 10 < e < 200, e
    assert abs(eps_max(0.0) ) < 1e-12

    g = [300. + i for i in range(1000)]
    y = spectrum(t, g)
    peak = g[max(range(len(g)), key=lambda i: y[i])]
    assert abs(peak - 395) < 15, peak          # 가장 센 전이 자리에 극대
    assert spectrum([(500., 0.0)], g) == [0.0] * len(g)   # fosc=0은 기여 없음
    print('self-check OK')


def main():
    outs = sorted(glob.glob('*_dd.out'))
    if not outs:
        sys.exit('*_dd.out 이 없습니다. 계산을 먼저 돌리세요 (--test 로 파서만 검증 가능).')

    data, rows = {}, []
    for f in outs:
        name = os.path.basename(f)[:-7]
        text = open(f, encoding='utf-8', errors='replace').read()
        bands = [(l, o) for l, o in transitions(text) if VIS[0] <= l <= VIS[1]]
        data[name] = bands
        ftot = sum(o for _, o in bands)
        print(f"\n=== {name} ===  전이 {len(bands)}개  |  Σf={ftot:.3e}"
              f"  eps_max≈{eps_max(max((o for _, o in bands), default=0.)):.1f} M-1cm-1")
        for lam, fosc in sorted(bands):
            print(f"   {lam:8.1f} nm  ({1e7/lam:8.0f} cm-1)  f={fosc:.2e}  eps≈{eps_max(fosc):7.1f}")
            rows.append(f"{name},{lam:.2f},{1e7/lam:.1f},{fosc:.6e},{eps_max(fosc):.3f}")
        for l in lft_params(text):
            print('   LFT|', l)

    open('dd_bands.csv', 'w').write(
        "species,wavelength_nm,energy_cm1,fosc,eps_max_M-1cm-1\n" + "\n".join(rows) + "\n")
    print(f"\n-> dd_bands.csv ({len(rows)}행)")

    # ---- 인텐시티 요약: 조성비 설명의 핵심 지표 ----
    print("\n" + "=" * 70)
    print("화학종별 인텐시티 (가시영역 350-800 nm 합)")
    print(f"{'화학종':<14}{'Σf':>12}{'eps_max':>11}{'Ni_aq6 대비':>13}  기하")
    vis = lambda bs: sum(o for l, o in bs if 350 <= l <= 800)
    ref = vis(data.get('Ni_aq6', []))
    for name, bands in data.items():
        v = vis(bands)
        rel = f"{v/ref:9.1f}x" if ref > 1e-12 else "    (기준≈0)"
        geo = 'Td (중심대칭 없음)' if ('NiCl4' in name or 'NiCl3' in name) else 'Oh (중심대칭)'
        print(f"{name:<14}{v:12.3e}{eps_max(v):11.1f}{rel:>13}  {geo}")
    print("""
[반드시 읽을 것] 팔면체(Oh) 화학종의 계산 세기는 숫자 그대로 쓰지 마세요.
  Oh는 중심대칭이라 d-d 전이가 Laporte 금지 -> 정적 계산에서 f≈0이 '정답'입니다.
  실험에서 보이는 Oh의 약한 흡수(eps 1~10)는 진동-전자 결합(vibronic coupling)
  때문인데 이 계산에는 그 항이 없습니다.
  => Td/Oh 세기비를 인용하지 말 것. "Td가 자릿수로 더 세다"는 대칭성(Laporte 규칙)
     논증으로 말하고, 계산은 Td 정량값 제공 + Oh≈0 확인 용도로만 쓰세요.""")

    try:
        import matplotlib; matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return
    grid = [VIS[0] + i for i in range(int(VIS[1] - VIS[0]))]
    plt.figure(figsize=(9, 4.5))
    for name, bands in data.items():
        plt.plot(grid, spectrum(bands, grid), lw=1.2, label=name)
    plt.axvline(395, color='k', ls=':', lw=.8)
    plt.axvline(725, color='k', ls=':', lw=.8)
    plt.xlabel('nm'); plt.ylabel('relative absorbance (fosc, broadened)')
    plt.title('predicted d-d spectra per Ni(II) species (CASSCF/NEVPT2)')
    plt.legend(fontsize=8); plt.tight_layout(); plt.savefig('dd_spectra.png', dpi=130)
    print('-> dd_spectra.png')


if __name__ == '__main__':
    test() if '--test' in sys.argv else main()
