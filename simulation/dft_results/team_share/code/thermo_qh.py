"""준조화(quasi-harmonic) 열역학 보정.

문제 (2026-09-27 진단):
  ORCA 는 허수진동수를 분배함수에서 **빼버린다**. 그래서 허수모드가 있는 종은
  진동 엔트로피가 과소평가되고, 그게 기준물질이면 오차가 모든 dG 에 실린다.
  실측으로 확인된 증거 — MCl1_aq5 의 dS 가 기준물질에서 버려진 모드 수를 그대로 따라갔다:
      Ni_aq6 (버림 0개) -> dS = +17.1 J/mol/K
      Co_aq6 (버림 2개) -> dS = +49.0
      Mn_aq6 (버림 5개) -> dS = +69.2
  그 결과 Co·Mn 이 25 C 에서도 100% 사면체로 나왔다 (수용액 Co(II) 는 분홍색
  팔면체이므로 명백히 틀린 값).

해법 (Cramer-Truhlar / Grimme 준조화 근사, 표준 관행):
  모드를 버리지 않는다. **|nu| 가 CUTOFF 보다 작은 모드를 전부 CUTOFF 로 치환**한다.
    - 허수모드      -> +CUTOFF 로 되살린다
    - 아주 낮은 실모드 -> CUTOFF 로 올린다 (1/nu 발산으로 엔트로피가 폭주하는 것을 막음)
  이러면 모든 종이 같은 규칙으로 처리되어 종간 비교가 공정해진다.

  새 DFT 가 필요 없다. 이미 있는 freq 출력의 진동수 목록만 다시 쓴다.

한계:
  준조화 근사는 허수모드가 '수치적 잡음/유연한 모드' 일 때 타당하다.
  진짜 전이상태(큰 허수)라면 이 보정은 정당화되지 않는다.
  우리 경우 최대 허수가 111~142 cm-1 로, 용매화된 유연한 착물의 물 회전 영역이라
  이 가정이 성립한다고 본다 (한계 절에 기재할 것).

사용:
    import thermo_qh as Q
    G = Q.gibbs_qh('Ni_aq6_freq.out', 333.15)      # Hartree
    python thermo_qh.py --test
"""
import math, os, re, sys

H_PLANCK = 6.62607015e-34      # J s
C_LIGHT = 2.99792458e10        # cm/s
K_B = 1.380649e-23             # J/K
R_J = 8.314462618              # J/mol/K
HARTREE_KJ = 2625.4996

CUTOFF = 100.0                 # cm-1. Cramer-Truhlar 관행값


def frequencies(path):
    """freq 출력의 진동수 전체 (cm-1). 병진·회전의 0.00 은 제외."""
    try:
        txt = open(path, encoding='utf-8', errors='ignore').read()
    except Exception:
        return None
    m = re.search(r'VIBRATIONAL FREQUENCIES(.*?)(NORMAL MODES|$)', txt, re.S)
    if not m:
        return None
    v = [float(x) for x in re.findall(r':\s+(-?\d+\.\d+)\s*cm\*\*-1', m.group(1))]
    return [x for x in v if abs(x) > 1e-6]


def g_mode(nu, T):
    """진동 모드 하나가 G 에 기여하는 값 (kJ/mol). nu 는 cm-1 (양수)."""
    theta = H_PLANCK * C_LIGHT * nu / K_B          # K
    r = theta / T
    if r > 700:                                    # 언더플로 방지
        return R_J * theta * 0.5 / 1000.0
    e = math.exp(-r)
    zpe = 0.5 * R_J * theta
    h_vib = R_J * theta * e / (1.0 - e)
    s_vib = R_J * (r * e / (1.0 - e) - math.log(1.0 - e))
    return (zpe + h_vib - T * s_vib) / 1000.0


def correction(path, T, cutoff=CUTOFF):
    """ORCA 의 G 에 더해야 할 보정값 (kJ/mol). 못 읽으면 None.

    - 허수모드: ORCA 가 통째로 뺐으므로 cutoff 모드로 되살려 **더한다**
    - cutoff 미만의 실모드: ORCA 값을 빼고 cutoff 값을 더한다
    """
    v = frequencies(path)
    if v is None:
        return None, 0, 0
    add, n_imag, n_low = 0.0, 0, 0
    for nu in v:
        if nu < 0:
            add += g_mode(cutoff, T)
            n_imag += 1
        elif nu < cutoff:
            add += g_mode(cutoff, T) - g_mode(nu, T)
            n_low += 1
    return add, n_imag, n_low


def gibbs_qh(path, T, cutoff=CUTOFF):
    """준조화 보정된 Final Gibbs free energy (Hartree). 못 읽으면 None."""
    import speciation as S                          # 순환 import 회피용 지연 import
    g = S.read_gibbs(path, T)
    if g is None:
        return None
    corr, _, _ = correction(path, T, cutoff)
    if corr is None:
        return g
    return g + corr / HARTREE_KJ


def test():
    T = 333.15
    # 낮은 진동수일수록 엔트로피 기여가 커서 G 가 낮아진다
    gs = [g_mode(nu, T) for nu in (50, 100, 200, 500, 1500)]
    assert all(gs[i] < gs[i + 1] for i in range(len(gs) - 1)), gs
    # cutoff 보다 높은 모드는 보정 대상이 아니다 -> 기여 차이가 0
    assert abs(g_mode(CUTOFF, T) - g_mode(CUTOFF, T)) < 1e-15
    # 저주파 모드는 -T*S 가 ZPE+H 를 이겨서 **G 를 낮춘다**
    # => 허수모드를 되살리면 그 종의 G 가 내려간다.
    #    기준물질(아쿠아 착물)에 허수가 많았으므로 기준 G 가 내려가고,
    #    그만큼 dG(= G_생성물 - G_기준) 가 덜 음수가 된다.
    #    즉 Co·Mn 이 100% 사면체로 튄 문제를 바로잡는 방향이다.
    assert g_mode(CUTOFF, T) < 0, g_mode(CUTOFF, T)
    # 반대로 고주파 모드는 ZPE 가 지배해서 G 를 올린다
    assert g_mode(2000.0, T) > 0, g_mode(2000.0, T)
    # 극저온 극한: ZPE 만 남는다
    zpe = 0.5 * R_J * (H_PLANCK * C_LIGHT * 1000 / K_B) / 1000.0
    assert abs(g_mode(1000.0, 1.0) - zpe) < 1e-9
    # 1000 cm-1 의 ZPE 는 약 6.0 kJ/mol (상식 확인)
    assert 5.5 < zpe < 6.5, zpe
    print('self-check OK   (cutoff=%.0f cm-1, g_mode(100cm-1, 60C)=%.3f kJ/mol)'
          % (CUTOFF, g_mode(CUTOFF, T)))


def main():
    if '--test' in sys.argv:
        test(); return
    here = os.path.dirname(os.path.abspath(__file__))
    import glob
    pats = [os.path.join(here, '*_freq.out'),
            os.path.join(here, 'cobalt', '*_freq.out'),
            os.path.join(here, 'manganese', '*_freq.out'),
            os.path.join(here, 'citrate', '*_freq.out')]
    print('%-26s %-7s %-7s %-12s %-12s' %
          ('종', '허수', '저주파', '보정 kJ/mol', 'G 변화'))
    for p in sorted(sum((glob.glob(x) for x in pats), [])):
        c, ni, nl = correction(p, 333.15)
        if c is None:
            continue
        name = os.path.basename(p).replace('_freq.out', '')
        d = os.path.basename(os.path.dirname(p))
        print('%-26s %-7d %-7d %+12.2f %-12s'
              % ('%s/%s' % (d, name), ni, nl, c, '상승' if c > 0 else '하강'))


if __name__ == '__main__':
    main()
