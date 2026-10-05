"""DES 의 pH 와 시트르산 탈양성자화 비율 추정.

측정 pH 가 없어도 추정 가능한 이유: 이 계에 염기가 없다. 시트르산이 유일한
산이고 양성자 받개는 물뿐이므로 pH 는 시트르산 농도로 결정된다.
ponytail: 1단계 해리(pKa1)만 본다. pKa2=4.76 은 pH~1 에서 무의미(<0.1%).

한계 (이게 이 스크립트의 결론을 뒤집을 수 있는 유일한 지점):
  중성이 95~99% 라는 건 "양이 많다" 는 뜻이지 "금속에 붙는 건 중성이다" 가 아니다.
  COO- 는 5% 뿐이어도 결합력이 세면 금속 배위권을 다 차지할 수 있다.
  뒤집히는 문턱 = 농도 불리 8.3 + 기존 격차 4.8 = 약 13 kJ/mol (333 K).
  그런데 COO- 착물을 DFT 로 직접 재려면 양성자 1개의 기준 자유에너지가 필요하고
  그 오차가 10~20 kJ/mol 이라 재려는 효과보다 크다 -> 계산으로는 판정 불가.
  대신 Ni/Co 앵커가 독립적으로 65.1 / 67.7 로 일치한 것이 경험적 방어다.
"""
import math

PKA1 = 3.13                    # 시트르산 1차 해리 (수용액, 25 C)
KA1 = 10 ** (-PKA1)
X_H2O = 0.80                   # 실험 물 몰분율 (ternary_water.py 와 동일)
MW = dict(h2o=18.015, ca=192.12, chcl=139.62)
RHO = dict(h2o=1.00, ca=1.665, chcl=1.10)      # g/mL

def volume_L(n):
    """부피가법 근사. 조성 간 '상대' 비교가 목적이라 초과부피는 무시."""
    return sum(n[k] * MW[k] / RHO[k] for k in n) / 1000.0

def h_plus(C):
    """Ka = [H+]^2/(C-[H+]) 를 정확히 푼다 (약산 근사는 이 농도에서 깨진다)."""
    return (-KA1 + math.sqrt(KA1 * KA1 + 4 * KA1 * C)) / 2.0

print(' 조성    x(ChCl)   [CA]/M    pH    H3Cit(중성) %')
print(' ' + '-' * 47)
for lab, x in (('1:19', 0.05), ('3:17', 0.15), ('1:4', 0.20),
               ('1:2 관행', 0.333), ('3:1', 0.75), ('19:1', 0.95)):
    n = dict(h2o=X_H2O, chcl=(1 - X_H2O) * x, ca=(1 - X_H2O) * (1 - x))
    C = n['ca'] / volume_L(n)
    h = h_plus(C)
    ph = -math.log10(h)
    f_neutral = 1.0 / (1.0 + 10 ** (ph - PKA1))
    print(' %-8s %5.2f   %6.2f   %5.2f      %5.1f' % (lab, x, C, ph, 100 * f_neutral))

# 자체검증: 묽어지면 pH 가 올라가고 중성분율이 떨어져야 한다
a, b = h_plus(3.0), h_plus(0.01)
assert a > b, 'more acid must give more H+'
assert 1.0 / (1 + 10 ** (-math.log10(a) - (-PKA1))) > 1.0 / (1 + 10 ** (-math.log10(b) - (-PKA1)))
print('\n self-check OK')
