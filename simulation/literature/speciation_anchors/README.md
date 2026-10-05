# 화학종 앵커 문헌 (2026-09-27 확보)

`dft_ni_complexes/` 의 화학종 지도를 **실측에 고정(anchor)** 하기 위해 모은 논문들.
전부 DOI 검증 완료(OpenAlex), 철회 없음.

## 왜 앵커가 필요한가

우리 계산은 Ni(II) 클로로/아쿠아 착물의 ΔG 를 DFT 로 내지만, 전하가 +2 → −2 로
바뀌는 반응이라 **암시적 용매 오차가 수십 kJ/mol** 이다. 실제로 보정값이
+67 kJ/mol per Cl⁻ 까지 나온다. 그래서 생 ΔG 를 그대로 못 쓰고,
**"어떤 조건에서 화학종 분율이 얼마"** 라는 실측 한 점에 맞춰야 한다.

## 앵커가 갖춰야 할 조건 (중요도순)

| 축 | 우리 조건 | 중요도 |
|---|---|---|
| 물 함량 | **몰분율 0.77~0.82** (DES 3000 mg + 물 1270 µL) | ★★★ |
| 금속 | Ni (그리고 Co, Mn) | ★★★ |
| 정량성 | 화학종 분율 숫자 (색 변화 서술로는 부족) | ★★★ |
| HBD | 시트르산 | ★★★ ← 아래 참조 |
| 기법 | XAS/EXAFS > UV-Vis | ★★ |

**HBD 가 ★★★ 인 이유 (2026-09-27 계산으로 확인)**: 같은 조성·같은 물함량에서
HBD 만 바꿔 a_ChCl 을 COSMO-RS 로 비교하니 **시트르산이 EG 보다 48.7배 낮았다**
(물 몰분율 0.045 기준). 물 효과(18~168배)와 같은 자릿수다.
→ "HBD 는 달라도 된다" 는 가정은 **기각**. 재현: `python dft_ni_complexes/hbd_effect.py`

다만 우리 모델은 **조성이 아니라 활동도에 앵커를 걸므로**, 앵커 논문의 계
(그 HBD, 그 물함량)에서 a_ChCl 을 계산해 쓰면 HBD·물 불일치가 동시에 처리된다.

---

## 파일별 용도

### ⭐ Busato2022_Ni-urea-MDES_water-scan_XAS.pdf
`10.1021/acs.inorgchem.2c00864` · Inorg. Chem. 2022 · Busato, Tofoni, Mannucci et al.
"On the Role of Water in the Formation of a DES Based on NiCl₂·6H₂O and Urea"

**Ni 앵커 1순위 후보.** NiCl₂·6H₂O + urea (1:3.5) 를 물/MDES 몰비 **W 를 바꿔가며**
MD·ab initio·UV-Vis·NIR·SAXS/WAXS·**XAS** 로 추적. 금속이 니켈이고 물이 변수라
우리 최대 불일치(물 함량)를 정면으로 메운다. Mannucci 2026 과 **같은 그룹**.

뽑을 것: W 별 Ni 화학종 분율, Ni–Cl 배위수, 전환이 일어나는 W 지점.

### ⭐ Iname2026_Co-Ni-speciation-DES_colorimetric.pdf
`10.1016/j.molliq.2026.129750` · J. Mol. Liq. 2026 · Iname, Vitry, Tapsoba
"Chemical speciation of Co²⁺ and Ni²⁺ ions in DES and colorimetric determination of
cobalt from lithium-ion battery samples"

**주제가 가장 많이 겹친다.** Ni 와 Co 를 한 논문에서 다루므로 앵커 두 개를 동시에
얻을 수 있다. 리튬이온전지 시료라 응용도 같다.
OpenAlex 에 초록이 없어 제목 외에는 원문을 봐야만 판단 가능.

뽑을 것: 두 금속의 화학종 분율, 사용한 DES 조성, 물 함량.

### Mannucci2026_Co-ChCl_water-dilution_XAS-MCR.pdf
`10.1021/acs.inorgchem.6c01344` · Inorg. Chem. 2026

**Co 앵커 (확보 완료).** ChCl/CoCl₂·6H₂O (1:2) 를 물로 W=0~50 희석하며 XAS+MCR 로 정량:
- W=0 (물 몰분율 **0.80** — 우리와 동일): Td [CoCl₄]²⁻ **40%** / Oh 클로로 40% / 완전수화 20%
- W=50 (물 0.95): 사면체 **검출 안 됨**, 완전수화 80%

→ 우리 계산(Co 100% 사면체)이 **틀렸다는 직접 증거**이자 고칠 앵커 두 점.

### Hartley2025_Ni-ChCl-EG_thermochromism_EXAFS.pdf
`10.1021/acs.jpcc.5c05771` · J. Phys. Chem. C 2025

**현재 쓰고 있는 Ni 앵커 1.** ChCl:EG **1:2** (x_ChCl = 0.333 — 우리 가정과 일치 ✅)
에서 Oh→Td 전환이 **90~100 °C**.
⚠ 물 함량이 **몰분율 0.045** (NiCl₂·6H₂O 0.1 M 에서 온 것뿐) — 우리(0.80)의 **1/18**.
⚠ HBD 가 EG — 시트르산 대비 a_ChCl 이 48.7배 다름.
→ 앵커 활동도를 **이 논문 자신의 계**에서 계산해야 한다. (a_ChCl = 1.83e−1 @95 °C)

### Hammond2023_Ce-ChCl-urea-water_neutron.pdf
`10.1021/acs.inorgchem.3c02205` · Inorg. Chem. 2023

금속이 Ce 라 직접 앵커는 아니지만, ChCl/urea/water 를 **w = 2, 5, 10** 으로
중성자회절+EPSR 로 측정해 **"물이 리간드를 밀어내는 정도"를 정량**한다.
"질량·부피 분율이 아니라 **몰비**가 결정변수" 라는 우리 논증의 직접 근거.

### Ma2018_water-effect-IL-DES_review.pdf
`10.1039/c8cs00325d` · Chem. Soc. Rev. 2018 · 595회 인용

IL/DES–물 계 6종의 물함량 전 영역 물성을 종합한 리뷰. 배경 근거용.

### Minerals2025_aqueous-ChCl-DES_CuCo-leaching.pdf
`10.3390/min15080815` · Minerals 2025

**"Aqueous ChCl-based DES"** 로 Cu–Co 광석을 침출. 화학종 논문은 아니지만
**물이 많은 DES 로 침출한다**는 운전 조건 선례. 우리 조건(물 80%)의 정당화에 쓸 수 있다.

---

## 검색으로 확인된 공백 (포스터에 쓸 것)

**ChCl : 카르복실산 계의 금속 화학종 연구를 못 찾았다.** 독립 검색 3회 전부 빈손:
- `carboxylic acid DES metal coordination speciation citric malic oxalic`
- `oxalic acid ChCl DES metal complex EXAFS coordination`
- `levulinic/malonic/glycolic ChCl DES transition metal speciation UV-Vis`

DES 화학종 연구는 거의 전부 릴라인(ChCl:urea)·에타린(ChCl:EG) 이고,
카르복실산 DES 는 **침출 수율 논문만** 있다.

→ 여기에 위 HBD 계산(시트르산이 a_ChCl 을 50배 억누름)을 붙이면:
**"표준 DES 의 화학종 지식을 카르복실산 DES 에 그대로 옮길 수 없는데,
그 계의 선행 화학종 연구가 없다"** — 갭 주장이 강화된다.

⚠ "못 찾았다" 는 OpenAlex 기준이며 존재하지 않는다는 뜻이 아니다.
Crossref/Scopus 교차검색은 미실시.
