# 조성비별 UV-Vis 설명용 DFT 설계 + 문헌근거 (2026-09-23)

> 이 문서는 **완성된 보고서 문장이 아니라 사실관계·수치·근거 모음**입니다.
> 보고서/발표 문장은 직접 작성하는 것을 전제로 합니다.

---

## 0. 전제 확인 필요 (민이 형에게)

| 항목 | 현재 가정 | 확인 필요 이유 |
|---|---|---|
| `1;4` 등의 표기 순서 | **ChCl : 유기산(CA)** 순 | 순서가 반대면 아래 해석의 Cl-rich/Cl-poor 방향이 통째로 뒤집힘 |
| 셀 광로길이 | 1 cm 가정 | 몰흡광계수 환산에 필요 |
| 레퍼런스(블랭크) | 물? 같은 조성 DES? | 같은 조성 DES여야 배경이 상쇄됨 |
| 침출 금속 | Ni 중심 가정 | Co가 섞이면 밴드 해석이 달라짐 (Co(II)는 Ni보다 훨씬 강한 색) |

가정대로라면 5개 시료는 **x(ChCl) 축 위의 5점**입니다:

| 시료 | 몰비 | x(ChCl) | 성격 |
|---|---|---|---|
| 1:19 | 1:19 | **0.05** | 산 과잉 극단 |
| 3:17 | 3:17 | **0.15** | 산 과잉 |
| 1:4 | 1:4 | **0.20** | 산 과잉 |
| 3:1 | 3:1 | **0.75** | 염 과잉 |
| 19:1 | 19:1 | **0.95** | 염 과잉 극단 |

> 샘플링 공백: **x = 0.20 ~ 0.75 사이가 통째로 비어 있음.** 문헌 표준 DES 조성인
> 1:2 (x = 0.333)와, 우리 COSMO-RS가 예측한 gE 최소점(x ≈ 0.35~0.40, `HANDOFF_NOTES.md` F절)이
> 둘 다 이 공백 안에 들어갑니다. 다음 실험 배치에서 x = 0.33 한 점은 꼭 추가할 것.

---

## 1. 문헌조사 — 침출효율을 지배하는 물성

### 1-1. 연구질문 해석

"DES/유기산 혼합 침출계에서 **금속 침출효율을 지배하는 물성**은 무엇이고, 그중
**조성비(몰비)에 의존하는 것**은 무엇이며, 그것이 **금속 화학종(speciation) → UV-Vis 스펙트럼**으로
어떻게 이어지는가." 대상 분야: 용매야금/이온야금(ionometallurgy), 2015년 이후 리뷰 및 1차 논문 우선.

### 1-2. 검색 전략

- 도구 순서: 로컬 `pyalex`(OpenAlex) → DOI 검증. 외부 폴백 미사용(전부 OpenAlex에서 해결).
- 질의 12건: `deep eutectic solvent metal leaching efficiency viscosity molar ratio`,
  `water content deep eutectic solvent leaching lithium ion battery cathode metal`,
  `deep eutectic solvent hydrogen bond donor acidity metal oxide dissolution mechanism`,
  `nickel/cobalt chloro complex speciation UV-Vis`,
  `speciation chloride water activity metal complex ionic liquid EXAFS UV-Vis deep eutectic`,
  `CASSCF NEVPT2 ab initio ligand field theory d-d transition`,
  `choline chloride urea ethaline viscosity water content conductivity molar ratio dependence`,
  `citric acid choline chloride deep eutectic solvent leaching NMC`,
  `COSMO-RS prediction deep eutectic solvent property screening activity coefficient` 외.
- 필터: `from_publication_date = 2015-01-01`, relevance 정렬, 상위 8건씩 검토.

### 1-3. 참고문헌 검증표

전부 DOI 조회로 검증, **철회(retraction) 없음(is_retracted = False)**, 출처는 모두 OpenAlex.

| # | 저자/제목(축약) | 연도 | DOI | 라벨 | 초록 |
|---|---|---|---|---|---|
| R1 | Chemical speciation of Co²⁺ and Ni²⁺ ions in deep eutectic solvents and highly efficient colorimetric determination | 2026 | 10.1016/j.molliq.2026.129750 | VERIFIED | OpenAlex에 초록 없음 |
| R2 | Solvometallurgical recovery of cobalt from LIB cathode materials using deep-eutectic solvents (ChCl–citric acid) | 2020 | 10.1039/d0gc00940g | VERIFIED | 확보 |
| R3 | Status and advances of deep eutectic solvents for metal separation and recovery | 2022 | 10.1039/d1gc03851f | VERIFIED | 확보(짧음) |
| R4 | Separation of nickel from cobalt and manganese in LIBs using deep eutectic solvents (oxalic acid:ChCl) | 2022 | 10.1039/d2gc00606e | VERIFIED | 확보 |
| R5 | Decomposition of Deep Eutectic Solvent Aids Metals Extraction in LIB Recycling | 2022 | 10.1002/cssc.202200966 | VERIFIED | 확보 |
| R6 | Model for Metal Extraction from Chloride Media with Basic Extractants: A Coordination Chemistry Approach | 2019 | 10.1021/acs.inorgchem.9b01782 | VERIFIED | 확보 |
| R7 | Assessment of TD-DFT and LF-DFT for study of d–d transitions in first row transition metal hexaaqua complexes | 2015 | 10.1063/1.4922111 | VERIFIED | 확보 |
| R8 | Overview of acidic deep eutectic solvents on synthesis, properties and applications | 2019 | 10.1016/j.gee.2019.03.002 | VERIFIED | 확보 |
| R9 | Synthesis and Dissolution of Metal Oxides in Ionic Liquids and Deep Eutectic Solvents | 2019 | 10.3390/molecules25010078 | VERIFIED | 확보 |
| R10 | A Comprehensive Study of Density, Viscosity, and Electrical Conductivity of Choline Halide-Based Eutectic Solvents | 2024 | 10.1021/acs.jced.4c00218 | VERIFIED | 확보 |
| R11 | Meta-analysis of viscosity of aqueous deep eutectic solvents and their components | 2020 | 10.1038/s41598-020-78101-y | VERIFIED | 확보 |
| R12 | Liquid Structure and Transport Properties of the Deep Eutectic Solvent Ethaline | 2020 | 10.1021/acs.jpcb.0c04058 | VERIFIED | 확보 |
| R13 | Review on Hydrometallurgical Recovery of Metals with Deep Eutectic Solvents | 2020 | 10.3390/suschem1030016 | VERIFIED | 확보 |
| R14 | On the Role of Water in the Formation of a Deep Eutectic Solvent Based on NiCl₂·6H₂O and Urea | 2022 | 10.1021/acs.inorgchem.2c00864 | VERIFIED | 확보 |
| R15 | Trace Water Changes Metal Ion Speciation in Deep Eutectic Solvents: Ce³⁺ Solvation and Nanoscale Water Clusters | 2023 | 10.1021/acs.inorgchem.3c02205 | VERIFIED | 확보 |
| R16 | The peculiar effect of water on ionic liquids and deep eutectic solvents | 2018 | 10.1039/c8cs00325d | VERIFIED | 확보 |
| R17 | ChCl–ethylene glycol based DESs as lixiviants for cobalt recovery from LIBs | 2022 | 10.1039/d2gc02075k | VERIFIED | 확보 |
| R18 | Everything You Wanted to Know about Deep Eutectic Solvents but Were Afraid to Be Told | 2023 | 10.1146/annurev-chembioeng-101121-085323 | VERIFIED | 확보 |
| R19 | Modeling the Physicochemical Properties of Natural Deep Eutectic Solvents | 2020 | 10.1002/cssc.202000286 | VERIFIED | 확보 |
| R20 | Are There Magic Compositions in Deep Eutectic Solvents? Effects of Composition and Water Content (AIMD, ChCl:EG) | 2020 | 10.1021/acs.jpcb.0c04844 | VERIFIED | 확보 (※ `literature/acs.jpcb.0c04844.pdf`로 이미 원문 보유) |
| R21 | Efficient Extraction of Li, Co, Ni from NMC Cathodes with ChCl–pyrogallol DES | 2025 | 10.3390/recycling10030088 | VERIFIED | 확보 |
| R22 | Sustainable leaching of critical metals from LIB black mass using citric acid + choline chloride DES | 2025 | 10.1007/s43621-025-02214-5 | VERIFIED | 확보 |
| R23 | Role of Oxidants in Metal Extraction from Sulfide Minerals in a Deep Eutectic Solvent | 2024 | 10.1021/acsomega.4c01052 | VERIFIED | 확보 |

### 1-4. 근거 종합 — 침출효율 지배 물성 5종

아래는 VERIFIED 문헌에서 확보한 초록 범위에서 **직접 지지되는 내용만** 정리한 것입니다.
(각 항목의 정량적 상관계수까지는 초록만으로 확정 불가 → 원문 확보 필요 항목은 §1-5에 표시)

**(P1) 점도 / 물질전달 — 조성비 의존 강함, 우리 COSMO-RS로 간접 접근 가능**
- R11(Sci Rep 2020): DES의 높은 상온 점도가 응용을 제한하며, **물 함량과 온도로 미세조정**해야 한다고 명시. 4종 DES + 성분 단독에 대한 메타분석.
- R10(JCED 2024): 콜린할라이드 + EG/Gly 계에서 **물 첨가가 밀도·점도·전기전도도에 미치는 영향**을 조성 범위 전체에서 실측. 과잉몰부피(excess molar volume)까지 제시.
- R12(JPCB 2020): 에타린(ChCl:EG 1:2)의 밀도·확산계수·점도·구조인자를 MD로 재현. 용매화 환경이 **수소결합 종류에 따라 동적으로 바뀜**.
- R13(Sustain Chem 2020): DES 침출의 **동역학(kinetics)** 이 남은 과제라고 리뷰 차원에서 지목.
- → 우리 파이프라인 연결점: 점도 자체는 opencosmorspy가 직접 못 내지만, `gE_over_RT`·`pm_E_hb`·σ-moment는 문헌 QSPR의 점도 예측 입력(`HANDOFF_NOTES.md` C-1의 Mohan 2024 R²=0.99)과 동일한 형식.

**(P2) HBD 산도(Brønsted/Lewis) — 금속산화물 용해의 1차 동력**
- R8(Green Energy Environ 2019): 산성 DES를 Brønsted/Lewis로 분류하고, **어는점·산도·밀도·점도·전도도**를 조성 설계 변수로 다룸. 용해·추출·금속전착 응용을 산성도 축으로 정리.
- R9(Molecules 2019): IL/DES에서의 금속산화물 **용해** 를 별도 리뷰 주제로 다룸(저온 조건).
- R2(Green Chem 2020): **ChCl–citric acid DES** 로 LiCoO₂에서 Co를 회수 — 우리 계와 가장 가까운 선행연구.
- R22(Discover Sustain 2025): citric acid + ChCl DES로 블랙매스 침출, **DES 비율(ratio)** 을 온도·고액비와 함께 최적화 대상으로 명시.
- R4(Green Chem 2022): oxalic acid:ChCl로 LiNMC에서 Co·Mn만 선택 침출, Ni은 잔사에 농축 → **산 종류가 선택성을 만든다**는 직접 증거.

**(P3) 물 함량 / 물 활동도 — 침출효율과 화학종을 동시에 흔드는 변수**
- R16(Chem Soc Rev 2018): IL/DES에 **미량의 물조차 물성을 크게 바꾼다**; 물은 불순물이자 점도 저감·가격 저감용 첨가제. 6개 대표계의 실험 관측을 미시 메커니즘에 매핑.
- R15(Inorg Chem 2023): ChCl/urea/water 계에서 **희석의 몰비(molar hydration ratio w)** 가 Ce³⁺ 배위껍질 조성을 좌우 — **질량·부피 분율이 아니라 몰비가 결정 변수**라고 명시. Cl⁻와 H₂O가 다른 리간드를 밀어냄.
- R14(Inorg Chem 2022): NiCl₂·6H₂O + urea MDES는 **무수염으로는 아예 공융이 안 만들어짐**. 물/MDES 몰비 W를 바꿔가며 MD·ab initio·**UV-Vis·NIR**·SAXS/WAXS·XAS로 Ni²⁺ 클러스터 구조 추적 → 우리와 같은 분광+계산 교차검증 구조.
- R20(JPCB 2020): ChCl:EG를 1:1 / 1:2 / 1:2:1(+water)로 AIMD 비교 — **"매직 조성"이 있는가**를 조성 축에서 직접 물은 논문. 우리 핵심가설("관행 몰비 = 공융점일 뿐")과 정면으로 맞물림.

**(P4) 염화물 활동도 / 금속 화학종(speciation) — UV-Vis 피크의 직접 원인**
- R1(J Mol Liq 2026): **DES 내 Co²⁺·Ni²⁺의 화학종을 규정하고 이를 이용한 비색(colorimetric) 정량** — 제목 자체가 "DES에서 Ni/Co 화학종 → 색"입니다. **우리 주제와 가장 직접 겹치는 논문이므로 원문 확보 최우선.**
- R6(Inorg Chem 2019): 염화물 매질에서의 금속 추출을 **배위화학 모델**로 재해석. 수화로 안정화가 덜 된(=전하밀도가 낮은) 화학종이 더 잘 추출된다는 가설. → Ni(H₂O)₆²⁺(전하밀도 높음) vs NiCl₄²⁻(낮음)의 대비가 추출능 차이를 만든다는 논리의 문헌 근거.
- R15: Ce³⁺ 배위껍질이 **주로 염화물**, 다음이 물 — DES에서 Cl⁻가 1차 리간드라는 직접 관측.
- → 우리 계 해석: x(ChCl)이 커질수록 Cl⁻ 활동도↑, 물 활동도↓ → Ni(II)가 **팔면체 아쿠아/클로로 → 사면체 클로로**로 이동. 이 전이가 UV-Vis에서 가장 극적인 신호(§3).

**(P5) DES 자체의 열/화학 안정성 — 고온 침출에서 "물성"이 아니라 "변질"이 효율을 만듦**
- R5(ChemSusChem 2022): ChCl–EG의 **분해생성물이 오히려 금속산화물 용해를 돕는다**. 즉 관측된 고효율이 DES 본래 물성이 아닐 수 있음 → 재사용성/경제성 훼손.
- R17(Green Chem 2022): 180 °C 고온 ChCl:EG 침출 연구들이 **ChCl:EG의 제한된 열안정성을 무시**하고 있다고 지적, 구조분석으로 확인.
- R18(Annu Rev 2023): "DES"라는 용어와 그 지속가능성·안정성·독성 주장에 쌓인 오해를 정리 — **eutectic과 deep eutectic을 열역학적으로 구분**할 것을 권고.
- → 우리 조건은 60 °C이므로 분해 위험은 낮지만, **보고서 한계 절에 반드시 명시할 것**. 특히 19:1(x=0.95) 같은 극단 조성은 애초에 균질한 DES가 맞는지 확인 필요(R18의 지적).

### 1-5. 한계

- **초록 수준 종합입니다.** R2·R3은 OpenAlex 초록이 한 문장짜리 그래픽 초록이라 정량값을 못 얻었고, R1은 초록이 아예 없습니다. §1-4의 각 주장은 "그 논문이 그 주제를 다룬다"는 수준까지만 검증됐고, **구체적 수치(예: 점도 X cP에서 침출률 Y%)는 원문 확보 후에만 인용 가능**합니다.
- OpenAlex 단일 소스, 2015년 이후, 영어권 색인 한정. Crossref/Scopus 교차검색 미실시.
- "침출효율 vs 단일물성"의 정량적 상관을 **직접 보고한** 논문은 이번 검색에서 못 찾았습니다. 대부분 온도·시간·고액비·몰비를 조작변수로 두고 효율을 최적화할 뿐, 물성을 매개변수로 놓지 않습니다 → **이게 우리 갭분석의 핵심 주장과 일치**(관행적 최적화 vs 물성 기반 설명).
- 검색에 잡히지 않은 것: Abbott 그룹의 초기 이온야금 원논문(2006–2014), 즉 2015 필터로 잘렸습니다. 필요하면 별도 조회할 것.

### 1-6. 제외/미검증 후보

| 후보 | 사유 |
|---|---|
| "Mechanisms of Photoredox Catalysis Featuring Nickel–Bipyridine Complexes" (10.1021/acscatal.4c02036) | VERIFIED이나 주제 불일치(광촉매). 검색 노이즈 |
| Abbott 2006~2014 원논문군 | 연도 필터로 미검색. NOT_SEARCHED |
| 점도-침출효율 정량 상관 논문 | 이번 검색으로는 NOT_FOUND. 존재하지 않는다는 뜻은 아님 |

### 1-7. 더 볼 만한 인접 개념

1. **Ionometallurgy** — DES/IL 기반 야금의 상위 용어. R3·R13이 진입점.
2. **Water activity (a_w) vs water content** — R15·R16이 몰비 기준이 옳다고 지적. 우리 COSMO-RS가 a_w를 직접 산출 가능.
3. **AILFT (ab initio ligand field theory)** — CASSCF에서 10Dq·Racah B를 뽑아 d-d 밴드를 "리간드장 언어"로 설명. 이미 `_dd.inp`의 `actorbs dorbs`가 이걸 켭니다.
4. **Type V DES / eutectic vs deep eutectic 구분** — R18. 극단 조성(19:1, 1:19)이 DES인지 그냥 혼합물인지 방어 논리에 필요.
5. **Excess molar volume / excess property** — R10. 조성비를 연속변수로 다룰 때 "비이상성"의 실험적 지표. 우리 `gE_over_RT`의 실측 대응물.
6. **Colorimetric metal determination in DES** — R1. UV-Vis를 정량 도구로 쓰는 방법론.

---

## 2. DFT 설계

### 2-1. 설명해야 할 인과사슬

```
조성비 x(ChCl)
   → [Cl⁻] 활동도 ↑ / [H₂O] 활동도 ↓      (P3, P4 / COSMO-RS로 계산 가능)
   → Ni(II) 화학종 분포 이동               (Oh 아쿠아 → Oh 클로로 → Td 클로로)
   → 리간드장 분열 10Dq 및 대칭성 변화      (CASSCF/NEVPT2 + AILFT)
   → d-d 전이 에너지·세기 변화             (← UV-Vis 피크)
```

**핵심: 몰비 축에는 DFT를 돌리지 않습니다.** 조성비마다 계산을 따로 도는 게 아니라,
**화학종마다 한 번씩** 계산하고, 조성비는 "그 화학종들의 혼합 비율"로 설명합니다
(`HANDOFF_NOTES.md` C-2와 같은 논리 — 계산비용은 몰비 축에서 거의 0).

### 2-2. 계산 화학종 5종 (`build_inputs.py`)

| 화학종 | 전하 | 기하 | 원자수 | def2-TZVP 기저함수* | Cl-활동도 구간 |
|---|---|---|---|---|---|
| `Ni_aq6` [Ni(H₂O)₆]²⁺ | +2 | Oh | 19 | ~303 | Cl 없음 (극단) |
| `NiCl1_aq5` [NiCl(H₂O)₅]⁺ | +1 | Oh | 17 | ~297 | 낮음 |
| `NiCl2_aq4` [NiCl₂(H₂O)₄] | 0 | Oh (trans) | 15 | ~291 | 중간 |
| `NiCl3_aq1` [NiCl₃(H₂O)]⁻ | −1 | Td | 7 | ~199 | 높음 |
| `NiCl4` [NiCl₄]²⁻ | −2 | Td | 5 | ~193 | 포화 (극단) |

\* Ni 45 + Cl 37 + O 31 + H 6 으로 합산(def2-TZVP 축약 기준). 계산비용 스케일링의 기준값.

전부 Ni(II) d⁸ → S=1 (mult 3). CAS(8,5)에서 삼중항 10개(³F+³P), 단일항 15개(¹D+¹G+¹S) —
`nroots 10,15`는 **d⁸ 다중항 전체를 정확히 덮습니다.** 원 인풋 설계가 이 점은 정확합니다.

### 2-3. 2단계 방법론

| 단계 | 방법 | 이유 |
|---|---|---|
| 1. 기하최적화 | `r2SCAN-3c Opt CPCM(SMD/water)` | 복합법이라 BSSE 보정 내장, 전이금속 착물 기하에 검증된 저비용 조합 |
| 2. d-d 전이 | `CASSCF(8,5)/SC-NEVPT2 def2-TZVP CPCM(water)` + AILFT | **TD-DFT는 전이금속 d-d에서 오차가 큼**. R7(J Chem Phys 2015)이 1주기 전이금속 헥사아쿠아 착물에서 TD-DFT가 d², d⁴, 저스핀 d⁶에서만 만족스럽다고 보고 → Ni(II) d⁸은 그 범주 밖 |

R7이 우리 방법 선택의 직접적 문헌 근거입니다. 발표에서 "왜 TD-DFT를 안 썼나" 질문이 나오면 이 논문.

### 2-4. 이번에 적용한 인풋 수정 (`build_inputs.py`)

| 변경 | 전 | 후 | 이유 |
|---|---|---|---|
| 1단계 진동계산 | `Opt Freq` | `Opt` | **UV-Vis 설명에 Freq 불필요.** 수치 Hessian이면 19원자 = 115회 gradient → 런타임 약 2배. 화학종 분포 ΔG까지 갈 때 되살릴 것(코드에 주석으로 명시) |
| RI 근사 | `RIJCOSX` + `def2/JK` | `RIJK` + `def2/JK` | 조합 불일치. RIJCOSX는 `def2/J`와, RIJK는 `def2/JK`와 짝. CASSCF/NEVPT2 + 이 크기(≤300 BF)면 RIJK가 정확도·안정성 모두 유리 |
| 메모리 | `maxcore 4000`, `--mem=8G` | `maxcore 8000`, `--mem=16G` | 25 roots NEVPT2가 4 GB/코어로는 디스크 스래싱 위험. RAM 251 GB이므로 여유 |
| `species.txt` 순서 | 생성 순 | **우선순위 순** | `--array=1-2`, `1-3`이 곧 §3의 우선순위 상위 N개가 되도록 |
| `run_array.sh` | — | 단계별 벽시계 → `timing.log` | 첫 잡으로 나머지를 외삽 (§5) |

**아직 안 고친 것 / 알려진 한계**
- CAS(8,5)는 d 오비탈만 → **double-shell effect 누락**으로 d-d 에너지가 보통 10~20% 과대평가됩니다.
  밴드 *순서*와 *상대 이동*은 맞지만 절대 파장은 오차가 있습니다.
  보정 옵션: CAS(8,10)(d + d') 또는 실험 한 점에 대한 스케일링. **시간 남으면 NiCl₄에만 CAS(8,10)을 추가로 돌려 보정계수를 뽑는 게 가장 싸게 먹힙니다.**
- 1단계는 SMD, 2단계는 일반 CPCM으로 용매모델이 불일치(CASSCF에 SMD 미정의). 정상적인 타협이지만 한계 절에 기재.
- 암시적 용매만 사용 → DES의 실제 2차 배위권(콜린 양이온, 시트르산 카르복실기)이 없습니다. R14·R15가 보여주듯 DES 내 실제 배위환경은 물+Cl⁻ 혼합이므로, **"물 연속체 근사"라는 점을 명시**해야 합니다.

---

## 3. 5개 중 무엇을 돌릴 것인가 (우선순위)

판단 기준: **UV-Vis 스펙트럼에서 서로 구별되는가**(설명력) ÷ **계산비용**(§5).

핵심 사실: 팔면체 Ni(II) 아쿠아/클로로 화학종들(`Ni_aq6`, `NiCl1_aq5`, `NiCl2_aq4`)은
리간드장 세기 차이가 작아 **스펙트럼이 서로 매우 비슷**합니다(약한 밴드, 완만한 적색이동).
반면 **팔면체 → 사면체 전환은 색이 통째로 바뀌는 사건**입니다(사면체 클로로니켈레이트는
몰흡광계수가 두 자릿수 이상 크고 가시영역에 강한 밴드). 따라서 정보량은 **양 끝점에 몰려 있습니다.**

| 순위 | 화학종 | 왜 이걸 먼저 | 설명 가능해지는 시료 |
|---|---|---|---|
| **1** | `NiCl4` | **가장 쌈**(5원자, 193 BF, 강직한 Td → opt 수렴 빠름). 파이프라인 검증 + §5 시간 보정의 기준점. Cl-포화 극단 | 19:1 (x=0.95), 3:1 (x=0.75) |
| **2** | `Ni_aq6` | Cl-없음 극단. **문헌 실측 스펙트럼이 가장 잘 알려진 계 → 우리 방법의 정확도 검증 앵커**(R7이 바로 이 계열을 다룸) | 1:19 (0.05), 3:17 (0.15), 1:4 (0.20) |
| **3** | `NiCl2_aq4` | Oh↔Td 교차 직전의 중간체. 중간 조성(x≈0.33 등 향후 시료)을 설명할 때 필요 | 미래 x=0.33 시료 |
| 4 | `NiCl3_aq1` | 스펙트럼이 `NiCl4`와 거의 겹침(둘 다 Td 클로로). 추가 설명력 낮음 | — |
| 5 | `NiCl1_aq5` | 스펙트럼이 `Ni_aq6`와 거의 겹침. 추가 설명력 가장 낮음 | — |

**2개만 돌린다면 → `NiCl4` + `Ni_aq6`.**
이 둘이 전체 계열의 **양 끝을 괄호로 묶습니다.** 측정된 5개 조성 전부를
"두 극단 화학종의 혼합비"로 설명할 수 있고, 이게 조성비→색 논리의 최소 완결형입니다.

**3개라면 → 위 둘 + `NiCl2_aq4`.** 중간 조성 예측력이 생기고, 비어 있는 x=0.20~0.75 구간에
대한 예측을 발표에서 제시할 수 있습니다(= 다음 실험 설계 근거).

> 실행 명령: `sbatch --array=1-2 run_array.sh` (또는 `1-3`, `1-5`). `species.txt`가 이미 이 순서입니다.

---

## 4. 실측 데이터의 한계 (보고서 "한계" 절에 넣을 사실)

재측정은 불가한 상황이므로, **가진 데이터로 말할 수 있는 것과 없는 것의 경계**만 못박아 둡니다.
DFT 설계는 이것과 무관하게 §2·§3·§5·§6대로 진행합니다.

`uvvis_sanity.py` / `ni_band2.py` 실행 결과 (수치는 재현 가능):

| 사실 | 수치 |
|---|---|
| 측정 범위 | 200–800 nm. **Ni(II) 최강 밴드 ν₁(~1176 nm)은 범위 밖** |
| 검출기 포화(A = 7.0) 구간 | 208–217 nm — 이 아래는 무효 |
| 엑셀의 "218–221 nm 피크" | 이웃 2–3 nm 평균과 +0.33 ~ +0.41 차이 → 단일점 스파이크 |
| 300–460 nm 요철 극대 | 5개 시료 **전부** 304 / 345 / 363 / 398 / 416 / 438 nm (±6 nm), 평균 간격 27 nm |
| 요철 진폭 vs 시료간 베이스라인 차 | 0.184 A vs **0.224 A** (오차가 더 큼) |
| 340 nm 신호와 600–700 nm 오프셋 상관 | Pearson **r = 0.894** |

→ **말할 수 있는 것**: 398 nm 부근에 요철이 존재하며, 그 위치는 [Ni(H₂O)₆]²⁺의 ν₃(~395 nm)와
  일치한다.
→ **말할 수 없는 것**: 그 요철이 니켈 d-d 밴드라는 것. 위치가 5개 시료에서 전혀 움직이지 않고,
  20~40 nm 간격 요철 6개 중 하나이며, 세기 차이가 시료 간 베이스라인 오차보다 작기 때문.
  ν₁(1176 nm)을 못 쟀으므로 교차확인 수단도 없음.

→ **따라서 이번 DFT의 역할은 "측정 피크 재현"이 아니라 "예측 제시"입니다.**
  *조성비 → 염화물 활동도 → Ni 화학종 → 예상 d-d 밴드*까지를 계산으로 완결하고,
  실측과는 **정성적 방향성**(염 과잉 쪽에서 사면체종 신호가 커야 함) 수준으로만 대조합니다.
  이 편이 방어 가능하고, 다음 배치 실험의 설계 근거로도 그대로 쓰입니다.

---

## 5. DFT 소요시간 — 추정과 실측 절차

### 5-1. 결론 먼저: **시간은 병목이 아닙니다**

서버는 **128스레드**이고 우리 잡은 **serial(1코어)** 입니다(`CLAUDE.md` 실측: 작은 분자는 MPI가 오히려 느림).
→ 5개를 **동시에** 던져도 5코어밖에 안 씁니다. **총 CPU시간이 아니라 가장 느린 한 잡의 벽시계가 소요시간**입니다.

| 화학종 | 기저함수 | 1단계 Opt 추정 | 2단계 NEVPT2 추정 | 합계(벽시계) |
|---|---|---|---|---|
| `NiCl4` | 193 | 10 ~ 25 분 | 15 ~ 45 분 | **0.5 ~ 1.2 h** |
| `NiCl3_aq1` | 199 | 20 ~ 45 분 | 20 ~ 50 분 | **0.7 ~ 1.6 h** |
| `NiCl2_aq4` | 291 | 1.5 ~ 3 h | 1 ~ 3 h | **2.5 ~ 6 h** |
| `NiCl1_aq5` | 297 | 2 ~ 4 h | 1 ~ 3 h | **3 ~ 7 h** |
| `Ni_aq6` | 303 | 2 ~ 4 h | 1 ~ 3 h | **3 ~ 7 h** |
| | | | **5개 동시 실행 벽시계** | **3 ~ 7 h** |

추정 근거:
- **실측 앵커**: ChCl 이온쌍(21원자, BP86/def2-TZVP, 같은 서버, serial) 전체 파이프라인이
  **3시간 58분** — 그것도 300 opt cycle까지 간 병적인 케이스 포함(`HANDOFF_NOTES.md` A절).
  `Ni_aq6`(19원자, 303 BF)는 그와 같은 급입니다.
- **Opt cycle 수**: `NiCl4`는 Td 강직체 → 10~15 cycle. `Ni_aq6`는 물 6개의 회전 자유도 → 40~80 cycle.
  이 차이가 원자수 차이보다 크게 작용합니다.
- **NEVPT2 스케일링**: CAS(8,5)는 CI 공간이 작아 비용은 적분변환이 지배 → 대략 N_bf⁴.
  (303/193)⁴ ≈ **6배**.
- 부정확한 부분: Ni²⁺ + CPCM의 SCF 수렴 난이도, 25 roots CASSCF의 root-following 실패 가능성.
  **이게 진짜 리스크입니다** (아래 5-3).

**마감이 9/29(월)이고 오늘이 9/23(화) → 약 6일. 계산시간으로는 5개 전부 여유입니다.**

### 5-2. 추정을 실측으로 바꾸는 절차 (오늘 바로)

```bash
cd ~/.../dft_ni_complexes
python build_inputs.py            # 수정본 반영
sbatch --array=1-1 run_array.sh   # NiCl4 하나만 (최저비용 교정 잡)
# 끝나면
cat timing.log                    # -> "NiCl4 opt=XXXs dd=YYYs"
```

그 다음 외삽:
```
Ni_aq6 opt  ≈ NiCl4_opt × (303/193)^3 × (60/13)      # 기저함수^3 × opt cycle 비
Ni_aq6 dd   ≈ NiCl4_dd  × (303/193)^4  ≈ NiCl4_dd × 6
```
`NiCl4` 합계가 1시간을 넘으면 → 3개로 축소. 30분 이내면 → 5개 전부 제출.

### 5-3. 진짜 병목은 CPU가 아니라 **수렴 실패 재시도**

실패 시 대응(각 1회당 반나절 소요로 예산 잡을 것):

| 증상 | 대응 |
|---|---|
| SCF 미수렴 (특히 `Ni_aq6` q=+2) | `! SlowConv` 추가, 또는 `%scf maxiter 500 end` |
| Opt가 300 cycle 소진 | XTB2 사전최적화 경유 (ChCl 때 쓴 그 수법, `HANDOFF_NOTES.md` A절) |
| CASSCF가 엉뚱한 오비탈 수렴 | `%casscf ... end` 앞에 `%moinp` 로 1단계 궤도 읽어오기 + `orca_plot`으로 활성공간 육안 확인 |
| 25 roots 메모리 초과 | `nroots`를 10,0 으로 낮춤 (삼중항만 — **가시영역 spin-allowed 밴드는 전부 삼중항이므로 실질 손실 없음**) |

> **마지막 항목이 사실상 안전장치입니다.** 정 급하면 `mult 3` / `nroots 10`만 돌려도
> 우리가 설명해야 할 밴드(³A₂g→³T₁g 등 spin-allowed)는 전부 나옵니다. 단일항 15개는
> 매우 약한 spin-forbidden 밴드용이라 UV-Vis 설명에는 없어도 됩니다.

### 5-4. 권장 일정

| 날짜 | 할 일 |
|---|---|
| 9/23(화) | `--array=1-1`(NiCl4) 제출 → `timing.log` 확보 |
| 9/24(수) | 실측 기반 재추정. 문제없으면 `--array=2-5` 제출(나머지 4종 동시) |
| 9/25(목)~9/26(금) | 결과 수집, AILFT의 10Dq/Racah B 추출, 화학종별 흡수 스펙트럼 작도 |
| 9/26(금) | **UV-Vis 재측정**(§4) — 이게 제일 리드타임이 긺. 민이 형과 일정 먼저 잡을 것 |
| 9/27~28(주말) | 계산 스펙트럼 vs 실측 대조, 조성비→화학종 분포 논리 정리 |
| 9/29(월) | 보고 |

**임계경로는 DFT가 아니라 UV-Vis 재측정입니다.**

---

## 6. 실행 절차 (이대로 따라가면 끝)

```bash
cd ~/.../dft_ni_complexes
python build_inputs.py                  # 인풋 재생성 (수정본 반영)

sbatch --array=1-1 run_array.sh         # 1) NiCl4 하나만 = 교정 잡
cat timing.log                          #    -> "NiCl4 opt=XXXs dd=YYYs"

sbatch --array=2-5 run_array.sh         # 2) 나머지 4종 동시 (또는 --array=2-3 만)

python extract_dd.py                    # 3) 결과 수확
```

`extract_dd.py`가 만드는 것:

| 산출물 | 내용 |
|---|---|
| 화면 출력 | 화학종별 d-d 전이 목록 (nm, cm⁻¹, 진동자세기) + AILFT의 10Dq·Racah B 원문 줄 |
| `dd_bands.csv` | 위 내용의 표 형태 — 영진 형 ML 쪽/보고서 표에 그대로 투입 |
| `dd_spectra.png` | 화학종별 예측 스펙트럼 5개 겹쳐 그린 그림 (395·725 nm 기준선 표시) |

파서는 ORCA 버전에 안 흔들리게 짰습니다: 흡수스펙트럼 표에서 **nm = 10⁷/cm⁻¹ 관계를
만족하는 인접 숫자쌍**을 찾는 방식이라 컬럼 순서가 바뀌어도 동작합니다.
`python extract_dd.py --test`로 ORCA 출력 없이도 파서 검증 가능(통과 확인함).

**발표에 쓸 그림 한 장**: `dd_spectra.png`에서 `Ni_aq6`(팔면체, 약한 밴드)와
`NiCl4`(사면체, 훨씬 강한 밴드)를 나란히 놓으면 "염이 많아지면 색이 강해진다"가
계산만으로 한눈에 보입니다. 이게 조성비→색 논리의 시각적 핵심입니다.

**주의 하나**: CAS(8,5)는 double-shell 누락으로 d-d 에너지를 보통 10~20% 과대평가합니다
(§2-4). 즉 계산 파장이 실측보다 **짧게** 나옵니다. 밴드 *순서*와 *화학종 간 상대 이동*은
신뢰할 수 있지만 절대 파장은 그대로 인용하지 말고, 이 사실을 그림 캡션에 적으세요.

---

## 7. 결과물 위치

| 파일 | 내용 |
|---|---|
| `build_inputs.py` | 인풋 생성 (이번에 수정) |
| `species.txt` | **우선순위 순** 화학종 목록 |
| `run_array.sh` | Slurm 잡배열 + `timing.log` 기록 |
| `*_opt.inp` / `*_dd.inp` | 1/2단계 ORCA 인풋 |
| `extract_dd.py` | **결과 후처리** — d-d 밴드표 + 예측 스펙트럼 생성 (`--test`로 자체검증) |
| `dd_bands.csv` / `dd_spectra.png` | 위 스크립트 산출물 (계산 후 생성) |
| `../uvvis_sanity.py` | UV-Vis 건전성 점검 (§4 한계 수치를 생산한 코드) |
| `../uvvis_sanity.png` | raw vs 베이스라인 보정 비교 그림 |
