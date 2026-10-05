# 원문에서 뽑은 앵커 수치 (2026-09-27)

`README.md` 는 "왜 이 논문들을 모았나", 이 파일은 **"원문에서 실제로 뽑은 숫자"** 다.
전부 본문에서 직접 확인한 값이며, 초록만 보고 적은 것은 없다.

---

## 1. Busato 2022 — Ni 앵커 (`10.1021/acs.inorgchem.2c00864`)

**계**: NiCl₂·6H₂O : urea = 1 : 3.5 (MDES), 물/MDES 몰비 **W** 를 0 → 26 으로 스캔
**측정**: 실온(room temperature). MD + CASSCF/NEVPT2 + UV-Vis/NIR + SAXS/WAXS + **Ni K-edge XAS**

| 조건 | 값 |
|---|---|
| W = 0 (순수 MDES) Ni–Cl 배위수 | **평균 1.3** (물 4.5, urea 0.3) |
| Cl 배위 순간분포 | 1개 **43.8%** / 2개 30.9% / 0개 20.6% |
| **사면체 분율** | **0** — 원문: *"direct proof of the fully octahedral coordination of the Ni²⁺ ion"* |
| W 증가 시 | 물이 배위권을 채워 **hexa-aquo** [Ni(H₂O)₆]²⁺ 로 수렴 |
| 물 몰분율 (W=0, 환산) | ≈ 0.57 (NiCl₂ 1 + H₂O 6 + urea 3.5 기준) |

### 우리 모델에 대한 함의

**매우 강한 제약이다.** 염화물이 극도로 풍부한(Ni 당 Cl 2개) 금속-DES 인데도
실온에서 Ni 가 **완전히 팔면체**다. 우리 모델은 25 °C 에서 f_Td ≈ 0.0003 을 주므로
**이 제약을 통과한다 ✅** (Ni 쪽은 문제없음).

### 방법론 인용 근거 (중요)

이들도 **CASSCF(8,5) 및 (14,8) / SC-NEVPT2** 를 `[Ni(H₂O)₆]²⁺, [NiCl(H₂O)₅]⁺,
[NiCl₂(H₂O)₄] cis·trans, [NiCl₃(H₂O)₃]⁻, [NiCl₄(H₂O)₂]²⁻` 에 적용했다.
**우리와 같은 이론수준·같은 화학종 계열.** "왜 CASSCF/NEVPT2 인가" 질문에
Inorg. Chem. 게재 선례로 답할 수 있다.
(차이: 이들은 전부 6배위로 두었고, 우리는 Td 종을 4배위로 둠)

### 검출 논증에 그대로 쓸 문장

> *"the absorption intensities connected to the octahedral coordination are orders of
> magnitude lower than those found in the tetrahedral case, due to the presence of an
> inversion center... The formation of even little tetrahedral coordination in solution
> would be therefore easily detectable due to the appearance of the well-known double
> band between 600 and 800 nm"*

→ 우리 "팔면체는 Laporte 금지라 안 보인다 / 사면체면 보였어야 한다" 논증의 직접 근거.

---

## 2. Iname 2026 — Ni·Co 동시 (`10.1016/j.molliq.2026.129750`)

**계**: ChCl 기반 DES 5종. **HBD 를 바꿔가며** Ni²⁺·Co²⁺ 화학종을 UV-Vis 로 판정.
**측정**: 280–780 nm, 1 cm 셀, PerkinElmer Lambda 365.
⚠ **"A blank was then run to eliminate interfering signals for each solvent"**
 — **용매마다 블랭크를 따로 잡았다.** 우리가 민이형께 요청한 바로 그 방식.

### Table 6 (원문) — Ni²⁺ 화학종

| DES | λmax / nm | 강도 | 색 | 기하 |
|---|---|---|---|---|
| ChCl–**lactic acid** (1:2) | 425.3 강 / 708.3 약 / 764.2 약 | 약 | — | **Oh** |
| ChCl–glycerol (1:2) | 405.8 강 / 680 약 / 746.7 약 | 약 | — | Oh |
| ChCl–EG (1:2) | 420.9 강 / 691.4 약 / 760.5 약 | 약 | 녹색 | Oh |
| ChCl–urea (1:2) | 417.8 강 / 680 약 / 753.2 약 | 약 | — | Oh |
| **ChCl–malonic acid (1:1)** | **685.2 강 / 708.1** | **강** | **파랑** | **Td [NiCl₄]²⁻** |

### Co²⁺

**5종 DES 전부에서 [CoCl₄]²⁻ (사면체).** 원문: *"Co²⁺ preferentially forms
tetrachlorocobaltate(II) complexes, [CoCl₄]²⁻, in all DESs systems studied"*

→ **우리 Co 계산 방향(Co 는 쉽게 사면체)이 지지된다 ✅**
(단 Mannucci 는 물을 많이 넣으면 Oh 로 돌아간다고 정량했으므로, 둘은 모순이 아니라
 "건조 DES 에서는 Td, 물이 많으면 Oh" 로 일관된다.)

### 농도 — 간접적이지만 결정적인 정보

| 계 | 쓴 Ni 농도 |
|---|---|
| Oh 나오는 DES 4종 | **10 ~ 50 mmol/kg** |
| **Td 나오는 ChCl–malonic** | **1 ~ 5 mmol/kg** (10배 낮춤) |

사면체가 훨씬 세게 흡수하므로 포화를 피하려고 10배 희석한 것. **ε 차이의 실측 증거.**

---

## 3. ⚠ 이전 주장 철회

**철회**: "ChCl : 카르복실산 계의 금속 화학종 연구를 못 찾았다" (2026-09-27 오전)

**틀렸다.** Iname 2026 이 **ChCl–락트산**과 **ChCl–말론산** 두 카르복실산 DES 의
Ni·Co 화학종을 정면으로 다룬다. OpenAlex 에 이 논문 초록이 없어서 제목만으로는
판별이 안 됐던 것이고, 원문을 보고서야 확인됐다.

→ 갭 주장을 **"선행연구가 없다"** 가 아니라
   **"시트르산(삼카르복실산) 계는 아직 없고, 이카르복실산에서는 니켈이 사면체가 된다는
   보고가 있어 우리 계를 그냥 유추할 수 없다"** 로 바꿔야 한다. (더 정확하고 여전히 유효)

---

## 4. 새로 생긴 위험과 새로 생긴 결과

### ⚠ 위험: 폴리카르복실산이면 Ni 가 사면체일 수 있다

말론산 = 이카르복실산 → Ni **사면체**. 시트르산 = **삼**카르복실산.
우리 지도는 "60 °C 에서 Ni 는 팔면체"라고 하는데, HBD 만으로 뒤집힐 가능성이 생겼다.

**다만 우리 편인 요소도 있다:**
- 우리 계는 **물 몰분율 0.80** — Iname 의 DES 는 사실상 무수. 물은 Oh 를 강하게 밀어준다
- ChCl–malonic 만 **1:1** 이고 나머지는 1:2 → x(ChCl) 0.50 vs 0.333 이라 **조성 교란**이 있다
- 저자들의 기전 설명("stronger-field ligand 라서 Td")은 배위장 이론과 방향이 반대라 신뢰도가 낮다

### ✅ 결과: 우리 음성 UV-Vis 가 이제 의미를 갖는다

Iname 은 **1~5 mmol/kg** 의 사면체 Ni 에서 685/708 nm 에 **강한** 밴드를 본다.
우리 시료는 총 [Ni] ≈ **34 mM** (100% 침출 시)로 **10배 이상 진하다.**

→ **우리 계에서 Ni 가 사면체였다면 밴드가 안 보일 수 없다. 안 보였다 = 사면체가 아니다.**

이건 계산이 아니라 **실측에 근거한 진술**이다. 우리 지도의 "60 °C 에서 Ni 는 팔면체"
예측과 일치한다. 위 ⚠ 위험에 대한 우리 쪽 반론이기도 하다.

### ✅ 블랭크 실험의 근거가 훨씬 강해졌다

Iname 은 **10~50 mmol/kg 의 팔면체 Ni 에서 680~764 nm 약한 밴드를 본다.**
우리 농도(34 mM)는 그 범위 안이다. 그런데 우리는 아무것도 못 봤다.
차이는 **블랭크뿐이다** — 그들은 용매별 블랭크, 우리는 물.

→ 블랭크를 고치면 **팔면체 Ni 밴드가 실제로 보일 가능성이 있다.**
   "그래프가 예뻐진다"가 아니라 **"없던 데이터가 생긴다"** 는 얘기다.
   민이형께 요청할 때 이 근거를 쓸 것.

---

## 5. 남은 것

- Mannucci 2026 의 Co 분율(물 0.80 에서 Td 40%) 을 Co 앵커로 코드에 반영
- Hartley 앵커 활동도를 **그 논문 자신의 계**(ChCl:EG, 물 0.045)에서 계산해 적용
  (계산해둔 값: a_ChCl = 1.83e−1 @ 95 °C)
- Busato 의 W=0 조건(f_Td=0)을 **Ni 쪽 소거 검증**으로 추가 — 이미 통과하지만 명시 가치 있음
