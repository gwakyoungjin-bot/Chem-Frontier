# ML 파이프라인 인계 자료 (영진 형 전달용 재료)

> 이 문서는 **완성된 전달 문장이 아니라 사실관계·수치·근거 모음**입니다.
> 실제 메시지/보고서 문장은 직접 작성하는 것을 전제로 합니다.
> 최종 갱신: 2026-09-08

---

## A. 데이터 변경 이력

| 버전 | 시점 | 내용 | 쌍 수 | 행 수 | 컬럼 수 |
|---|---|---|---|---|---|
| v1 | 2026-09-01 | 중성분자 22종 전체 쌍 조합 | 231 | 21,945 | 62 |
| v2 | 2026-09-08 | `is_ionic_hba` 플래그 컬럼 추가 (전부 False) | 231 | 21,945 | 63 |
| **v3** | **2026-09-09** | **콜린클로라이드(ChCl) × 22종 추가 — 현재 전달본 (`features_all.csv`)** | **253** | **24,035** | **23물질** |

> **v3 계산 경위**: ChCl은 이온쌍이라 가스상 기하최적화 수렴에 애를 먹어 두 번 실패했습니다.
> ① XTB2 사전최적화 도입(N–Cl 6.07 → 3.65 Å 교정)으로 수렴 지표가 100배 개선(RMS gradient 0.0108 → 0.00032),
> ② 사이클 상한을 66 → 300으로 상향(`%geom MaxIter 300`)하여 해결. 최종 ORCA 계산 3시간 58분 소요.
> 최종 구조에서 이온쌍이 정상 유지됨(N–Cl 4.65 Å, O–Cl 3.06 Å 수소결합 유지).

**대상 물질 22종의 출처**
- 랩실 재고 기반 11종: ethylene glycol, citric acid, diethylene glycol, glycerol, decanoic acid, cinnamic acid, urea, menthol, boric acid, imidazole, cyclohexanol
- 문헌 기반 11종 (6편에서 확인): ascorbic acid, lactic acid, malonic acid, oxalic acid, malic acid, tartaric acid, glycine, D-glucose, D-fructose, sorbitol, trifluoroacetamide
  - Amusa 2025 / Zurob 2020 / Lu 2022 / Thompson 2022 / Pal & Jadeja 2019 / Lemaoui 2020

**변수 그리드 (쌍당 95개 조합)**
- 몰분율 `x_hba`: 0.05 ~ 0.95, 0.05 간격 (19개)
- 온도 `T`: 298.15 / 318.15 / 338.15 / 358.15 / 373.15 K (= 25/45/65/85/100 °C, 5개)

---

## B. 컬럼 스키마

| 그룹 | 컬럼 | 내용 |
|---|---|---|
| 식별 | `hba_name`, `hbd_name` | 두 성분의 이름 (※ 화학적 역할 아님 — D-1 참조) |
| 조건 | `x_hba`, `T` | 몰분율, 온도(K) |
| 구분 | `is_ionic_hba` | 이온성 HBA(ChCl 등) 포함 여부 |
| 활동도계수 | `lng_*_tot` | 전체 ln γ (핵심 값) |
| | `lng_*_enth` | 잔여항 — 분자 간 실제 상호작용 기여 |
| | `lng_*_comb` | 조합항 — 분자 크기·모양 차이 기여 |
| 초과물성 | `gE_over_RT` | G^E/RT = Σ xᵢ·ln γᵢ. 0에서 벗어날수록 비이상적 |
| 에너지 분해 | `*_pm_E_hb` | 수소결합 기여 (J/mol) |
| | `*_pm_E_mf` | misfit(정전기 불일치) 기여 |
| | `*_pm_A_int` | 총 상호작용 자유에너지 |
| 공융점 | `T_eutectic_pred`, `x_eutectic_pred`, `delta_x_from_eutectic` | **현재 전부 비어있음** (D-3 참조) |
| 분자 descriptor | `*_sigma_profile_bin0~10` | σ-profile 11구간 적분값 — **ML 핵심 입력** |
| | `*_sigma_moment0~6` | σ-profile 압축본 7개 |
| | `*_area`, `*_volume`, `*_dipole_moment_mag` | 표면적(Å²), 부피(Å³), 쌍극자모멘트 |

`hba_*` / `hbd_*` 접두어가 붙은 descriptor는 **분자 고유값**이라 같은 물질이면 몰비·온도와 무관하게 동일합니다.

---

## C. 설계 결정과 근거

**C-1. σ-profile을 ML 입력으로 쓰는 이유**
- 벤치마크 논문들이 쓴 입력이 ln γ가 아니라 σ-profile임: Lemaoui 2020(R²=0.985, 전도도), Mohan 2024(R²=0.99, 점도), Wang 2021
- σ-profile 구간을 −0.025~+0.025 e/Å² 폭 0.005로 자르는 것도 이들과 동일한 스킴

**C-2. 몰비를 촘촘하게 스캔할 수 있는 이유**
- Lemaoui 2020 Eq.(2): DES descriptor는 성분 descriptor의 몰분율 가중합
- → 화합물당 COSMO 계산 **1회**면 몰비는 후처리만으로 무한히 세분화 가능. 계산비용이 몰비 축에서는 거의 0
- 조합폭발은 "어떤 화학종을 쓸지" 축에서만 발생 (Velez & Acevedo 2022)

**C-3. NaCl-EG를 시뮬레이션에서 제외한 이유**
- 이온 용액은 장거리 정전기(Debye-Hückel/Pitzer) 항이 필요한데, 오픈소스 `opencosmorspy`에는 없음. 이를 갖춘 확장판(COSMO-RS-ES, COSMO-RS-PDHS)은 전부 상용 COSMOtherm 전용
- 실제로 Na⁺/Cl⁻를 분리 이온으로 넣어 테스트한 결과가 비물리적이었음:

  | 값 | 결과 | 중성쌍 정상 범위 |
  |---|---|---|
  | `lng_anion_tot` | −80.06 (γ≈10⁻³⁵) | 대략 −1 ~ +1 |
  | `pm_E_hb` | −60,650 J/mol | 수소결합 하나에 60 kJ/mol은 비현실적 |
  | `gE_over_RT` | −12.03 | 대략 −0.1 |

- NaCl-EG는 습식실험 트랙에서만 진행 (NMC811 직접침출 선행연구로 이미 정당화되어 COSMO-RS 근거가 불필요)

**C-4. ChCl은 이온쌍 1성분으로 처리하는 이유**
- Lemaoui 2020 Eq.(2) 원문: "NC is the total number of components in the DES mixture, and j represents a particular component **including the HBA, the HBD**, or any other components"
- 그들의 HBA 4종(BTPPCl, MTPPBr, ChCl, DEACl)이 전부 염인데도 각각 **1개 성분**으로 취급
- → 우리도 동일하게 처리 → 출력 스키마가 기존과 완전히 동일 → 한 파일로 합칠 수 있음

**C-5. 삼성분계(EG+CA+Ethanol) 미실시** — 팀 논의로 제외 결정

---

## D. ML 학습 시 주의사항 (중요)

**D-1. ⚠️ `hba_name`/`hbd_name`은 화학적 역할이 아니라 위치 슬롯입니다 — 전처리 필수**

조합 생성 시 리스트 순서대로 넣어서 슬롯 배정이 치우쳐 있습니다:

| 물질 | slot1(`hba_*`) | slot2(`hbd_*`) |
|---|---|---|
| ethylene_glycol | 21 | **0** |
| citric_acid | 20 | 1 |
| glycine | 4 | 17 |
| sorbitol | 1 | 20 |
| trifluoroacetamide | **0** | 21 |

문제 두 가지:
- **커버리지 공백**: EG의 descriptor가 slot2에 한 번도 등장하지 않음 → 그런 입력이 오면 외삽
- **대칭성 미학습**: (A,B,x=0.3)과 (B,A,x=0.7)은 같은 물리계인데 다른 벡터로 보임

해결책 (택 1, 전처리 단계). **아래 두 코드는 실제 `features_broad.csv`로 실행 검증했습니다.**

---

**방법 1 — Lemaoui 2020 방식 가중합 (권장)**

`S_DES = Σ xⱼ·Sⱼ` (성분 descriptor의 몰분율 가중합). 가중합은 순서에 무관하므로 슬롯 문제가 원천 소멸하고, 벤치마크 논문과 형식이 일치해 비교가능성도 확보됩니다.

```python
import pandas as pd

df = pd.read_csv('features_broad.csv')
x1 = df['x_hba'].to_numpy(float)
x2 = 1.0 - x1

# hba_/hbd_ 짝이 맞는 '수치형' 컬럼만 추출 (hba_name 같은 문자열 제외)
suf = [c[4:] for c in df.columns
       if c.startswith('hba_') and ('hbd_' + c[4:]) in df.columns
       and pd.api.types.is_numeric_dtype(df['hba_' + c[4:]])]

des = pd.DataFrame({
    'x_hba': x1,
    'T': df['T'],
    'is_ionic_hba': df['is_ionic_hba'],
    'pair': df['hba_name'] + '|' + df['hbd_name'],   # 분할용 그룹키 (D-2)
})
for s in suf:
    des['des_' + s] = x1 * df['hba_' + s].to_numpy(float) + x2 * df['hbd_' + s].to_numpy(float)
```

- 검증 결과: 24개 descriptor가 합쳐져 **21,945행 × 28열**, 결측 0
- 대상 컬럼: `sigma_profile_bin0~10`, `sigma_moment0~6`, `area`, `volume`, `dipole_moment_mag`, `pm_E_hb`, `pm_E_mf`, `pm_A_int`
- `lng_*`(활동도계수)는 짝 구조가 아니라 성분별 출력값이므로 이 변환 대상이 아닙니다 — 목적변수로 쓰거나 별도 처리하세요

---

**방법 2 — 미러 증강**

각 행에 (성분 스왑 + x → 1−x) 행을 추가해 대칭성을 명시적으로 학습시킵니다. descriptor를 분리된 채로 두고 싶을 때 사용합니다.

```python
import pandas as pd

df = pd.read_csv('features_broad.csv')

# hba_* <-> hbd_* 컬럼명 스왑 매핑
ren = {}
for c in [c for c in df.columns if c.startswith('hba_')]:
    s = c[4:]
    if ('hbd_' + s) in df.columns:
        ren[c] = 'hbd_' + s
        ren['hbd_' + s] = c

mirrored = df.rename(columns=ren).copy()
mirrored['x_hba'] = 1.0 - df['x_hba'].to_numpy(float)

aug = pd.concat([df, mirrored[df.columns]], ignore_index=True)
```

- 검증 결과: 21,945 → **43,890행**, 슬롯 편향 완전 해소 (EG 21/21, trifluoroacetamide 21/21)
- 주의: 증강 후에도 D-2의 그룹 분할은 그대로 필요하며, **미러 쌍이 train/test로 찢어지지 않도록** 그룹키를 정렬된 쌍 이름으로 만드세요:
  ```python
  aug['pair'] = ['|'.join(sorted([a, b])) for a, b in zip(aug.hba_name, aug.hbd_name)]
  ```

---

**방법 3 — 정준 순서 고정**: 분자량 등 결정적 기준으로 항상 정렬. 임의성은 없애지만 대칭성을 가르치진 못해 위 두 방법보다 약합니다.

> 원본 CSV가 틀린 게 아니라, 성분별로 분리 보관해서 **정보량이 더 많은 형태**입니다. 위 방법 모두 이 원본에서 유도 가능합니다. 다만 **무처리 학습은 안 됩니다.**

**D-2. ⚠️ train/test split은 반드시 "쌍 단위(GroupKFold)"로**

한 쌍당 95개 행은 독립 샘플이 아니라 같은 화학계의 격자점입니다(몰분율 19 × 온도 5). 행 단위로 무작위 분할하면 **같은 쌍이 train과 test에 동시에** 들어가서, 모델이 "그 쌍을 이미 봤기 때문에" 맞히게 됩니다 → **데이터 누수**, 성능 과대평가.

일반화 성능을 판단하는 실질 표본 수는 21,945가 아니라 **231(→253)** 입니다.

**방법 A — 쌍 단위 분할 (최소 요건)**

```python
from sklearn.model_selection import GroupKFold

groups = df['hba_name'] + '|' + df['hbd_name']

gkf = GroupKFold(n_splits=5)
for train_idx, test_idx in gkf.split(df, y, groups=groups):
    X_tr, X_te = X.iloc[train_idx], X.iloc[test_idx]
    y_tr, y_te = y.iloc[train_idx], y.iloc[test_idx]
    ...
```

검증 결과 (5-fold): 각 fold가 train 약 185쌍 / test 약 46쌍으로 나뉘고, **train-test 중복 쌍 0개** 확인.

| fold | train | test | 중복 쌍 |
|---|---|---|---|
| 0 | 17,480행 / 184쌍 | 4,465행 / 47쌍 | 0 |
| 1~4 | 17,575행 / 185쌍 | 4,370행 / 46쌍 | 0 |

**방법 B — 분자 단위 분할 (더 엄격, QSPR 일반화 주장을 하려면 권장)**

방법 A는 쌍은 다르지만 **같은 분자**가 train과 test에 (다른 짝으로) 등장할 수 있습니다. "처음 보는 분자에도 통하는가"를 검증하려면 분자 단위로 빼야 합니다.

```python
test_mols = {'sorbitol', 'glycine', 'menthol'}   # 예시
is_test = df.hba_name.isin(test_mols) | df.hbd_name.isin(test_mols)

X_tr, y_tr = X[~is_test], y[~is_test]
X_te, y_te = X[is_test],  y[is_test]
```

검증 결과: train 16,245행 / test 5,700행, **test 분자가 train에 전혀 등장하지 않음** 확인.

> 논문·보고서에 "일반화된다"고 쓰려면 방법 B 결과를 함께 제시하는 게 설득력이 훨씬 큽니다. 방법 A만으로는 "본 적 있는 분자들의 새로운 조합"까지만 검증됩니다.

**D-3. 공융점 컬럼 3개는 현재 비어있음**
- `T_eutectic_pred` / `x_eutectic_pred` / `delta_x_from_eutectic`
- 계산에 성분별 융해엔탈피(ΔH_fus)·융점(T_fus) 문헌값이 필요한데 아직 수집 안 됨
- solver(`eutectic.py`)는 작동 검증 완료 (이상용액 자체검증 통과)

**D-4. 전부 계산값입니다** — 실측이 아닙니다. 실측 데이터는 민이 형 실험의 침출효율뿐이고, 이 CSV는 그 정답(y)과 연결할 입력(X)입니다.

**D-5. ChCl 행의 용도** — 학습 데이터 증량보다 **문헌 대조 검증 앵커**로서의 가치가 큽니다
- ChCl:EG(에타린), ChCl:urea(릴라인), ChCl:glycerol(글리셀린)은 문헌 실측값이 가장 풍부한 3대 DES
- 우리 데이터셋에 citric/malonic/oxalic acid가 이미 있어서, ChCl 추가 시 Peeters 2020·Lu 2022·Thompson 2022 논문 시스템이 즉시 재현됨
- `is_ionic_hba == True` 로 필터링해서 검증셋으로 분리 가능

---

## E. 한계 (보고서에 명시할 것)

1. **이온쌍 근사** — 염 HBA를 중성 이온쌍 1성분으로 취급. 실제 ChCl DES는 전기가 통하므로(=이온이 해리·이동) 물리적으로는 부정확한 근사. 문헌 표준을 따른 것이며 QSPR 목적에는 검증된 방식(Lemaoui 2020 R²=0.985)
2. **콘포머 1개 제약** — 라이브러리가 분자당 콘포머 1개만 지원(`molecules.py:28`). 이온쌍은 Cl⁻의 병진 자유도까지 있어 이 제약이 특히 크게 작용
3. **MMFF가 이온쌍 구조를 구별하지 못함** — 콜린클로라이드 콘포머 5개의 Cl⁻ 위치가 2.02~7.94 Å로 제각각인데 MMFF 에너지가 전부 39.06 kcal/mol로 동일. 이온 간 정전기를 계산하지 않는다는 뜻 → 시작구조를 수소결합 모티프(O···Cl=3.10 Å)로 직접 구성해 사용
4. **시뮬레이션 대상 ≠ 실험 대상** — 실험 후보 2번(NaCl-EG)은 시뮬레이션 불가라 제외. 따라서 "이 물질의 성능"이 아니라 "DES 조성계 전반의 구조-물성 관계"를 학습하는 구조
5. **231쌍 중 상당수는 실제 공융물 형성 여부 미검증** — 넓은 화학공간 학습용으로 생성한 조합이며, 전부가 실재하는 DES는 아님

---

## F. 데이터 완결성 검증 결과 (2026-09-08 실시)

**계산 단계 완주 확인** — 22종 전부 최종 단계까지 완료됨
- 각 분자의 `.orcacosmo` 헤더 method가 전부 `DFT_CPCM_BP86_def2-TZVP+def2-TZVPD_SP` = 파이프라인 **마지막 단계** 산출물
- CPCM 보정 섹션 / 표면 세그먼트 / 좌표 / 쌍극자모멘트 전부 포함 확인 (22/22)
- 중간 단계 작업폴더가 전부 정리된 상태 = 정상 완료 신호 (실패 시에는 단계 폴더가 남음)

**CSV 무결성 확인**

| 검사 항목 | 결과 |
|---|---|
| 형상 | 21,945행 × 63열 (기대치 일치) |
| 고유 쌍 수 | 231 (기대치 일치) |
| 쌍당 행 수 | 전부 정확히 95행 (누락된 격자점 없음) |
| 고유 물질 수 | 22 |
| 몰분율 / 온도 격자 | 19개 / 5개 (설계대로) |
| 핵심 계산값 13개 컬럼 | 결측 0, 무한대 0 |
| descriptor 42개 컬럼 | 결측 0 |
| `ln γ` 범위 | −6.47 ~ +5.91 (유기 혼합물로서 타당) |
| `gE_over_RT` 범위 | −2.03 ~ +1.39 |
| `pm_E_hb` 범위 | −39.8 ~ +4.7 kJ/mol (수소결합 에너지로서 타당) |

> 참고: 제외한 NaCl 이온계 테스트에서는 ln γ = −80(γ≈10⁻³⁵), pm_E_hb = −60,650 J/mol 같은
> 비물리적 값이 나왔습니다. 위 범위는 그와 대조적으로 정상 범위입니다.

**ChCl(v3) 추가분 검증 — 2026-09-09**

구조: 2,090행 × 63열, 22쌍 전부 95행, `is_ionic_hba=True`, 컬럼 순서까지 v2와 동일, 결측·무한대·중복 0.

**물리적 타당성 — 세 가지가 독립적으로 확인됨**

① **열역학 극한 거동이 정확함** (ChCl:urea, 298 K 예시)
- lnγ(ChCl): 무한희석 x=0.05에서 −6.48 → 순수극한 x=0.95에서 −0.002로 매끄럽게 0 수렴
- lnγ(urea): 정확히 거울상. `gE_over_RT`는 단일 최소를 갖는 정상적인 초과자유에너지 곡선

② **알려진 DES 조성을 계산이 스스로 재현** (문헌 3종 모두 1:2 몰비 = x_ChCl 0.333)

| 시스템 | 계산된 gE 최소점 | 문헌 조성 |
|---|---|---|
| ChCl:urea (릴라인) | 0.35 | 0.333 |
| ChCl:glycerol (글리셀린) | 0.40 | 0.333 |
| ChCl:ethylene glycol (에타린) | 0.40 | 0.333 |

격자 간격이 0.05이므로 셋 다 **한 칸 이내** 일치. 조성 정보를 입력한 적 없음.

③ **문헌상 LIB 침출에 실제 쓰인 조합이 상위 랭크** (22개 중 상호작용 강도 순)
- oxalic acid (Thompson 2022) **1위** (gE/RT −5.85)
- citric acid (Peeters 2020) **2위** (−4.69)
- malonic acid (Lu 2022) 6위 (−3.62)
- 반대로 menthol(+0.07), cyclohexanol(+0.02) 등 부피 큰 소수성 알코올은 최하위 = 염과 DES를 못 만든다는 화학적 상식과 일치

**값 범위가 중성 22종보다 넓은 이유 (정상)**

| 지표 | ChCl | 중성 22종 |
|---|---|---|
| `lng_hba_tot` | −29.65 ~ 3.87 | −6.47 ~ 5.91 |
| `gE_over_RT` | −5.85 ~ 0.99 | −2.03 ~ 1.39 |
| `hba_pm_E_hb` | −149,750 ~ 1,341 J/mol | −39,766 ~ 4,731 J/mol |

극값은 전부 **oxalic acid의 x=0.05, 즉 무한희석 지점**에서 발생합니다. 부분몰량은 무한희석에서 발산하는 것이 정상이며, 옥살산은 이 조합 중 가장 강한 산(pKa₁ 1.25)이라 염화물과의 상호작용이 가장 강한 것이 화학적으로 타당합니다.
제외한 NaCl 분리이온 방식은 이와 달리 **일반 조성에서도** 값이 붕괴했습니다(lnγ −80 등) — 성격이 다릅니다.

**미완성 항목 (의도된 것)**
- 공융점 3컬럼(`T_eutectic_pred`, `x_eutectic_pred`, `delta_x_from_eutectic`)은 전부 비어 있음 — 성분별 ΔH_fus/T_fus 문헌값 미수집 (D-3 참조)

---

## G. 파일 위치

| 파일 | 내용 |
|---|---|
| `features_broad.csv` | 최종 데이터셋 (21,945행 × 63컬럼) |
| `extract_features.py` | 데이터 생성 스크립트 |
| `eutectic.py` | 공융점 solver (SLE 곡선 교점) |
| `pairs_broad_server.csv` | 231쌍 조합 목록 |
| `chcl_prep/choline_chloride.inp` | ChCl ORCA 계산 입력 |
| `chcl_prep/choline_chloride_start.xyz` | ChCl 이온쌍 시작구조 (검증 완료) |
| `chcl_prep/pairs_chcl_server.csv` | ChCl 22쌍 목록 (`is_ionic_hba=True`) |
| `orcacosmo/` | 분자 24종의 COSMO 계산 결과 |
