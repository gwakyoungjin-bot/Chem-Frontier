# Chem Frontier

Choline chloride-citric acid(ChCl-CA) 기반 수계 침출 조성(DES 후보 formulation)의 몰비 탐색을 위해
COSMO-RS 계산, UV-Vis 실험, Gaussian process 회귀를 순차적으로 연결한 프로젝트입니다.

이 README와 `docs/PROJECT_SUMMARY.md`는 그중 **COSMO-UV-ML 작업 흐름에서 함께 수행한 내용만** 정리합니다. 저장소의 `simulation/` 아래에 별도로 올라온 DFT·분광 작업 결과는 이 문서에서 재해석하거나 검증하지 않습니다.

## 정리 문서

- [COSMO-UV-ML 작업 요약](docs/PROJECT_SUMMARY.md): 연구 목적, 데이터, COSMO 기반 초기설계, 1-2차 실험, ML 검증, 한계와 별도 DFT 작업의 연결 범위

## 현재 단계

1. ChCl-CA의 19개 조성격자에 대해 COSMO-RS descriptor 공간을 구성했습니다.
2. 실험 정답이 없는 상태에서 조성거리와 COSMO 거리를 이용해 1차 5몰비를 선정했습니다.
3. 1차 UV-Vis 결과 5개로 고정 RBF-kernel Gaussian process 모델을 학습했습니다.
4. 2차 결과를 보기 전에 모델과 19격자 예측을 고정하고, 신규 5몰비의 실측값과 먼저 비교했습니다.
5. 이후 1-2차의 10개 고유 몰비로 모델을 갱신하고 교차검증, y-scrambling, 적용영역 및 독립 수치 재계산을 수행했습니다.

> 현재 모델의 목표값은 300-400 nm 최대흡광도인 상대적 UV proxy입니다. 검량된 Ni 농도 또는 침출률(%) 모델이 아니며, COSMO descriptor의 추가 예측 이점도 현재 데이터에서는 입증되지 않았습니다.

## 연구 흐름

```text
COSMO-RS descriptor space
        -> label-free initial composition design
        -> round-1 UV-Vis experiment
        -> fixed Gaussian process surrogate
        -> frozen prospective prediction check
        -> round-2 experiment
        -> 10-composition model update
```

## 핵심 원칙

- 비이상성 지표 자체를 침출 성능으로 간주하지 않습니다.
- 파장점을 독립 실험으로 부풀리지 않고, 몰비 하나를 독립 조건 하나로 셉니다.
- 새 실험값은 기존 예측을 먼저 평가한 뒤 학습에 편입합니다.
- 코드 재현성과 화학적 예측력을 구분합니다.
- 실패한 검증지표와 모델의 한계를 숨기지 않습니다.
