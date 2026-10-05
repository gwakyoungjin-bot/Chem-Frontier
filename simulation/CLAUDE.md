# 케미프론티어 프로젝트 — Claude Code 참고 문서

이 문서는 새 대화 스레드를 시작해도 지금까지의 진행 상황을 이어서 이해할 수 있도록 정리한 것입니다. 작업 전에 반드시 이 문서를 먼저 읽고 시작하세요.

## 프로젝트 개요

2026 Chem Frontier 공모전 (팀 "전지전능", 4인: 동현·영진·민이·민규). 주제: DES(공융용매)/이종 유기산 혼합 침출계의 **조성비(몰비)를 관행적 고정값이 아니라 연속변수로 다뤄 성능(금속 침출효율 등)을 직접 최적화**하는 방법론 — ML 예측·계산화학 시뮬레이션·실제 실험 3자 교차검증. 본선 발표 10/06 BEXCO.

핵심 참고 파일:
- `케미프론티어_워크플로우노트.pdf`, `scratch_workflow_note.html` — 8/5~8/6 워크플로우/역할분담/일정
- `DES_조성비_최적화_선행연구_종합보고서.docx` — 문헌조사 종합본 (37편 검증, 핵심가설: 관행적 몰비=공융점일 뿐 성능최적점 아님)
- `갭분석_지형도_v2.svg` — 문헌 갭 분석 시각화
- `화학물질리스트.xlsx` — (민이 형 쪽 실험 가능 물질 리스트로 추정, 아직 본문 미검토)
- `literature/` — 원문 확보한 PDF들 + 서버 세팅용 오프라인 설치 번들
- 상세 설치 계획 원본: `C:\Users\Donghyeon\.claude\plans\warm-hugging-whale.md`

**2026-09-23 추가 — 새 스레드는 이 두 개를 먼저 읽을 것:**
- `조성비_논리흐름.md` — **발표/보고의 뼈대.** 조성비→염화물 활동도→Ni 화학종→색→성능
  논리 흐름, 단계별 근거 수치, 재현 명령, 인용 문헌, 예상 질문 대응, 남은 작업.
- `dft_ni_complexes/DFT_PLAN.md` — Ni 착물 DFT 설계: 문헌조사 23편 검증표, 화학종 5종 중
  무엇을 돌릴지의 우선순위(`NiCl4`+`Ni_aq6` 우선), 소요시간 추정, 실행 절차, 결과 후처리.

관련 스크립트: `bridge_composition.py`(로컬, 조성→활동도→화학종),
`uvvis_sanity.py`(로컬, UV-Vis 건전성), `dft_ni_complexes/extract_dd.py`(서버, DFT 결과 수확).

## 화학종 앵커 문헌 (2026-09-27 확보)

- `literature/speciation_anchors/` — 화학종 지도를 실측에 고정하기 위한 논문 7편 + **README.md**.
  README 에 논문별 용도, 우리 조건과의 불일치, 뽑아야 할 수치가 정리돼 있음. 앵커 작업 전 먼저 읽을 것.
- **HBD 효과는 가정이 아니라 확인된 사실**: 같은 조성·같은 물함량에서 시트르산이 EG 보다
  a_ChCl 을 **48.7배** 억누른다(물 몰분율 0.045 기준). 물 효과(18~168배)와 같은 자릿수라
  "HBD 는 달라도 된다" 가정은 **기각**. 재현: `python dft_ni_complexes/hbd_effect.py`
- **앵커는 조성이 아니라 활동도에 건다**: 앵커 논문의 계(그 HBD, 그 물함량)에서 a_ChCl 을
  계산해 쓰면 HBD·물 불일치가 동시에 처리된다. Hartley 조건의 값 = 1.83e-1 (95 °C).
- **검색 공백**: ChCl:카르복실산 계 금속 화학종 연구를 못 찾음(독립 검색 3회). 갭 주장 근거로 쓸 것.

## 포스터 주장 범위 (2026-09-27 확정 — 이 선을 넘지 말 것)

> **우리는 화학종 설계 지도를 제시하고, 성능 검증을 다음 단계로 제안한다.**

- 조성비 -> 염화물 활동도 -> 금속 화학종 **까지만** 주장한다.
- 화학종 -> 침출효율 연결은 **제언(future work)** 으로 넘긴다. 효율 수치를 말하지 않는다.
- 이유: 침출효율 실측이 없고, 주제문의 "성능 최적화"와 산출물 사이 격차가
  가장 큰 취약점이었음. 데이터를 늘리는 대신 주장 범위를 좁혀 격차를 없앤다.
- 발표/포스터 문구가 이 선을 넘는지 매번 확인할 것.

**미해결 TODO (중요도순)**
1. **침출효율 실측(ICP 등) 확보** — 현재 이 폴더엔 UV-Vis만 있고 y값이 없음. 계산으로 대체 불가
2. DFT 실행 (서버 접속 시. `DFT_PLAN.md` §6대로)
3. ~~`1;4` 표기가 ChCl:CA 순서인지 확인~~ → **2026-09-30 확정: ChCl:CA 순서가 맞음.**
   즉 `1:19` = ChCl 1 : 시트르산 19 (x_ChCl = 0.05), `19:1` = 염 많은 쪽(x_ChCl = 0.95).
   기존 해석·그림·원고 전부 이 배정이었으므로 수정 불필요.

## 포스터 원고 검증 이력 (2026-09-30)

`포스터_내파트_원고.md` 를 서로 다른 기준 3종으로 병렬 검증하고 전면 개정했다.
재검증할 때 같은 함정을 다시 밟지 않도록 **결론만** 남긴다.

- **원고 범위는 "R&D 섹션 전체"가 아니라 "그 안의 시뮬레이션 블록 하나"다.**
  배경·주제문·문헌 37편은 Introduction 담당, 실험 방법·UV-Vis 품질·DES 제조는 민이형 담당.
- **폐기된 논거 3개** (다시 쓰지 말 것):
  1. "Ni·Co 보정값 65.1 / 67.7 이 일치하니 검증됐다" — 각각 앵커 1점에 파라미터 1개라
     잔차가 구조적으로 0. 자유도 없음. 시트르산 결합상수도 Ni 5.4 / Co 5.0 으로 비슷해
     배위를 빠뜨려도 두 값이 똑같이 밀린다. **검증력 0.**
  2. "Co 흡광계수가 Ni 의 200배라 색은 Co 다" — 비교 상대인 '팔면체 Ni' 이 우리 결론이라
     **순환논증**. 200배의 출처도 없음(저장소 값은 150/3 = 50배).
  3. "5개 중 2개 맞출 확률 1/10" — 염화물 최다 쪽 두 개는 누구나 찍는다. 무작위 가정 불성립.
- **50 % 문턱으로는 x=0.75 가 음이온이 아니다** (T_Co = 68 °C > 실험 60 °C).
  그래서 "음이온이 된다"가 아니라 **"사면체 분율이 0 에서 유의하게 올라간다"** 로 말한다.
  60 °C 실측 5조성: 0 / 0 / 0 / **0.2424** / **0.9614** (`dft_results/ftd_60C.csv`).
- **Mn 전환 순서 주장은 삭제했다.** 고스핀 d⁵ 는 CFSE = 0 이라 교과서상 Td 전환이 가장 쉬운데
  우리 계산은 Co < Mn 이 나온다. Mn 은 미보정이므로 순서를 주장하지 않는다.
- **pH 는 소수점으로 인용하지 않는다.** `ph_estimate.py` 는 수용액 25 °C pKa·활동도 1 가정이라
  3~5 M·60 °C 에는 적용범위 밖. "1~2 수준"까지만. 단 결론("대부분 중성")은 pKa ±1 에도 견고.
- **시트르산 : 금속 비는 1000배가 아니라 약 100배다** (3.5 M vs 34 mM).
- 새로 만든 재현 스크립트: `dft_ni_complexes/ph_estimate.py`, `dft_ni_complexes/vial_color.py`.
  `anchored_map.py` 는 `ftd_60C.csv` · `shifts.csv` 를 저장하도록 패치했다(수치 근거 보존용).

동현의 역할: 파이프라인 총괄 + 계산화학 시뮬레이션(COSMO-RS 스크리닝) 담당.

## 계산 서버 환경 (핵심)

- **접속** (2026-09-24 정정): 실제 주소는 **`<SERVER_IP>`, 포트 `<PORT>`**, 계정 `navy030303`, 호스트명 `anode0` (**sudo 없음**).
  - 로컬 `~/.ssh/config`에 `anode0` 항목이 등록되어 있어 **`ssh anode0` 한 줄로 키 인증 접속**됨 (`~/.ssh/anode0_key`).
  - ⚠ MobaXterm 세션 이름이 `192.168.0.137`로 되어 있는데 **이건 세션 라벨일 뿐 실제 주소가 아님.** 집 공유기가 같은 사설대역(`192.168.0.x`)을 쓰면 엉뚱한 기기에 붙으므로 그 주소로 접속 시도하지 말 것.
  - 서버 홈은 `/home/users/navy030303` (`/home/navy030303` 아님).
- **⚠ Slurm에 `--mem` 쓰지 말 것**: 이 노드는 `RealMemory=1`로 설정되어 메모리를 관리하지 않음. `--mem`을 넣으면 `Memory specification can not be satisfied`로 **제출 자체가 거부됨**. 메모리 상한은 ORCA `%maxcore`로만 조절한다. 파티션 `long`의 TIMELIMIT은 `infinite`라 `-t`는 넉넉히 줘도 됨.
- 실제로는 **Docker 컨테이너**이며 (`/etc/resolv.conf`가 Docker Engine 생성), **완전히 인터넷이 막혀있음** (DNS도 안 됨, 프록시도 없음). 모든 설치는 로컬(Windows)에서 미리 받아 SFTP로 업로드하는 방식으로 진행했음.
- 하드웨어: CPU 128스레드(Xeon Gold 6530 ×2), RAM 251GB, GPU 2장(RTX 4090 24GB + RTX PRO 5000 Blackwell 48GB), 홈 디스크 쿼터 없음(654GB+ 여유)
- Slurm 단일 노드(`long` 파티션, 노드 자기 자신 하나뿐)로 잡 큐잉
- `/opt/anaconda3/2024.10`에 관리자가 깔아둔 공용 Anaconda(Python 3.12.2)가 있어 이걸 base로 venv를 만듦 (전체 인터넷 차단이라 conda 자체 설치/패키지 설치는 불가능해서 conda는 "이미 있는 python 바이너리 제공용"으로만 활용)
- 참고: 이 서버는 Materials Studio(BIOVIA, 다른 랩원이 씀 — `~/ms_jobs/`에 Forcite 고분자 가교 시뮬레이션 흔적 있음)와 Schrödinger Suite(`/opt/schrodinger/2026-1`, Jaguar/Desmond 등)도 이미 쓰이고 있는 공용 워크스테이션. COSMOtherm/TURBOMOLE/TmoleX는 이 서버에도, 동현 데스크탑의 Materials Studio 라이선스에도 **확인된 바 없음** → 오픈소스 대안(openCOSMO-RS)으로 확정.

## 지금까지 설치·검증 완료된 것

### 1. Python 가상환경
```bash
source ~/venvs/cosmors/bin/activate
```
- `/opt/anaconda3/2024.10/bin/python3` (3.12) 기반 `venv` (conda create 아님 — 오프라인이라 venv로 만듦)
- 설치된 패키지(전부 오프라인 wheel로 설치, `pip install --no-index --find-links=...`): rdkit, numpy, scipy, matplotlib, py3Dmol, pandas, plotly, scikit-image + 의존성 전체
- `opencosmorspy` (openCOSMO-RS_py, TUHH-TVT) — import명이 `opencosmors`가 아니라 **`opencosmorspy`**임에 주의
- 검증 명령: `python -c "import rdkit, opencosmorspy; print('OK')"` → 통과 확인됨

### 2. ORCA 6.1.1 (양자화학 엔진)
- 위치: `~/software/orca/`
- ORCA Forum에서 동현이 직접 계정 가입 후 개인 다운로드(재배포 금지 라이선스 — 다른 팀원도 필요하면 각자 받아야 함), MobaXterm SFTP로 업로드해서 `tar xf ... --strip-components=1`로 압축 해제
- 빌드: `orca_6_1_1_linux_x86-64_shared_openmpi418.tar.xz` (OpenMPI 4.1.8 전용 빌드)
- PATH 설정: `~/.bashrc`에 `export PATH="$HOME/software/orca:$PATH"` 및 `LD_LIBRARY_PATH` 추가됨 (영구 반영)

### 3. OpenMPI 4.1.8 (ORCA 병렬계산용)
- 시스템 모듈은 `ompi/5.0.6`으로 **ORCA가 요구하는 4.1.8과 버전이 안 맞아서** 소스에서 직접 컴파일
- 소스: `~/build/openmpi-4.1.8/`, 설치 위치: `~/software/openmpi/`
- 빌드는 `srun --cpus-per-task=8 --pty bash`로 세션 확보 후 `./configure --prefix=$HOME/software/openmpi && make -j8 && make install`로 진행 (안전원칙: 무거운 작업은 로그인쉘에서 직접 안 함)
- **전역 `.bashrc`엔 안 넣음** (시스템 openmpi 5.0.6과 충돌 방지, GROMACS 등 다른 모듈과의 간섭 방지). 대신 필요할 때만 아래로 활성화:
```bash
source ~/software/activate_orca_mpi.sh
# 내용: PATH/LD_LIBRARY_PATH에 ~/software/orca, ~/software/openmpi/bin,lib 추가
```
- **중요 설정**: `~/.openmpi/mca-params.conf`에 `rmaps_base_oversubscribe = true` 추가 필요 — 이게 없으면 Slurm 세션 안에서 mpirun이 "슬롯 부족" 에러로 즉시 실패함 (Slurm이 `--cpus-per-task`를 "슬롯 1개"로만 인식하는 문제, 오버서브스크라이브 허용으로 우회)

### 4. openCOSMO-RS_conformer_pipeline (RDKit→ORCA→COSMO 파일 생성 스크립트)
- 위치: `~/software/cosmors_offline/openCOSMO-RS_conformer_pipeline/`
- **저장소 자체에 버그 3곳 있었음(패치 완료)** — `find_balloon_executable`, ORCA 경로 탐색, `otool_xtb` 경로 탐색 함수에서 `try:` 블록에 `except`가 누락되어 있었음 (원본 GitHub 저장소의 문제, 우리 쪽 전송 오류 아님). 전부 `except Exception: pass`로 패치함.
- **아직 안 고친 버그**: `_setup_folder()` 함수가 "이미 100% 완료된 계산을 재개하려 할 때" `last_step_dir` UnboundLocalError로 죽음. **우회법**: 이미 완료된 분자를 같은 input 파일에 다시 넣지 말 것 — 새 분자는 항상 그 분자만 담긴 새 `.inp` 파일로 실행.
- 사용법 (탭 구분 4컬럼: name, SMILES, xyz파일(빈칸 가능), charge):
```bash
printf "이름\tSMILES\t\t0\n" > my.inp
python ConformerGenerator.py --structures_file my.inp --cpcm_radii_file cpcm_radii.inp --n_cores 1
```
- 방법론: RDKit 콘포머 생성 → ORCA XTB2(ALPB) 예비 → BP86/def2-TZVP(-f) → BP86/def2-TZVP → CPCM_fast → CPCM_final. **BP86/def2-TZVP는 지난 문헌조사(Zurob/Abranches/Amusa 등)의 상용 TmoleX/COSMOthermX 워크플로우와 동일 이론수준** — 방법론 연속성 확보됨.
- 출력: `<분자명>/COSMO_TZVPD/<분자명>_c000.orcacosmo` 파일이 opencosmorspy의 입력이 됨

### 5. openCOSMO-RS_py (활동도계수 계산)
- 위치: `~/software/cosmors_offline/openCOSMO-RS_py/` (pip로 이미 venv에 설치되어 있음, 소스는 참고/예제용)
- 사용 예시:
```python
import numpy as np
from opencosmorspy import COSMORS
crs = COSMORS(par='default_orca')
crs.add_molecule(['water/COSMO_TZVPD/water_c000.orcacosmo'])
crs.add_molecule(['methanol/COSMO_TZVPD/methanol_c000.orcacosmo'])
crs.add_job(np.array([0.5, 0.5]), 298.15, refst='pure_component')
results = crs.calculate()
print(results['tot']['lng'])  # 전체 ln(활동도계수)
```
- 실제 검증 완료: water-methanol 이성분계 결과가 물리적으로 타당함(양의 편차, 문헌과 정성적으로 일치)

## 핵심 실측 데이터 및 결론

| 분자 | 코어수 | 소요시간 |
|---|---|---|
| water | 1 | 20.6초 |
| methanol | 1 | 54.4초 |
| ethanol | 2 (MPI) | 141.3초 (더 느려짐!) |

**결론: 우리 프로젝트 분자들(작은 유기분자)은 MPI 병렬화(`n_cores>1`)가 오히려 손해다.** 통신 오버헤드가 계산량보다 커짐. → **`n_cores=1`(serial)을 기본으로 하고, 여러 분자를 Slurm 잡 배열(`--array=1-N%64`)로 동시에 여러 개 띄우는 "태스크 병렬화" 방식을 쓸 것.** (OpenMPI 설치 자체는 완료되어 있으니 나중에 훨씬 큰 분자를 다루게 되면 다시 꺼내 쓸 수 있음.)

## 다음 단계 (아직 안 한 것)

1. 지난 문헌조사에서 확정한 실제 HBA/HBD 후보 분자 리스트(`화학물질리스트.xlsx` 확인 필요)로 `molecules.inp` 작성
2. Slurm 잡 배열로 전체 스크리닝 실행 (`%64` 동시실행 제한, 안전원칙 준수)
3. 결과 `.orcacosmo` 파일들을 모아 `opencosmorspy`로 몰비 스윕 + 활동도계수 계산 → 영진 형 ML 파이프라인과 연동

## 서버 안전 원칙 (계획서 원문, `warm-hugging-whale.md` 참고)

1. 전부 `$HOME` 내부에만 설치 (sudo 없음, 시스템 경로 접근 안 함)
2. 무거운 계산은 반드시 `sbatch`/`srun`으로 제출 (로그인쉘에서 직접 무거운 잡 돌리지 않기)
3. 코어 수 명시적 제한, 실행 전 `squeue`로 다른 사용자 작업 확인
4. ORCA 라이선스 준수 — tar 파일 공유 금지, 각자 개인 다운로드
5. 문제 생기면 `rm -rf ~/venvs/cosmors ~/software/orca ~/software/openmpi ~/software/cosmors_offline`로 완전 원복 가능 (전부 격리 설치라 시스템에 영향 없음)

## 로컬(Windows) 쪽 참고

- `literature/offline_bundle/`에 서버로 전송했던 wheel/저장소 로컬 사본이 남아있음 (재전송 필요시 재사용 가능)
- `literature/` 폴더에 원문 확보한 PDF 다수 (Zurob, Amusa, Qin, Velez&Acevedo, Lemaoui, Almashjary, Naseem, Kollau, Martins, Alizadeh 등) — 전부 `DES_조성비_최적화_선행연구_종합보고서.docx`에 검증상태 반영 완료됨
