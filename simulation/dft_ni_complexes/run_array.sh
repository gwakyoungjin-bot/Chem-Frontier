#!/bin/bash
#SBATCH -J ni_dft
#SBATCH -p long
#SBATCH --array=1-7%7
#SBATCH -c 1
#SBATCH -t 2-00:00:00
# --mem 은 쓰지 않는다: 이 서버 Slurm 은 RealMemory=1 이라 제출이 거부된다.
# 메모리 상한은 ORCA 의 %maxcore(6000 MB) 가 결정한다.
# 파티션 TIMELIMIT 은 infinite 이므로 -t 2일은 안전하게 받아들여진다.
#SBATCH -o logs/slurm-%A_%a.out
#SBATCH --requeue
#
# 무인 실행 안전판. 원칙: 어떤 단계가 실패해도 잡을 죽이지 않고 다음으로 넘어간다.
#   - set -e 를 쓰지 않는다 (의도적). 실패는 기록하고 계속한다.
#   - 단계별 .done 플래그 -> 재제출하면 끝난 단계는 건너뛴다 (--requeue 대비).
#   - 단계마다 재시도 사다리: 기본 -> 수렴옵션 완화 -> 최후수단.
#   - 기하가 하나도 없으면 직전 단계 결과로 대체한다. 절대 빈손으로 끝내지 않는다.
# 안전원칙(CLAUDE.md): serial(코어1) x 잡배열. 작은 분자는 MPI가 오히려 느림.

cd "$SLURM_SUBMIT_DIR" || exit 1
[ -f ~/software/activate_orca_mpi.sh ] && source ~/software/activate_orca_mpi.sh

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" species.txt)
N=${LINE%% *}
STEPS=${LINE#* }
[ -z "$N" ] && exit 0

say() { echo "[$(date '+%m-%d %H:%M')] $N $*" >> STATUS.txt; }

# orca 가 PATH 에 없으면 즉시 기록하고 종료 (조용히 실패하는 것만 막는다)
command -v orca >/dev/null || { say "FATAL orca not on PATH"; exit 1; }

# 디스크 여유 확인 (NEVPT2 가 수 GB 쓴다)
FREE=$(df -Pk . | awk 'NR==2{print int($4/1048576)}')
[ "$FREE" -lt 20 ] && say "WARN disk only ${FREE}GB free"

# 단계별 시간상한(초). 한 단계가 멎어서 walltime 을 다 먹고 뒤 단계를 굶기는 것을 막는다.
# 최악의 경우 합계가 Slurm -t(2일=172800s) 아래가 되도록 잡았다:
#   xtb 2x1800 + opt 4x14400 + freq 2x32400 + dd 4x9000 = 162000s < 172800s
# 정상 경로는 훨씬 짧으므로(총 7~15h) 이 상한은 병적인 경우에만 걸린다.
TMAX=0

# $1=단계이름 $2=인풋베이스 $3...=추가키워드(재시도 사다리)
# 성공하면 0, 사다리 전부 실패하면 1. 어느 쪽이든 호출부는 계속 진행한다.
attempt() {
  local step=$1 base=$2; shift 2
  [ -f "$base.inp" ] || { say "$step SKIP (no $base.inp)"; return 1; }
  local try=0
  for extra in "$@"; do
    try=$((try+1))
    cp "$base.inp" "${base}_run.inp"
    [ -n "$extra" ] && sed -i "1i ! $extra" "${base}_run.inp"
    local t0=$SECONDS
    timeout -k 60 "$TMAX" orca "${base}_run.inp" > "${base}.out" 2>&1
    [ $? -eq 124 ] && say "$step TIMEOUT try=$try (상한 ${TMAX}s 초과 -> 중단하고 다음으로)"
    local dt=$((SECONDS-t0))
    if grep -q "ORCA TERMINATED NORMALLY" "${base}.out"; then
      say "$step OK try=$try ${dt}s ${extra:+[$extra]}"
      cp "${base}_run.inp" "${base}.used.inp"
      rm -f "${base}"*.tmp*
      return 0
    fi
    say "$step FAIL try=$try ${dt}s ${extra:+[$extra]}"
    # 실패해도 도중 기하가 남아있으면 다음 시도의 출발점으로 승계
    [ -f "${base}_run.xyz" ] && cp "${base}_run.xyz" "${N}_geo.xyz"
    rm -f "${base}"*.tmp*
  done
  return 1
}

has() { case ",$STEPS," in *",$1,"*) return 0;; *) return 1;; esac; }

# ---- 시작 구조 배치 -------------------------------------------------
[ -f "${N}_geo.xyz" ] || cp "${N}_start.xyz" "${N}_geo.xyz"

# ---- 0단계: XTB2 사전최적화 (실패해도 무해, 원본 구조로 계속) -------
if has xtb && [ ! -f "${N}.xtb.done" ]; then
  TMAX=1800
  if attempt xtb "${N}_xtb" "" "SlowConv"; then
    [ -f "${N}_xtb_run.xyz" ] && cp "${N}_xtb_run.xyz" "${N}_geo.xyz"
    touch "${N}.xtb.done"
  else
    say "xtb 실패 -> 손으로 만든 초기구조 그대로 사용 (계속 진행)"
  fi
fi

# ---- 단일원자(Cl-)는 단일점만 -------------------------------------
if has sp; then
  TMAX=1800
  [ -f "${N}.sp.done" ] || { attempt sp "${N}_sp" "" "SlowConv" && touch "${N}.sp.done"; }
  say "DONE (sp only)"; exit 0
fi

# ---- 1단계: 기하최적화 --------------------------------------------
# 사다리: TightOpt -> SCF 완화 -> 수렴기준 완화 -> 최후수단
if [ ! -f "${N}.opt.done" ]; then
  TMAX=14400
  attempt opt "${N}_opt" "" "SlowConv" "SlowConv NormalOpt" "VerySlowConv NormalOpt KDIIS"
  if [ -f "${N}_opt_run.xyz" ]; then
    cp "${N}_opt_run.xyz" "${N}_opt.xyz"
    [ -f "${N}_opt_run.gbw" ] && cp "${N}_opt_run.gbw" "${N}_opt.gbw"
    grep -q "ORCA TERMINATED NORMALLY" "${N}_opt.out" && touch "${N}.opt.done" \
      || say "opt 미수렴 -> 마지막 기하를 그대로 사용 (한계로 기록할 것)"
  else
    cp "${N}_geo.xyz" "${N}_opt.xyz"
    say "opt 결과 없음 -> XTB 기하로 대체 (한계로 기록할 것)"
  fi
fi

# ---- 2단계: 진동수 -> 열역학 G ------------------------------------
# 여기가 제일 오래 걸린다(수치 Hessian = 6N+1회 gradient). 실패해도 dd 는 진행.
if has freq && [ ! -f "${N}.freq.done" ]; then
  TMAX=32400
  attempt freq "${N}_freq" "" "SlowConv" && touch "${N}.freq.done" \
    || say "freq 실패 -> 이 종은 dG 없음. 나머지 종만으로 사다리 일부 구성 가능"
fi

# ---- 3단계: d-d 전이 (NiCl4, Ni_aq6 만) ---------------------------
# 1차: MORead + 삼중항10/단일항15.  2차: MORead 없이 삼중항 10개만.
if has dd && [ ! -f "${N}.dd.done" ]; then
  TMAX=9000
  if attempt dd "${N}_dd" "" "SlowConv"; then
    touch "${N}.dd.done"
  elif attempt dd2 "${N}_dd2" "" "SlowConv"; then
    say "dd 축소판 성공(삼중항만) - 가시영역 밴드는 전부 포함됨"
    cp "${N}_dd2.out" "${N}_dd.out"
    touch "${N}.dd.done"
  else
    say "dd 실패 -> eps 예측 없음. dG 결과는 영향 없음"
  fi
fi

say "DONE"
