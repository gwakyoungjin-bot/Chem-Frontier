#!/bin/bash
#SBATCH -J mcl3
#SBATCH -p long
#SBATCH --array=1-6%6
#SBATCH -c 1
#SBATCH -t 2-00:00:00
# --mem 금지: 이 노드는 RealMemory=1 이라 --mem 을 넣으면 제출이 거부된다.
# 메모리 상한은 ORCA 의 %maxcore(6000 MB) 가 결정한다.
#SBATCH -o logs/slurm-%A_%a.out
#SBATCH --requeue

cd "$SLURM_SUBMIT_DIR" || exit 1
[ -f ~/software/activate_orca_mpi.sh ] && source ~/software/activate_orca_mpi.sh
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" species.txt)
N=${LINE%% *}; STEPS=${LINE#* }
[ -z "$N" ] && exit 0
say() { echo "[$(date '+%m-%d %H:%M')] $N $*" >> STATUS.txt; }
command -v orca >/dev/null || { say "FATAL orca not on PATH"; exit 1; }

TMAX=0
attempt() {
  local step=$1 base=$2; shift 2
  [ -f "$base.inp" ] || { say "$step SKIP"; return 1; }
  local try=0
  for extra in "$@"; do
    try=$((try+1))
    cp "$base.inp" "${base}_run.inp"
    [ -n "$extra" ] && sed -i "1i ! $extra" "${base}_run.inp"
    local t0=$SECONDS
    timeout -k 60 "$TMAX" orca "${base}_run.inp" > "${base}.out" 2>&1
    [ $? -eq 124 ] && say "$step TIMEOUT try=$try"
    local dt=$((SECONDS-t0))
    if grep -q "ORCA TERMINATED NORMALLY" "${base}.out"; then
      say "$step OK try=$try ${dt}s ${extra:+[$extra]}"; rm -f "${base}"*.tmp*; return 0
    fi
    say "$step FAIL try=$try ${dt}s ${extra:+[$extra]}"
    # 실패해도 마지막 기하를 물려주면 다음 시도가 더 가까운 곳에서 출발한다
    [ -f "${base}_run.xyz" ] && cp "${base}_run.xyz" "${N}_geo.xyz"
    rm -f "${base}"*.tmp*
  done
  return 1
}

[ -f "${N}_geo.xyz" ] || cp "${N}_start.xyz" "${N}_geo.xyz"

TMAX=1800
if [ ! -f "${N}.xtb.done" ]; then
  attempt xtb "${N}_xtb" "" "SlowConv" \
    && { [ -f "${N}_xtb_run.xyz" ] && cp "${N}_xtb_run.xyz" "${N}_geo.xyz"; touch "${N}.xtb.done"; } \
    || say "xtb 실패 -> 초기구조로 진행"
fi

TMAX=21600
if [ ! -f "${N}.opt.done" ]; then
  attempt opt "${N}_opt" "" "SlowConv" "SlowConv NormalOpt" "VerySlowConv NormalOpt KDIIS"
  if [ -f "${N}_opt_run.xyz" ]; then
    cp "${N}_opt_run.xyz" "${N}_opt.xyz"
    [ -f "${N}_opt_run.gbw" ] && cp "${N}_opt_run.gbw" "${N}_opt.gbw"
    grep -q "ORCA TERMINATED NORMALLY" "${N}_opt.out" && touch "${N}.opt.done" \
      || say "opt 미수렴 -> 마지막 기하 사용"
  else
    cp "${N}_geo.xyz" "${N}_opt.xyz"; say "opt 결과 없음"
  fi
fi

TMAX=36000
if [ ! -f "${N}.freq.done" ]; then
  attempt freq "${N}_freq" "" "SlowConv" && touch "${N}.freq.done" || say "freq 실패"
fi
say "DONE"
