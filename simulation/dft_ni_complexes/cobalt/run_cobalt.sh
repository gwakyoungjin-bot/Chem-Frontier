#!/bin/bash
#SBATCH -J co_dft
#SBATCH -p long
#SBATCH --array=1-5%5
#SBATCH -c 1
#SBATCH -t 2-00:00:00
# --mem 금지 (RealMemory=1). 메모리는 %maxcore 로만 제한.
#SBATCH -o logs/slurm-%A_%a.out
#SBATCH --requeue
# 사다리 구조는 ../run_array.sh 와 동일. 실행 중 파일을 건드리지 않으려 분리했다.

cd "$SLURM_SUBMIT_DIR" || exit 1
[ -f ~/software/activate_orca_mpi.sh ] && source ~/software/activate_orca_mpi.sh
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" species.txt)
N=${LINE%% *}; STEPS=${LINE#* }
[ -z "$N" ] && exit 0
say() { echo "[$(date '+%m-%d %H:%M')] $N $*" >> STATUS.txt; }
command -v orca >/dev/null || { say "FATAL orca not on PATH"; exit 1; }
has() { case ",$STEPS," in *",$1,"*) return 0;; *) return 1;; esac; }

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

TMAX=9000
if has dd && [ ! -f "${N}.dd.done" ]; then
  if attempt dd "${N}_dd" "" "SlowConv"; then touch "${N}.dd.done"
  elif attempt dd2 "${N}_dd2" "" "SlowConv"; then
    cp "${N}_dd2.out" "${N}_dd.out"; touch "${N}.dd.done"; say "dd MORead 없이 성공"
  else say "dd 실패 -> eps 없음 (dG 는 영향 없음)"; fi
fi
say "DONE"
