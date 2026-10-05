#!/bin/bash
#SBATCH -J cofreq
#SBATCH -p long
#SBATCH --array=1-4%4
#SBATCH -c 1
#SBATCH -t 0-06:00:00
# --mem 금지 (RealMemory=1)
#SBATCH -o logs/slurm-%A_%a.out
#SBATCH --requeue

cd "$SLURM_SUBMIT_DIR" || exit 1
[ -f ~/software/activate_orca_mpi.sh ] && source ~/software/activate_orca_mpi.sh
N=$(sed -n "${SLURM_ARRAY_TASK_ID}p" species.txt)
[ -z "$N" ] && exit 0
say(){ echo "[$(date '+%H:%M:%S')] $N $*" >> STATUS.txt; }
command -v orca >/dev/null || { say "FATAL orca 없음"; exit 1; }

run(){ local s=$1 t=$2
  [ -f "${N}.$s.done" ] && return 0
  local t0=$SECONDS
  timeout -k 30 "$t" orca "${N}_${s}.inp" > "${N}_${s}.out" 2>&1
  if grep -q "ORCA TERMINATED NORMALLY" "${N}_${s}.out"; then
    say "$s OK $((SECONDS-t0))s"; touch "${N}.$s.done"; rm -f "${N}_${s}"*.tmp*; return 0
  fi
  say "$s 실패 $((SECONDS-t0))s"; rm -f "${N}_${s}"*.tmp*; return 1
}

cp "${N}_start.xyz" "${N}_geo.xyz"
run opt 5400 && [ -f "${N}_opt.xyz" ] || cp "${N}_geo.xyz" "${N}_opt.xyz"
run freq 10800
say DONE
