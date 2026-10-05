#!/bin/bash
#SBATCH -J mixed
#SBATCH -p long
#SBATCH --array=1-3%3
#SBATCH -c 1
#SBATCH -t 0-04:00:00
# --mem 금지 (RealMemory=1). 메모리는 %maxcore 로만.
#SBATCH -o logs/slurm-%A_%a.out
#SBATCH --requeue

cd "$SLURM_SUBMIT_DIR" || exit 1
[ -f ~/software/activate_orca_mpi.sh ] && source ~/software/activate_orca_mpi.sh
N=$(sed -n "${SLURM_ARRAY_TASK_ID}p" species.txt)
[ -z "$N" ] && exit 0
say() { echo "[$(date '+%H:%M:%S')] $N $*" >> STATUS.txt; }
command -v orca >/dev/null || { say "FATAL orca 없음"; exit 1; }

run() {  # run <단계> <타임아웃초>
  local s=$1 t=$2
  [ -f "${N}.$s.done" ] && return 0
  local t0=$SECONDS
  timeout -k 30 "$t" orca "${N}_${s}.inp" > "${N}_${s}.out" 2>&1
  local dt=$((SECONDS-t0))
  if grep -q "ORCA TERMINATED NORMALLY" "${N}_${s}.out"; then
    say "$s OK ${dt}s"; touch "${N}.$s.done"; rm -f "${N}_${s}"*.tmp*; return 0
  fi
  say "$s 미완 ${dt}s (시간초과 또는 실패)"; rm -f "${N}_${s}"*.tmp*; return 1
}

cp "${N}_start.xyz" "${N}_geo.xyz"

# 1) xTB 로 빠르게 다듬는다
run xtb 600 && [ -f "${N}_xtb.xyz" ] && cp "${N}_xtb.xyz" "${N}_geo.xyz"

# 2) xTB 기하로 먼저 SP — 여기서 이미 '답' 이 하나 확보된다.
#    뒤 단계가 시간에 걸려 죽어도 빈손이 아니게 하는 장치.
cp "${N}_geo.xyz" "${N}_sp.xyz"
run sp 1200 && cp "${N}_sp.out" "${N}_sp_xtbgeom.out" && rm -f "${N}.sp.done"

# 3) DFT 로 느슨하게 최적화 (1시간 상한)
if run opt 3600; then
  [ -f "${N}_opt.xyz" ] && cp "${N}_opt.xyz" "${N}_sp.xyz"
else
  # 미수렴이어도 마지막 기하가 xTB 보다는 낫다
  [ -f "${N}_opt.xyz" ] && cp "${N}_opt.xyz" "${N}_sp.xyz" && say "미수렴 기하 사용"
fi

# 4) 개선된 기하로 다시 SP
run sp 1200

# 5) UV-Vis (345 nm 확인). 실패해도 에너지 결론에는 영향 없음.
run tda 2400 || say "tda 생략 -> UV-Vis 비교 불가"
say DONE
