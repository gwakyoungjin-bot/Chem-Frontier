#!/bin/bash
#SBATCH -J co_watch
#SBATCH -p long
#SBATCH -c 1
#SBATCH -t 1-00:00:00
# --mem 금지 (RealMemory=1). 이 잡은 대부분 sleep 이라 부하가 거의 없다.
#SBATCH -o logs/watch-%j.out
#SBATCH --requeue
#
# co_dft 가 끝나면 검증 -> (필요시)재시도 -> 후처리까지 자동으로 이어간다.
# --requeue 로 노드 장애 시 재시작돼도 watch_cobalt.py 가 현재 상태를 다시 보고 이어간다.

cd "$SLURM_SUBMIT_DIR" || exit 1
source ~/venvs/cosmors/bin/activate 2>/dev/null
[ -f ~/software/activate_orca_mpi.sh ] && source ~/software/activate_orca_mpi.sh
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4

python watch_cobalt.py
exit 0
