#!/bin/bash
#SBATCH -J ni_night
#SBATCH -p long
#SBATCH -c 1
#SBATCH -t 12:00:00
# --mem 금지 (이 노드는 RealMemory=1). 이 잡은 대부분 sleep 이라 부하도 거의 없다.
#SBATCH -o logs/overnight-%j.out
#SBATCH --requeue
#
# overnight.py 를 Slurm 안에서 돌린다. 이유:
#   - SSH 세션이 끊겨도 살아남는다 (nohup 과 달리 스케줄러가 관리)
#   - --requeue 로 노드 장애 시 자동 재시작. overnight.py 는 현재 상태를 다시
#     점검하고 이어가므로 재시작해도 안전하다.
# 감독 스크립트가 sbatch 로 재계산을 다시 던지므로, 이 잡 자체는 1코어만 쓴다.

cd "$SLURM_SUBMIT_DIR" || exit 1
source ~/venvs/cosmors/bin/activate 2>/dev/null
[ -f ~/software/activate_orca_mpi.sh ] && source ~/software/activate_orca_mpi.sh

python overnight.py
exit 0          # 감독이 무엇을 보고하든 잡 자체는 성공으로 끝낸다
