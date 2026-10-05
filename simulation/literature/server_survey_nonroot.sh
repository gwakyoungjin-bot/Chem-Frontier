#!/bin/bash
# COSMO-RS/QM 계산 SW 조사 스크립트 (sudo 없는 일반 계정용)
# 실행: bash server_survey_nonroot.sh 2>&1 | tee ~/cosmo_survey_result.txt

echo "===== [0] 홈 디렉토리 내 계산화학 관련 흔적 확인 (사이드바에서 이미 보인 폴더들) ====="
echo "-- .schrodinger 폴더 --"
ls -la ~/.schrodinger 2>/dev/null
find ~/.schrodinger -maxdepth 2 2>/dev/null
echo "-- ms_jobs / ms_test 폴더 (Materials Studio 잡 폴더로 추정) --"
ls -la ~/ms_jobs ~/ms_test 2>/dev/null
find ~/ms_jobs ~/ms_test -maxdepth 2 2>/dev/null | head -30
echo "-- .scripts 폴더 --"
ls -la ~/.scripts 2>/dev/null
echo "-- .pbs_spool (PBS 스케줄러 사용 흔적) --"
ls -la ~/.pbs_spool 2>/dev/null
echo "-- .bash_profile / .bashrc 내 module load / 계산SW 관련 라인 --"
grep -iE "module|cosmo|turbo|schrodinger|gaussian|orca|nwchem|gamess|materials|SCHRODINGER|TURBODIR" ~/.bash_profile ~/.bashrc 2>/dev/null

echo ""
echo "===== [1] Schrodinger Suite 확인 (홈에 .schrodinger 있으므로 최우선 확인) ====="
echo "SCHRODINGER=$SCHRODINGER"
which schrodinger jaguar desmond maestro 2>/dev/null
[ -n "$SCHRODINGER" ] && ls -la "$SCHRODINGER" 2>/dev/null

echo ""
echo "===== [2] COSMOtherm/TURBOMOLE/TmoleX 표준 경로 확인 (world-readable 여부만 확인 가능) ====="
for d in /opt/COSMOlogic* /opt/COSMOthermX* /opt/turbomole* /opt/TURBOMOLE* \
         /opt/biovia* /opt/BIOVIA* /usr/local/COSMOlogic* /usr/local/turbomole* \
         /apps/COSMOlogic* /apps/turbomole* /software/COSMOlogic* /software/turbomole* \
         /shared/COSMOlogic* /shared/turbomole*; do
  [ -e "$d" ] && echo "FOUND: $d" && ls -la "$d" 2>/dev/null
done
echo "TURBODIR=$TURBODIR"
echo "COSMOTHERM_DIR=$COSMOTHERM_DIR"
which ridft dscf jobex TmoleX cosmotherm COSMOtherm 2>/dev/null
env | grep -iE "cosmo|turbo|biovia"

echo ""
echo "===== [3] HPC 모듈 시스템 (Lmod/environment modules) 확인 ====="
if command -v module >/dev/null 2>&1; then
  module avail 2>&1
else
  echo "module 명령 없음"
fi

echo ""
echo "===== [4] 잡 스케줄러 확인 (PBS 흔적 있음 → 상세 확인) ====="
command -v qstat >/dev/null 2>&1 && { echo "-- qstat -Q (큐 목록) --"; qstat -Q; echo "-- qstat -f (내 잡) --"; qstat -f 2>/dev/null | head -50; }
command -v pbsnodes >/dev/null 2>&1 && { echo "-- pbsnodes -a (노드/자원) --"; pbsnodes -a 2>/dev/null | head -80; }
command -v sinfo >/dev/null 2>&1 && { echo "-- Slurm sinfo --"; sinfo; }

echo ""
echo "===== [5] 대체/경쟁 QM 엔진 확인 (Gaussian, ORCA, NWChem, GAMESS) ====="
which g16 g09 orca nwchem gamess.x rungms 2>/dev/null

echo ""
echo "===== [6] 오픈소스 COSMO 대안(파이썬 패키지) 설치 여부 확인 ====="
pip list 2>/dev/null | grep -iE "cosmo"
conda list 2>/dev/null | grep -iE "cosmo"
python3 -c "import cosmo" 2>&1 | head -1

echo ""
echo "===== [7] 접근 가능한 범위 내 파일시스템 검색 (권한 없는 곳은 조용히 건너뜀) ====="
find /opt /usr/local /apps /software /shared "$HOME" -maxdepth 4 \
  \( -iname "*cosmo*therm*" -o -iname "*tmolex*" -o -iname "*turbomole*" -o -iname "*opencosmo*" -o -iname "*cosmosac*" \) \
  2>/dev/null | head -50

echo ""
echo "===== [8] 하드웨어 자원 확인 (계산 자원 판단용, 일반계정에서도 조회 가능) ====="
echo "-- CPU --"; lscpu 2>/dev/null | grep -E "Model name|CPU\(s\)|Socket|Thread" || cat /proc/cpuinfo | grep -m1 "model name"
echo "-- MEM --"; free -h
echo "-- GPU --"; command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi --query-gpu=name,memory.total --format=csv || echo "GPU 없음/nvidia-smi 없음"
echo "-- 홈 디스크 쿼터 --"; df -h "$HOME" 2>/dev/null; quota -s 2>/dev/null

echo ""
echo "===== 조사 완료. 출력 전체를 저에게 붙여넣어 주세요. ====="
