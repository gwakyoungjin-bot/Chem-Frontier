#!/bin/bash
# COSMO-RS/QM 계산 SW 서버 조사 스크립트
# 실행: sudo bash server_survey.sh 2>&1 | tee ~/cosmo_survey_result.txt

echo "===== [1] 표준 설치 경로 확인 (COSMOlogic/BIOVIA/TURBOMOLE 공식 문서 기준 경로) ====="
for d in /opt/COSMOlogic* /opt/COSMOthermX* /opt/turbomole* /opt/TURBOMOLE* \
         /opt/biovia* /opt/BIOVIA* /usr/local/COSMOlogic* /usr/local/turbomole* \
         /apps/COSMOlogic* /apps/turbomole* /software/COSMOlogic* /software/turbomole*; do
  [ -e "$d" ] && echo "FOUND: $d" && ls -la "$d" 2>/dev/null
done

echo ""
echo "===== [2] PATH / 환경변수 확인 ====="
echo "TURBODIR=$TURBODIR"
echo "COSMOTHERM_DIR=$COSMOTHERM_DIR"
which ridft dscf jobex TmoleX cosmotherm COSMOtherm 2>/dev/null
env | grep -iE "cosmo|turbo|biovia"

echo ""
echo "===== [3] HPC 모듈 시스템 (Lmod/environment modules) 확인 ====="
if command -v module >/dev/null 2>&1; then
  module avail 2>&1 | grep -iE "cosmo|turbomole|gaussian|orca|nwchem|gamess|materials"
else
  echo "module 명령 없음 (Lmod 미사용 환경일 수 있음)"
fi

echo ""
echo "===== [4] 잡 스케줄러 확인 (Slurm/PBS/SGE - 계산자원 할당량 파악용) ====="
command -v sinfo >/dev/null 2>&1 && echo "Slurm 감지:" && sinfo && squeue -u $USER
command -v qstat >/dev/null 2>&1 && echo "PBS/SGE 감지:" && qstat -Q

echo ""
echo "===== [5] FlexLM 라이선스 서버 확인 (BIOVIA 제품군 표준 라이선싱) ====="
ps aux | grep -iE "lmgrd|lmadmin|flex" | grep -v grep
find / -maxdepth 6 -iname "license.dat" -o -iname "*.lic" 2>/dev/null | grep -v "^find:"
sudo netstat -tlnp 2>/dev/null | grep -E ":270(0[0-9])" # FlexLM 기본 포트대 27000-27009

echo ""
echo "===== [6] 대체/경쟁 QM 엔진 확인 (Gaussian, ORCA, NWChem, GAMESS) ====="
which g16 g09 orca nwchem gamess.x rungms 2>/dev/null
find / -maxdepth 5 \( -iname "*gaussian*" -o -iname "*orca*" -o -iname "*nwchem*" -o -iname "*gamess*" \) -type d 2>/dev/null | grep -v "^find:"

echo ""
echo "===== [7] 오픈소스 COSMO 대안 설치 여부 확인 ====="
pip list 2>/dev/null | grep -iE "cosmo"
conda list 2>/dev/null | grep -iE "cosmo"
find / -maxdepth 6 -iname "*opencosmo*" -o -iname "*cosmosac*" -o -iname "*cosmo-sac*" 2>/dev/null | grep -v "^find:"

echo ""
echo "===== [8] 전체 파일시스템 넓게 재검색 (sudo 필요, 시간 걸림) ====="
sudo find / -iname "*cosmo*therm*" -o -iname "*tmolex*" 2>/dev/null | grep -v "^find:" | head -50

echo ""
echo "===== [9] 하드웨어 자원 확인 (계산 자원 할당량 판단용) ====="
echo "-- CPU --"; lscpu | grep -E "Model name|CPU\(s\)|Socket|Thread"
echo "-- MEM --"; free -h
echo "-- GPU --"; command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi --query-gpu=name,memory.total --format=csv || echo "GPU 없음/nvidia-smi 없음"
echo "-- DISK --"; df -h /opt /home 2>/dev/null

echo ""
echo "===== 조사 완료. cosmo_survey_result.txt 파일을 저에게 다시 붙여넣어 주세요. ====="
