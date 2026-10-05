"""Ni(II) 아쿠아/클로로 착물 ORCA 인풋 생성 (무인 실행 안전판, 2026-09-23).

목적이 바뀌었음: UV-Vis에 Ni 신호가 없으므로 "d-d 밴드 재현"이 아니라
  (1) 화학종 안정성 dG -> 조성-온도 가이드라인 지도
  (2) 흡광계수 eps 예측 -> "안 보인 게 정상"의 정량 증명
따라서 Freq(열역학)가 필수가 되고, d-d(NEVPT2)는 양 끝 2종만 돌린다.

토요일까지 무인 실행이므로 속도보다 "중간에 안 죽는 것"이 우선.
단계를 4개로 쪼개서 하나가 실패해도 뒤 단계가 이어지게 한다:
  xtb(사전최적화) -> opt -> freq -> dd
"""
import math, os, textwrap

R_NI_O, R_NI_CL_OH, R_NI_CL_TD = 2.06, 2.40, 2.27   # Angstrom, 실험 결정구조 전형값
R_OH, ANG_HOH = 0.97, 104.5

def water(u, d):
    """Ni에서 방향 u(단위벡터), 거리 d에 물 배위. O 고립전자쌍이 Ni를 향하도록 H를 바깥으로."""
    ux, uy, uz = u
    O = (ux*d, uy*d, uz*d)
    p = (1., 0., 0.) if abs(ux) < 0.9 else (0., 1., 0.)
    vx, vy, vz = p[1]*uz-p[2]*uy, p[2]*ux-p[0]*uz, p[0]*uy-p[1]*ux
    n = math.sqrt(vx*vx+vy*vy+vz*vz); vx, vy, vz = vx/n, vy/n, vz/n
    h = math.radians(ANG_HOH/2)
    out = []
    for s in (+1, -1):
        hx = ux*math.cos(h) + s*vx*math.sin(h)
        hy = uy*math.cos(h) + s*vy*math.sin(h)
        hz = uz*math.cos(h) + s*vz*math.sin(h)
        out.append(('H', (O[0]+hx*R_OH, O[1]+hy*R_OH, O[2]+hz*R_OH)))
    return [('O', O)] + out

OCT = [(1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1)]
s3 = 1/math.sqrt(3)
TET = [(s3,s3,s3), (s3,-s3,-s3), (-s3,s3,-s3), (-s3,-s3,s3)]

def build(dirs, ligs, d_cl):
    at = [('Ni', (0., 0., 0.))]
    for u, L in zip(dirs, ligs):
        if L == 'Cl':
            at.append(('Cl', (u[0]*d_cl, u[1]*d_cl, u[2]*d_cl)))
        else:
            at += water(u, R_NI_O)
    return at

# (이름, 전하, 스핀다중도, 원자리스트, 단계)
# Ni(II) d8 -> S=1 (mult 3).  H2O/Cl- 는 닫힌껍질 (mult 1).
# 참조종 H2O/Cl- 가 필요한 이유: 배위 교환 반응식의 양변을 맞춰야 dG가 나온다.
#   NiCl_n(H2O)_{6-n} + Cl- -> NiCl_{n+1}(H2O)_{m} + k H2O
# Td 종(NiCl3_aq1, NiCl4)은 물 개수가 한 번에 여러 개 빠지므로 계수가 일정하지 않다
# -> 문헌 앵커 하나로 상수 하나만 잡는 방식이 성립하려면 두 참조종 다 필요.
SPECIES = [
    ('NiCl4',     -2, 3, build(TET, ['Cl']*4,              R_NI_CL_TD), 'xtb,opt,freq,dd'),
    ('Ni_aq6',     2, 3, build(OCT, ['O']*6,               R_NI_CL_OH), 'xtb,opt,freq,dd'),
    ('NiCl2_aq4',  0, 3, build(OCT, ['Cl','Cl']+['O']*4,   R_NI_CL_OH), 'xtb,opt,freq'),
    ('NiCl3_aq1', -1, 3, build(TET, ['Cl','Cl','Cl','O'],  R_NI_CL_TD), 'xtb,opt,freq'),
    ('NiCl1_aq5',  1, 3, build(OCT, ['Cl']+['O']*5,        R_NI_CL_OH), 'xtb,opt,freq'),
    ('H2O',        0, 1, water((1.,0.,0.), 0.0),                        'xtb,opt,freq'),
    ('Cl_ion',    -1, 1, [('Cl', (0., 0., 0.))],                        'sp'),
]

HEAD = """%pal nprocs 1 end
%maxcore 6000
"""
# 이 서버 Slurm 은 RealMemory=1 로 설정돼 있어 메모리를 관리하지 않는다
# (--mem 을 쓰면 "Memory specification can not be satisfied" 로 제출 자체가 거부됨).
# 따라서 메모리 상한은 오로지 %maxcore 가 결정한다.
# 7잡 x 6 GB = 42 GB, 노드 여유 224 GB -> 다른 사용자에게 지장 없음.
# Slurm 이 OOM 으로 잡아주지 않으므로 maxcore 를 보수적으로 잡는 것이 안전장치다.

# 0단계: XTB2 사전최적화. 기본 단계로 승격 — 초기구조가 손으로 만든 거라
# 물 배향이 나쁘다. 몇 초 만에 DFT opt 사이클을 크게 줄여주고 실패 위험도 낮춘다.
XTB = """! XTB2 Opt ALPB(water)
%geom MaxIter 500 end
""" + HEAD + """* xyzfile {q} {m} {name}_geo.xyz
"""

# 1단계: 기하최적화. TightOpt — 허수진동수 방지용(Freq 앞단계라 조임).
OPT = """! r2SCAN-3c TightOpt CPCM(water)
%geom MaxIter 400 end
%cpcm smd true SMDsolvent "water" end
""" + HEAD + """* xyzfile {q} {m} {name}_geo.xyz
"""

# 2단계: 진동수 -> 열역학(G). 60 C = 333.15 K. 25 C도 같이 뽑아 온도의존성 확인용.
FREQ = """! r2SCAN-3c NumFreq CPCM(water)
%freq Temp 298.15, 333.15 end
%cpcm smd true SMDsolvent "water" end
""" + HEAD + """* xyzfile {q} {m} {name}_opt.xyz
"""

# 참조종 Cl- : 단원자라 opt/freq 자체가 없다. 단일점만.
SP = """! r2SCAN-3c CPCM(water)
%cpcm smd true SMDsolvent "water" end
""" + HEAD + """* xyzfile {q} {m} {name}_geo.xyz
"""

# 3단계: d-d 전이. TD-DFT는 전이금속 d-d에서 오차 큼 -> CASSCF(8,5)/NEVPT2 + AILFT.
# MORead 로 1단계 궤도를 읽어 활성공간 수렴 실패를 줄인다(가장 흔한 사고).
DD = """! def2-TZVP def2-TZVP/C def2/JK RIJK CPCM(water) MORead
%moinp "{name}_opt.gbw"
""" + HEAD + """%casscf
  nel 8
  norb 5
  mult 3,1
  nroots 10,15
  actorbs dorbs
  PTMethod SC_NEVPT2
end
* xyzfile {q} {m} {name}_opt.xyz
"""

# 3단계 재시도: 삼중항만(10 roots). 가시영역 spin-allowed 밴드는 전부 삼중항이라
# 실질 손실 없음. 메모리/root-following 실패 시의 안전판. MORead 도 뺀다.
DD2 = """! def2-TZVP def2-TZVP/C def2/JK RIJK CPCM(water)
""" + HEAD + """%casscf
  nel 8
  norb 5
  mult 3
  nroots 10
  actorbs dorbs
  PTMethod SC_NEVPT2
end
* xyzfile {q} {m} {name}_opt.xyz
"""

here = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(here, 'logs'), exist_ok=True)

def w(fn, txt):
    # encoding/newline 명시 필수: Windows 기본값(cp949 + CRLF)으로 쓰면
    # 서버에서 한글 주석이 깨지고 bash가 \r 때문에 오작동한다.
    open(os.path.join(here, fn), 'w', encoding='utf-8', newline='\n').write(txt)

rows = []
for name, q, mult, at, steps in SPECIES:
    # 시작 구조. bash 가 {name}_geo.xyz 로 복사해서 쓴다 (단계 간 인계 지점).
    w('%s_start.xyz' % name,
      '%d\n%s\n' % (len(at), name) +
      ''.join('%-3s %14.8f %14.8f %14.8f\n' % (e, *xyz) for e, xyz in at))
    f = dict(q=q, m=mult, name=name)
    if 'xtb'  in steps: w(name+'_xtb.inp',  XTB.format(**f))
    if 'opt'  in steps: w(name+'_opt.inp',  OPT.format(**f))
    if 'freq' in steps: w(name+'_freq.inp', FREQ.format(**f))
    if 'sp'   in steps: w(name+'_sp.inp',   SP.format(**f))
    if 'dd'   in steps:
        w(name+'_dd.inp',  DD.format(**f))
        w(name+'_dd2.inp', DD2.format(**f))
    rows.append('%s %s' % (name, steps))
    print('%-11s q=%+d mult=%d  %2d atoms  [%s]' % (name, q, mult, len(at), steps))

w('species.txt', '\n'.join(rows)+'\n')

w('run_array.sh', textwrap.dedent(r"""
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
""").lstrip())

print('\nwrote *_start.xyz, *_xtb/opt/freq/sp/dd/dd2.inp, species.txt, run_array.sh, logs/')
print('제출: sbatch run_array.sh      진행확인: cat STATUS.txt')
