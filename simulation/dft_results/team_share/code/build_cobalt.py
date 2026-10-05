"""Co(II) 아쿠아/클로로 착물 ORCA 인풋 생성 (2026-09-27).

왜 코발트인가:
  Nature Commun 2021 (10.1038/s41467-021-26814-7) — 염화물 농도로 Co 는 음이온성
  CoCl4(2-) 가 되고 Ni 는 양이온 [NiCl(H2O)5]+ 로 남는다. 이 '갈라짐'이
  Co/Ni 분리(음이온교환·용매추출)의 물리적 근거다.
  => 두 금속이 서로 반대 방향을 원하므로 운전창에 '양쪽 경계'가 생긴다:
       Co 전환선은 넘되 Ni 전환선은 넘지 않는 구간.
  Ni 만으로는 단방향 '피하라'밖에 못 만든다. Co 를 계산해야 창이 닫힌다.

전자구조:
  Co(II) d7 고스핀 -> S = 3/2 -> mult 4.
  CAS(7,5) 의 사중항은 4F(7) + 4P(3) = 10개 -> nroots 10 이 d7 사중항 전체를 덮는다.
  (d8 의 삼중항 10개와 입자-홀 대응)

** 대칭 안장점 교훈 반영 **
  build_inputs.py 는 물 H 를 결정론적으로 배치해 고대칭 구조를 만들었고,
  Ni_aq6/NiCl2_aq4 가 안장점에 수렴해 허수진동수가 10개/8개 나왔다.
  여기서는 처음부터 물 배향을 무작위로 틀고 미세 변위를 준다.

실행:
    python build_cobalt.py
    sbatch cobalt/run_cobalt.sh
"""
import math, os, random, sys, textwrap

# 금속별 고스핀 결합거리(실험 결정구조 전형값)와 스핀다중도.
#   Co(II) d7 고스핀 S=3/2 -> mult 4
#   Mn(II) d5 고스핀 S=5/2 -> mult 6 (반충만각이라 SCF 가 오히려 잘 수렴한다)
# Mn 은 d-d 가 전부 스핀금지라 분광값이 무의미하므로 dd 단계를 넣지 않는다.
METALS = {
    'Co': dict(sym='Co', nel=7, mult=4, r_o=2.08, r_cl_oh=2.45, r_cl_td=2.28,
               out='cobalt',    dd=True),
    'Mn': dict(sym='Mn', nel=5, mult=6, r_o=2.18, r_cl_oh=2.50, r_cl_td=2.37,
               out='manganese', dd=False),
}
METAL = 'Co'
if '--metal' in sys.argv:
    METAL = sys.argv[sys.argv.index('--metal') + 1]
    if METAL not in METALS:
        sys.exit('금속은 %s 중 하나' % ', '.join(METALS))
M = METALS[METAL]
# 잡 이름은 둘 다 co_dft: 이미 돌고 있는 watch_cobalt.py 가 JOBNAME='co_dft' 를
# 기다리므로, 같은 이름이면 Mn 까지 끝난 뒤 후처리가 돈다 (감독 수정 불필요).
JOBNAME = 'co_dft'
R_CO_O, R_CO_CL_OH, R_CO_CL_TD = M['r_o'], M['r_cl_oh'], M['r_cl_td']
R_OH, ANG_HOH = 0.97, 104.5
NOISE = 0.04                 # 전 원자 무작위 변위 (Angstrom)
SEED = 20260927

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, M['out'])
# 시드 교체 가능: 허수진동수가 남으면 watch_cobalt.py 가 다른 시드로 다시 부른다.
if '--seed' in sys.argv:
    SEED = int(sys.argv[sys.argv.index('--seed') + 1])
RNG = random.Random(SEED)


def sub(a, b): return tuple(x - y for x, y in zip(a, b))
def add(a, b): return tuple(x + y for x, y in zip(a, b))
def mul(a, s): return tuple(x * s for x in a)
def dot(a, b): return sum(x * y for x, y in zip(a, b))
def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
def norm(a): return math.sqrt(dot(a, a))
def unit(a): return mul(a, 1.0 / norm(a))


def water(u, d):
    """금속에서 방향 u, 거리 d 에 물 배위. 축 둘레 회전각을 무작위로 준다(대칭 깨기)."""
    O = mul(u, d)
    p = (1., 0., 0.) if abs(u[0]) < 0.9 else (0., 1., 0.)
    e1 = unit(cross(p, u))
    e2 = unit(cross(u, e1))
    phi = RNG.uniform(0, 2 * math.pi)                 # <- 핵심: 물마다 다른 배향
    v = add(mul(e1, math.cos(phi)), mul(e2, math.sin(phi)))
    h = math.radians(ANG_HOH / 2)
    out = [('O', O)]
    for s in (+1, -1):
        hv = add(mul(u, math.cos(h)), mul(v, s * math.sin(h)))
        out.append(('H', add(O, mul(hv, R_OH))))
    return out


OCT = [(1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1)]
s3 = 1/math.sqrt(3)
TET = [(s3,s3,s3), (s3,-s3,-s3), (-s3,s3,-s3), (-s3,-s3,s3)]


def build(dirs, ligs, d_cl):
    at = [(M['sym'], (0., 0., 0.))]
    for u, L in zip(dirs, ligs):
        at += [('Cl', mul(u, d_cl))] if L == 'Cl' else water(u, R_CO_O)
    # 잔여 대칭 제거
    return [(s, tuple(c + RNG.uniform(-NOISE, NOISE) for c in p)) for s, p in at]


# (이름, 전하, 기하, 리간드, Cl거리, 단계).  전부 Co(II) d7 고스핀 -> mult 4
_S, _E = M['sym'], (',dd' if M['dd'] else '')
SPECIES = [
    (_S+'Cl4',     -2, TET, ['Cl']*4,            R_CO_CL_TD, 'xtb,opt,freq'+_E),
    (_S+'_aq6',     2, OCT, ['O']*6,             R_CO_CL_OH, 'xtb,opt,freq'+_E),
    (_S+'Cl2_aq4',  0, OCT, ['Cl','Cl']+['O']*4, R_CO_CL_OH, 'xtb,opt,freq'),
    (_S+'Cl3_aq1', -1, TET, ['Cl']*3+['O'],      R_CO_CL_TD, 'xtb,opt,freq'),
    (_S+'Cl1_aq5',  1, OCT, ['Cl']+['O']*5,      R_CO_CL_OH, 'xtb,opt,freq'),
]
MULT = M['mult']

HEAD = "%pal nprocs 1 end\n%maxcore 6000\n"
XTB = """! XTB2 Opt ALPB(water)
%geom MaxIter 600 end
""" + HEAD + """* xyzfile {q} {m} {name}_geo.xyz
"""
OPT = """! r2SCAN-3c TightOpt CPCM(water)
%geom MaxIter 500 end
%cpcm smd true SMDsolvent "water" end
""" + HEAD + """* xyzfile {q} {m} {name}_geo.xyz
"""
FREQ = """! r2SCAN-3c NumFreq CPCM(water)
%freq Temp 298.15, 333.15 end
%cpcm smd true SMDsolvent "water" end
""" + HEAD + """* xyzfile {q} {m} {name}_opt.xyz
"""
# d7 사중항 10개. Ni 때 Ni_aq6 의 mult 3,1 조합이 수렴 실패했으므로 처음부터 단일 다중도로 간다.
DD = """! def2-TZVP def2-TZVP/C def2/JK RIJK CPCM(water) MORead
%moinp "{name}_opt.gbw"
""" + HEAD + """%casscf
  nel {nel}
  norb 5
  mult {m}
  nroots 10
  actorbs dorbs
  PTMethod SC_NEVPT2
end
* xyzfile {q} {m} {name}_opt.xyz
"""
DD2 = """! def2-TZVP def2-TZVP/C def2/JK RIJK CPCM(water)
""" + HEAD + """%casscf
  nel {nel}
  norb 5
  mult {m}
  nroots 10
  actorbs dorbs
  PTMethod SC_NEVPT2
end
* xyzfile {q} {m} {name}_opt.xyz
"""


def w(fn, txt):
    open(os.path.join(OUT, fn), 'w', encoding='utf-8', newline='\n').write(txt)


def sanity(name, at):
    n = len(at)
    worst = min(norm(sub(at[i][1], at[j][1])) for i in range(n) for j in range(i+1, n))
    if worst < 0.85:
        sys.exit('%s: 원자 겹침 %.2f A' % (name, worst))
    co = next(i for i, a in enumerate(at) if a[0] == M['sym'])
    d = sorted(norm(sub(at[co][1], a[1])) for k, a in enumerate(at)
               if k != co and a[0] in ('O', 'Cl'))
    ncoord = len([x for x in d if x < 2.7])
    print('  %-11s %2d atoms  최단 %.2f A  %s 배위 %d개 (%s)'
          % (name, n, worst, M['sym'], ncoord, ', '.join('%.2f' % x for x in d[:6])))
    if ncoord not in (4, 6):
        sys.exit('  -> 배위수가 4도 6도 아님')


def main():
    os.makedirs(os.path.join(OUT, 'logs'), exist_ok=True)
    rows = []
    print('구조 생성 (물 배향 무작위 = 대칭 깨기):')
    for name, q, dirs, ligs, dcl, steps in SPECIES:
        at = build(dirs, ligs, dcl)
        sanity(name, at)
        w('%s_start.xyz' % name, '%d\n%s\n' % (len(at), name) +
          ''.join('%-3s %14.8f %14.8f %14.8f\n' % (e, *xyz) for e, xyz in at))
        f = dict(q=q, m=MULT, name=name, nel=M['nel'])
        w(name+'_xtb.inp', XTB.format(**f))
        w(name+'_opt.inp', OPT.format(**f))
        w(name+'_freq.inp', FREQ.format(**f))
        if 'dd' in steps:
            w(name+'_dd.inp', DD.format(**f))
            w(name+'_dd2.inp', DD2.format(**f))
        rows.append('%s %s' % (name, steps))
    w('species.txt', '\n'.join(rows)+'\n')

    runner = textwrap.dedent(r"""
        #!/bin/bash
        #SBATCH -J JOBNAME_PLACEHOLDER
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
    """).lstrip().replace('JOBNAME_PLACEHOLDER', JOBNAME)
    w('run_%s.sh' % M['out'], runner)
    print('\nwrote %s/' % OUT)
    print('제출: sbatch %s/run_%s.sh' % (M['out'], M['out']))


if __name__ == '__main__':
    main()
