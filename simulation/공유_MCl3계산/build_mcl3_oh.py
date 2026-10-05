"""MCl3(H2O)3 6배위 착물 ORCA 인풋 생성 — 빠진 화학종 후보 검증용 (2026-09-30).

왜 이걸 돌리나 (민이형 지적):
  DFT 는 우리가 넣어준 후보 중에서만 답을 고른다. 기존 사다리는

      M(H2O)6(2+) -> MCl(H2O)5(+) -> MCl2(H2O)4(0) -> MCl3(H2O)(-) -> MCl4(2-)
                        6배위          6배위            ** 4배위 **      4배위

  Cl 2개에서 3개로 갈 때 배위수가 6에서 4로 한 번에 떨어진다. 그런데 그 사이의
  **MCl3(H2O)3 (6배위)** 를 후보에 넣은 적이 없다. 따라서 "Cl 3개째에서 사면체로
  꺾인다" 는 결론이 계산으로 찾은 것인지, 후보를 그렇게 짜서 나온 것인지 구분이 안 된다.
  이 계산이 그 구분을 만든다.

  6배위가 더 안정 -> 전환점이 오른쪽으로 밀린다. 지도를 고쳐야 한다.
  4배위가 더 안정 -> 기존 결론이 후보 선택의 결과가 아님을 보인 것이다.

이성질체 둘 다 돌린다:
  fac (facial)    Cl 3개가 한 면에 모임. 서로 전부 cis.
  mer (meridional) Cl 3개가 자오선에 늘어섬. 둘은 trans, 하나는 cis.
  어느 쪽이 안정한지는 미리 알 수 없고, 둘 중 낮은 쪽을 그 화학종의 에너지로 쓴다.

안전 설계 (기존 run_array.sh 와 동일한 사다리):
  * 단계마다 .done 플래그 -> 재제출해도 끝난 단계는 건너뛴다
  * 수렴 실패 시 SlowConv -> NormalOpt -> KDIIS 순으로 재시도
  * timeout 으로 무한루프 차단, --requeue 로 노드 문제 시 자동 재투입
  * --mem 금지 (이 노드는 RealMemory=1 이라 --mem 을 넣으면 제출 자체가 거부됨)
  * 기존 폴더를 건드리지 않는 별도 출력 폴더 -> 끝난 계산이 망가질 일이 없다

실행:
    python build_mcl3_oh.py          # 인풋 생성 (로컬에서 해도 됨)
    sbatch mcl3_oh/run_mcl3.sh       # 서버에서 제출
    squeue -u $USER                  # 상태 확인
"""
import math
import os
import random
import sys
import textwrap

# 시작 결합거리는 이미 수렴된 기존 계산값에서 가져왔다.
# 좋은 초기추정이 곧 안전장치다 — 수렴 실패와 재시도가 줄어든다.
METALS = {
    'Ni': dict(nel=8, mult=3, r_o=2.11, r_cl=2.37),   # d8 고스핀 S=1
    'Co': dict(nel=7, mult=4, r_o=2.16, r_cl=2.39),   # d7 고스핀 S=3/2
    'Mn': dict(nel=5, mult=6, r_o=2.28, r_cl=2.43),   # d5 고스핀 S=5/2
}
CHARGE = -1                      # MCl3(H2O)3 : M(2+) + 3 Cl(-) = 1-
R_OH, ANG_HOH = 0.97, 104.5
NOISE = 0.04                     # 전 원자 무작위 변위 (A). 대칭 안장점 방지
SEED = 20260930

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'mcl3_oh')
JOBNAME = 'mcl3'

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
    """금속에서 방향 u, 거리 d 에 물 배위. 축 둘레 회전각 무작위 = 대칭 깨기.

    결정론적으로 배치하면 고대칭 구조에 수렴해 허수진동수가 무더기로 나온다
    (Ni_aq6 에서 10개 나왔던 전례). 물마다 다른 각도를 준다.
    """
    O = mul(u, d)
    p = (1., 0., 0.) if abs(u[0]) < 0.9 else (0., 1., 0.)
    e1 = unit(cross(p, u))
    e2 = unit(cross(u, e1))
    phi = RNG.uniform(0, 2 * math.pi)
    v = add(mul(e1, math.cos(phi)), mul(e2, math.sin(phi)))
    h = math.radians(ANG_HOH / 2)
    out = [('O', O)]
    for s in (+1, -1):
        hv = add(mul(u, math.cos(h)), mul(v, s * math.sin(h)))
        out.append(('H', add(O, mul(hv, R_OH))))
    return out


# 팔면체 여섯 방향. fac 은 Cl 이 한 면(+x,+y,+z), mer 은 자오선(+x,-x,+y).
OCT = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
ISOMERS = {
    'fac': [0, 2, 4],        # +x, +y, +z  — 서로 전부 cis
    'mer': [0, 1, 2],        # +x, -x, +y  — 둘은 trans
}


def build(sym, cl_idx, r_o, r_cl):
    at = [(sym, (0., 0., 0.))]
    for i, u in enumerate(OCT):
        at += [('Cl', mul(u, r_cl))] if i in cl_idx else water(u, r_o)
    return [(s, tuple(c + RNG.uniform(-NOISE, NOISE) for c in p)) for s, p in at]


def sanity(name, at, sym):
    """구조가 말이 안 되면 제출 전에 여기서 죽는다. 서버 시간을 버리지 않기 위해."""
    n = len(at)
    worst = min(norm(sub(at[i][1], at[j][1])) for i in range(n) for j in range(i+1, n))
    if worst < 0.85:
        sys.exit('%s: 원자 겹침 %.2f A' % (name, worst))
    m = next(i for i, a in enumerate(at) if a[0] == sym)
    d = sorted(norm(sub(at[m][1], a[1])) for k, a in enumerate(at)
               if k != m and a[0] in ('O', 'Cl'))
    ncoord = len([x for x in d if x < 2.7])
    if ncoord != 6:
        sys.exit('%s: 배위수가 6이 아님 (%d) — 6배위를 검증하려는 계산이다' % (name, ncoord))
    ncl = len([a for a in at if a[0] == 'Cl'])
    no = len([a for a in at if a[0] == 'O'])
    if (ncl, no) != (3, 3):
        sys.exit('%s: Cl %d개 / O %d개 — MCl3(H2O)3 이 아니다' % (name, ncl, no))
    print('  %-14s %2d atoms  최단 %.2f A  배위 6개 (%s)'
          % (name, n, worst, ', '.join('%.2f' % x for x in d[:6])))


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
# 298.15 와 333.15 둘 다 뽑는다: van't Hoff 로 dH/dS 를 분리하려면 두 온도가 필요하다.
FREQ = """! r2SCAN-3c NumFreq CPCM(water)
%freq Temp 298.15, 333.15 end
%cpcm smd true SMDsolvent "water" end
""" + HEAD + """* xyzfile {q} {m} {name}_opt.xyz
"""


def w(fn, txt):
    # encoding/newline 명시 필수: Windows 기본값(cp949 + CRLF)으로 쓰면
    # 서버에서 한글 주석이 깨지고 bash 가 \r 때문에 오작동한다.
    open(os.path.join(OUT, fn), 'w', encoding='utf-8', newline='\n').write(txt)


def main():
    os.makedirs(os.path.join(OUT, 'logs'), exist_ok=True)
    names = []
    print('구조 생성 (물 배향 무작위 = 대칭 깨기):')
    for sym, M in METALS.items():
        for iso, cl_idx in ISOMERS.items():
            name = '%sCl3_aq3_%s' % (sym, iso)
            at = build(sym, cl_idx, M['r_o'], M['r_cl'])
            sanity(name, at, sym)
            w('%s_start.xyz' % name, '%d\n%s\n' % (len(at), name) +
              ''.join('%-3s %14.8f %14.8f %14.8f\n' % (e, *xyz) for e, xyz in at))
            f = dict(q=CHARGE, m=M['mult'], name=name)
            w(name + '_xtb.inp', XTB.format(**f))
            w(name + '_opt.inp', OPT.format(**f))
            w(name + '_freq.inp', FREQ.format(**f))
            names.append(name)

    w('species.txt', ''.join('%s xtb,opt,freq\n' % n for n in names))

    runner = textwrap.dedent("""
        #!/bin/bash
        #SBATCH -J JOBNAME
        #SBATCH -p long
        #SBATCH --array=1-NJOBS%NJOBS
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
          attempt xtb "${N}_xtb" "" "SlowConv" \\
            && { [ -f "${N}_xtb_run.xyz" ] && cp "${N}_xtb_run.xyz" "${N}_geo.xyz"; touch "${N}.xtb.done"; } \\
            || say "xtb 실패 -> 초기구조로 진행"
        fi

        TMAX=21600
        if [ ! -f "${N}.opt.done" ]; then
          attempt opt "${N}_opt" "" "SlowConv" "SlowConv NormalOpt" "VerySlowConv NormalOpt KDIIS"
          if [ -f "${N}_opt_run.xyz" ]; then
            cp "${N}_opt_run.xyz" "${N}_opt.xyz"
            [ -f "${N}_opt_run.gbw" ] && cp "${N}_opt_run.gbw" "${N}_opt.gbw"
            grep -q "ORCA TERMINATED NORMALLY" "${N}_opt.out" && touch "${N}.opt.done" \\
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
    """).lstrip().replace('JOBNAME', JOBNAME).replace('NJOBS', str(len(names)))
    w('run_mcl3.sh', runner)

    print('\nwrote %s/  (%d 종)' % (OUT, len(names)))
    print('제출:  sbatch mcl3_oh/run_mcl3.sh')
    print('확인:  squeue -u $USER   /   tail mcl3_oh/STATUS.txt')

    # 자체검증: 제출 전에 인풋이 온전한지 확인한다. 서버에서 발견하면 하루를 버린다.
    for n in names:
        for step in ('xtb', 'opt', 'freq'):
            p = os.path.join(OUT, '%s_%s.inp' % (n, step))
            assert os.path.exists(p), '인풋 누락: %s' % p
            body = open(p, encoding='utf-8').read()
            assert '{' not in body, '치환 안 된 자리표시자: %s' % p
            assert 'maxcore' in body, 'maxcore 누락: %s' % p
    sh = open(os.path.join(OUT, 'run_mcl3.sh'), encoding='utf-8').read()
    # 주석 안의 '--mem' 글자가 아니라 실제 SBATCH 지시문만 본다
    directives = [l for l in sh.splitlines() if l.strip().startswith('#SBATCH')]
    assert not any('--mem' in l for l in directives), \
        '#SBATCH --mem 이 들어가면 제출이 거부된다'
    assert '\r' not in sh, 'CRLF 가 섞이면 bash 가 오작동한다'
    assert 'array=1-%d' % len(names) in sh, '잡 배열 개수 불일치'
    print('self-check OK: 인풋 %d개, --mem 없음, LF 줄끝, 배열 %d개'
          % (len(names) * 3, len(names)))


if __name__ == '__main__':
    main()
