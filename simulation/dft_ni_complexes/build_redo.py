"""대칭 안장점에 갇힌 종을 대칭을 깨서 재최적화.

왜 필요한가 (2026-09-26 진단):
  build_inputs.py 의 water() 는 물 H 를 결정론적 규칙으로 배치한다 -> 초기구조가 고대칭.
  최적화는 대칭을 스스로 깨지 못하므로(대칭 좌표 방향 기울기가 0) 안장점에 수렴한다.
  실제로 허수진동수 개수가 배위수 물 개수를 따라갔다:
      Ni_aq6(물6) 허수10개 / NiCl2_aq4(물4) 8개 / NiCl3_aq1(물1),NiCl4(물0) 0개
      시트르산이 대칭을 깨는 NiCit_aq4(물4) 는 0개  <- 진단의 결정적 근거

  ORCA 열역학은 허수모드를 분배함수에서 제외한다. 즉 기준물질 Ni_aq6 가 모드 10개를
  잃어 엔트로피가 과소평가되고, 그 오차가 모든 dS 에 그대로 실린다.
  (재계산 전 값: Td 종 dS = +430~452 J/mol/K, 보정 shift +71 kJ/mol — 둘 다 비현실적)

무엇을 하나:
  기존 최적화 기하에서 출발해 배위수 물을 각각 Ni-O 축 둘레로 무작위 회전시키고,
  전체에 작은 무작위 변위를 준 뒤 다시 xtb -> opt -> freq.
  Ni-O 거리는 그대로 두므로 좋은 기하를 버리지 않고 대칭만 깬다.

실행:
    cd ~/dft_ni_complexes
    source ~/venvs/cosmors/bin/activate
    python build_redo.py
    sbatch redo/run_redo.sh
"""
import math, os, random, sys, textwrap

REDO = ['Ni_aq6', 'NiCl2_aq4']          # 허수진동수가 큰 종만
CHARGE = {'Ni_aq6': 2, 'NiCl2_aq4': 0}
NOISE = 0.05                             # 전 원자 무작위 변위 (Angstrom)
SEED = 20260926

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'redo')


def sub(a, b): return tuple(x - y for x, y in zip(a, b))
def add(a, b): return tuple(x + y for x, y in zip(a, b))
def mul(a, s): return tuple(x * s for x in a)
def dot(a, b): return sum(x * y for x, y in zip(a, b))
def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
def norm(a): return math.sqrt(dot(a, a))
def unit(a): return mul(a, 1.0 / norm(a))


def rotate(p, origin, axis, ang):
    """Rodrigues 회전: origin 을 지나고 방향이 axis 인 축 둘레로 p 를 ang 만큼."""
    v = sub(p, origin)
    k = unit(axis)
    c, s = math.cos(ang), math.sin(ang)
    r = add(add(mul(v, c), mul(cross(k, v), s)), mul(k, dot(k, v) * (1 - c)))
    return add(origin, r)


def read_xyz(p):
    lines = [l for l in open(p).read().splitlines()[2:] if l.strip()]
    return [(x.split()[0], tuple(map(float, x.split()[1:4]))) for x in lines]


def break_symmetry(atoms, rng):
    """배위수 물을 Ni-O 축 둘레로 무작위 회전 + 전체 미세 변위."""
    ni = next(i for i, a in enumerate(atoms) if a[0] == 'Ni')
    P = [list(a[1]) for a in atoms]
    sym = [a[0] for a in atoms]

    # 배위 산소: Ni 에서 2.5 A 이내의 O
    ox = [i for i, a in enumerate(atoms)
          if a[0] == 'O' and norm(sub(a[1], atoms[ni][1])) < 2.5]
    nrot = 0
    for io in ox:
        # 그 산소에 붙은 수소 (1.3 A 이내)
        hs = [i for i, a in enumerate(atoms)
              if a[0] == 'H' and norm(sub(a[1], atoms[io][1])) < 1.3]
        if len(hs) != 2:
            continue                                  # 물이 아니면 건드리지 않음
        axis = sub(atoms[io][1], atoms[ni][1])         # Ni->O 축
        ang = rng.uniform(math.radians(40), math.radians(140)) * rng.choice([1, -1])
        for ih in hs:
            P[ih] = list(rotate(tuple(P[ih]), atoms[io][1], axis, ang))
        nrot += 1

    for i in range(len(P)):                            # 잔여 대칭까지 확실히 제거
        for k in range(3):
            P[i][k] += rng.uniform(-NOISE, NOISE)
    return [(sym[i], tuple(P[i])) for i in range(len(P))], nrot, len(ox)


HEAD = "%pal nprocs 1 end\n%maxcore 6000\n"
XTB = """! XTB2 Opt ALPB(water)
%geom MaxIter 800 end
""" + HEAD + """* xyzfile {q} 3 {name}_geo.xyz
"""
OPT = """! r2SCAN-3c TightOpt CPCM(water)
%geom MaxIter 600 end
%cpcm smd true SMDsolvent "water" end
""" + HEAD + """* xyzfile {q} 3 {name}_geo.xyz
"""
FREQ = """! r2SCAN-3c NumFreq CPCM(water)
%freq Temp 298.15, 333.15 end
%cpcm smd true SMDsolvent "water" end
""" + HEAD + """* xyzfile {q} 3 {name}_opt.xyz
"""


def w(fn, txt):
    open(os.path.join(OUT, fn), 'w', encoding='utf-8', newline='\n').write(txt)


def main():
    os.makedirs(os.path.join(OUT, 'logs'), exist_ok=True)
    # 시드를 바꿔가며 여러 번 시도할 수 있게 한다 (overnight.py 가 라운드마다 다른 시드로 호출).
    seed = SEED
    if '--seed' in sys.argv:
        seed = int(sys.argv[sys.argv.index('--seed') + 1])
    print('seed = %d' % seed)
    rng = random.Random(seed)
    rows = []
    for name in REDO:
        src = os.path.join(HERE, '%s_opt.xyz' % name)
        if not os.path.exists(src):
            sys.exit('없음: %s' % src)
        at = read_xyz(src)
        new, nrot, nox = break_symmetry(at, rng)

        # 건전성: 원자 겹침 없고 Ni-O 거리가 보존됐는지
        worst = min(norm(sub(new[i][1], new[j][1]))
                    for i in range(len(new)) for j in range(i + 1, len(new)))
        ni = next(i for i, a in enumerate(new) if a[0] == 'Ni')
        d = sorted(norm(sub(new[ni][1], a[1])) for k, a in enumerate(new)
                   if k != ni and a[0] in ('O', 'Cl'))[:6]
        if worst < 0.85:
            sys.exit('%s: 원자 겹침 %.2f A' % (name, worst))
        print('  %-11s 물 %d/%d개 회전  최단거리 %.2f A  Ni 배위 %s'
              % (name, nrot, nox, worst, ', '.join('%.2f' % x for x in d)))

        w('%s_start.xyz' % name,
          '%d\n%s (대칭깨짐)\n' % (len(new), name) +
          ''.join('%-3s %14.8f %14.8f %14.8f\n' % (e, *xyz) for e, xyz in new))
        f = dict(q=CHARGE[name], name=name)
        w(name + '_xtb.inp', XTB.format(**f))
        w(name + '_opt.inp', OPT.format(**f))
        w(name + '_freq.inp', FREQ.format(**f))
        rows.append(name)
    w('species.txt', '\n'.join(rows) + '\n')

    w('run_redo.sh', textwrap.dedent(r"""
        #!/bin/bash
        #SBATCH -J ni_redo
        #SBATCH -p long
        #SBATCH --array=1-2%2
        #SBATCH -c 1
        #SBATCH -t 2-00:00:00
        # --mem 금지 (RealMemory=1). 메모리는 %maxcore 로만 제한.
        #SBATCH -o logs/slurm-%A_%a.out
        #SBATCH --requeue

        cd "$SLURM_SUBMIT_DIR" || exit 1
        [ -f ~/software/activate_orca_mpi.sh ] && source ~/software/activate_orca_mpi.sh
        N=$(sed -n "${SLURM_ARRAY_TASK_ID}p" species.txt)
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
            || say "xtb 실패 -> 대칭깨진 초기구조 그대로 진행"
        fi

        TMAX=21600
        if [ ! -f "${N}.opt.done" ]; then
          attempt opt "${N}_opt" "" "SlowConv" "SlowConv NormalOpt"
          if [ -f "${N}_opt_run.xyz" ]; then
            cp "${N}_opt_run.xyz" "${N}_opt.xyz"
            grep -q "ORCA TERMINATED NORMALLY" "${N}_opt.out" && touch "${N}.opt.done" \
              || say "opt 미수렴 -> 마지막 기하 사용"
          else
            cp "${N}_geo.xyz" "${N}_opt.xyz"; say "opt 결과 없음 -> 초기구조 사용"
          fi
        fi

        TMAX=43200
        if [ ! -f "${N}.freq.done" ]; then
          attempt freq "${N}_freq" "" "SlowConv" && touch "${N}.freq.done" \
            || say "freq 실패"
        fi
        say "DONE"
    """).lstrip())
    print('\nwrote %s/\n제출: sbatch redo/run_redo.sh' % OUT)


if __name__ == '__main__':
    main()
