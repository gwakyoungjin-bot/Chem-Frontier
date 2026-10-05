"""시트르산이 Cl- 를 몇 개까지 강하게 붙잡는가 — 1:1 임계의 계산적 근거.

왜:
  Shafie 2019 (J Mol Liq 288, 111081) Table 2 가 무수 ChCl:CA DES 의 C=O 를
  조성별로 준다:
      1:3  1708.78     1:2  1709.88     1:1  1718.97     2:1  1719.85   3:1  1720.07
  1:2 -> 1:1 구간에서만 +9.09 로 튀고 나머지는 +0.2~1.1 이다. 즉 **계단이 x = 0.5**,
  다시 말해 **ChCl : CA = 1 : 1** 에 있다. 물이 없는 계에서 나온 값이므로
  희석 효과가 아니다. 우리 계(물 0.80)의 계단, 그리고 우리 f_Td 상승 시작점
  (x = 0.45~0.50)과 같은 자리다.

가설:
  시트르산 1분자가 Cl- 1개를 강하게 붙잡는다. 그 이상은 결합이 급격히 약해진다.
  그래서 Cl-/CA = 1 을 넘는 순간 **잉여 Cl-** 가 생기고, 그것이
    (a) 금속의 물을 밀어내 사면체 착물을 만들고
    (b) 시트르산 수소결합 네트워크를 포화시켜 C=O 를 밀어올린다.
  하나의 임계가 두 관측을 만든다.

검증 방법:
  단계별 결합에너지를 본다.
    dE1 = E(CA·Cl-)  - E(CA) - E(Cl-)
    dE2 = E(CA·2Cl-) - E(CA·Cl-) - E(Cl-)
  가설이 참이면 |dE2| << |dE1| 이어야 한다. 두 번째 Cl- 가 첫 번째만큼 잘 붙으면
  1:1 이 특별할 이유가 없어져 가설이 기각된다.

  22~23 원자라 opt+freq 없이 opt+SP 로 1시간 안에 끝난다.
  (엔트로피가 빠지지만 같은 반응형(분자+이온->복합체) 둘의 **비교**라 상당 부분 상쇄된다)

실행:
    python build_clsat.py
    sbatch clsat/run_clsat.sh
"""
import math
import os
import sys
import textwrap

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'clsat')
SRC = os.path.join(os.path.dirname(HERE), 'dft_results', 'geometries',
                   'citrate', 'H3Cit_opt.xyz')      # 이미 DFT 최적화된 시트르산
R_HCL = 2.05                                        # H...Cl- 초기거리 (A)


def read_xyz(p):
    L = [l.split() for l in open(p, encoding='utf-8') if l.strip()]
    n = int(L[0][0])
    return [(r[0], tuple(float(v) for v in r[1:4])) for r in L[2:2 + n]]


def d(a, b): return math.dist(a, b)
def sub(a, b): return tuple(x - y for x, y in zip(a, b))
def add(a, b): return tuple(x + y for x, y in zip(a, b))
def mul(a, s): return tuple(x * s for x in a)
def unit(a):
    n = math.sqrt(sum(x * x for x in a)); return mul(a, 1.0 / n)


def acidic_hydrogens(at):
    """카르복실 O 에 붙은 H. O-H 거리 < 1.1 이고 그 O 가 C 와 1.2~1.6 인 것."""
    out = []
    for i, (e, p) in enumerate(at):
        if e != 'H':
            continue
        for j, (e2, q) in enumerate(at):
            if e2 == 'O' and d(p, q) < 1.15:
                # 그 산소가 카르복실인지 (C 에 붙어 있고, 그 C 에 또 다른 O 가 이중결합)
                for k, (e3, s) in enumerate(at):
                    if e3 == 'C' and d(q, s) < 1.6:
                        n_o = sum(1 for m, (e4, t) in enumerate(at)
                                  if e4 == 'O' and d(s, t) < 1.6)
                        if n_o >= 2:
                            out.append((i, j))
                        break
                break
    return out


def place_cl(at, oh_pairs, k):
    """k 번째 산성 H 의 O-H 축 연장선에 Cl- 를 놓는다."""
    hi, oi = oh_pairs[k]
    u = unit(sub(at[hi][1], at[oi][1]))
    return ('Cl', add(at[hi][1], mul(u, R_HCL)))


HEAD = "%pal nprocs 1 end\n%maxcore 6000\n"
OPT = """! r2SCAN-3c LooseOpt CPCM(water)
%geom MaxIter 150 end
%cpcm smd true SMDsolvent "water" end
""" + HEAD + """* xyzfile {q} 1 {name}_geo.xyz
"""
SP = """! r2SCAN-3c CPCM(water)
%cpcm smd true SMDsolvent "water" end
""" + HEAD + """* xyzfile {q} 1 {name}_sp.xyz
"""


def w(fn, txt):
    open(os.path.join(OUT, fn), 'w', encoding='utf-8', newline='\n').write(txt)


def main():
    if not os.path.exists(SRC):
        sys.exit('H3Cit_opt.xyz 없음: %s' % SRC)
    base = read_xyz(SRC)
    pairs = acidic_hydrogens(base)
    print('시트르산 %d원자, 산성 H %d개 검출' % (len(base), len(pairs)))
    if len(pairs) < 2:
        sys.exit('산성 H 를 2개 이상 못 찾았다')

    os.makedirs(os.path.join(OUT, 'logs'), exist_ok=True)
    species = {
        'H3Cit': (list(base), 0),
        'H3Cit_Cl1': (list(base) + [place_cl(base, pairs, 0)], -1),
        'H3Cit_Cl2': (list(base) + [place_cl(base, pairs, 0),
                                    place_cl(base, pairs, 1)], -2),
        'Cl_ion': ([('Cl', (0.0, 0.0, 0.0))], -1),
    }
    names = []
    for name, (at, q) in species.items():
        n = len(at)
        if n > 1:
            worst = min(d(at[i][1], at[j][1]) for i in range(n) for j in range(i + 1, n))
            if worst < 0.85:
                sys.exit('%s: 원자 겹침 %.2f' % (name, worst))
            print('  %-11s %2d원자  전하 %+d  최단 %.2f A' % (name, n, q, worst))
        else:
            print('  %-11s %2d원자  전하 %+d' % (name, n, q))
        w('%s_start.xyz' % name, '%d\n%s\n' % (n, name) +
          ''.join('%-3s %12.6f %12.6f %12.6f\n' % (e, *p) for e, p in at))
        f = dict(q=q, name=name)
        w('%s_opt.inp' % name, OPT.format(**f))
        w('%s_sp.inp' % name, SP.format(**f))
        names.append(name)
    w('species.txt', ''.join(n + '\n' for n in names))

    runner = textwrap.dedent("""
        #!/bin/bash
        #SBATCH -J clsat
        #SBATCH -p long
        #SBATCH --array=1-NJ%NJ
        #SBATCH -c 1
        #SBATCH -t 0-03:00:00
        # --mem 금지 (RealMemory=1)
        #SBATCH -o logs/slurm-%A_%a.out
        #SBATCH --requeue

        cd "$SLURM_SUBMIT_DIR" || exit 1
        [ -f ~/software/activate_orca_mpi.sh ] && source ~/software/activate_orca_mpi.sh
        N=$(sed -n "${SLURM_ARRAY_TASK_ID}p" species.txt)
        [ -z "$N" ] && exit 0
        say(){ echo "[$(date '+%H:%M:%S')] $N $*" >> STATUS.txt; }
        command -v orca >/dev/null || { say FATAL; exit 1; }
        run(){ local s=$1 t=$2
          [ -f "${N}.$s.done" ] && return 0
          local t0=$SECONDS
          timeout -k 30 "$t" orca "${N}_${s}.inp" > "${N}_${s}.out" 2>&1
          grep -q "ORCA TERMINATED NORMALLY" "${N}_${s}.out" \\
            && { say "$s OK $((SECONDS-t0))s"; touch "${N}.$s.done"; rm -f "${N}_${s}"*.tmp*; return 0; }
          say "$s 실패 $((SECONDS-t0))s"; rm -f "${N}_${s}"*.tmp*; return 1; }

        cp "${N}_start.xyz" "${N}_geo.xyz"
        cp "${N}_geo.xyz" "${N}_sp.xyz"
        run sp 900                                  # 먼저 답 하나 확보
        if run opt 2400; then
          [ -f "${N}_opt.xyz" ] && cp "${N}_opt.xyz" "${N}_sp.xyz" && rm -f "${N}.sp.done" && run sp 900
        fi
        say DONE
    """).lstrip().replace('NJ', str(len(names)))
    w('run_clsat.sh', runner)
    print('\nwrote %s/  (%d종)' % (OUT, len(names)))
    sh = open(os.path.join(OUT, 'run_clsat.sh'), encoding='utf-8').read()
    assert not any('--mem' in l for l in sh.splitlines() if l.strip().startswith('#SBATCH'))
    assert '\r' not in sh
    print('self-check OK')


if __name__ == '__main__':
    main()
