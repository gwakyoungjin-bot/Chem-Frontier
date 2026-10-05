"""C=O 를 내리는 게 무엇인지 가리는 모형화합물 계산 — DFT 만으로 인과를 판정한다.

무엇을 가르려 하나:
  실측 FT-IR 은 CA 가 많을수록 C=O 가 **낮다**(1713.9)
                ChCl 이 많을수록 **높다**(1729.2).  차이 15 cm-1.
  이 이동을 만드는 파트너가 누구냐에 따라 인과 판정이 갈린다:

    Cl- 가 C=O 를 내린다면 -> ChCl 쪽이 낮아야 하는데 실측은 반대 -> 염화물 가설 반증
                              동시에 'FT-IR 은 염화물을 못 본다' 가 확정된다
    이량체가 내린다면      -> CA 쪽이 낮다 -> 실측과 일치
                              FT-IR 은 CA 자기회합을 보는 것이고,
                              화학종 전이와는 인과가 아니라 공통 조성의존이다

  즉 이 계산 하나로 "FT-IR 이 염화물의 지표인가" 에 답이 난다. 실험이 필요 없다.

모형화합물:
  시트르산(21원자) 대신 아세트산(8원자)을 COOH 대리로 쓴다. C=O 수소결합 이동은
  국소적이라 곁사슬에 둔감하다. 8~16원자라 opt+freq 가 몇 분이면 끝난다.
  방향과 크기만 보면 되므로(정확도 불필요) 이 근사가 정당하다.

계산 종:
  AcOH            자유 단량체 — 기준
  AcOH_Cl         Cl- 가 O-H 에 수소결합       (ChCl 많은 환경)
  AcOH_dimer      COOH...COOH 고리형 이량체    (CA 많은 환경)
  AcOH_w          물 1개가 O-H 에 수소결합     (항상 존재하는 배경)

판정은 Δν(C=O) = ν(착물) - ν(자유 AcOH) 의 부호와 크기로 한다.

실행:
    python build_cofreq.py
    sbatch cofreq/run_cofreq.sh
"""
import math
import os
import sys
import textwrap

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'cofreq')

# 아세트산 (평면 syn 배좌). C=O 1.21, C-OH 1.34, O-H 0.97 A.
ACOH = [
    ('C', (0.000, 0.000, 0.000)),      # 메틸 C
    ('H', (-0.510, 0.940, 0.190)),
    ('H', (-0.530, -0.900, 0.280)),
    ('H', (0.180, -0.070, -1.070)),
    ('C', (1.320, 0.030, 0.720)),      # 카르복실 C
    ('O', (1.430, 0.060, 1.930)),      # C=O
    ('O', (2.390, 0.020, -0.110)),     # C-OH
    ('H', (3.190, 0.040, 0.430)),      # 산성 H
]
H_IDX, OH_O, CO_O = 7, 6, 5            # 산성 H, 하이드록실 O, 카르보닐 O


def sub(a, b): return tuple(x - y for x, y in zip(a, b))
def add(a, b): return tuple(x + y for x, y in zip(a, b))
def mul(a, s): return tuple(x * s for x in a)
def norm(a): return math.sqrt(sum(x * x for x in a))
def unit(a): return mul(a, 1.0 / norm(a))


def along_OH(d):
    """산성 H 의 O-H 축을 연장한 위치. 수소결합 받개를 여기에 놓는다."""
    o, h = ACOH[OH_O][1], ACOH[H_IDX][1]
    u = unit(sub(h, o))
    return add(h, mul(u, d))


def species():
    out = {}
    out['AcOH'] = (list(ACOH), 0, 1)

    # Cl- 가 산성 H 에 수소결합 (H...Cl 약 2.0 A)
    out['AcOH_Cl'] = (list(ACOH) + [('Cl', along_OH(2.05))], -1, 1)

    # 물 1개가 산성 H 를 받음 (H...O 약 1.7 A), 물은 대충 배향
    ow = along_OH(1.75)
    out['AcOH_w'] = (list(ACOH) + [('O', ow),
                                   ('H', add(ow, (0.55, 0.78, 0.10))),
                                   ('H', add(ow, (0.55, -0.60, 0.55)))], 0, 1)

    # 고리형 이량체: 두 번째 AcOH 를 카르복실 평면에서 180도 회전 + 평행이동
    # C=O ... H-O 가 서로 맞물리도록 O...O 약 2.7 A 간격으로 배치한다
    # 반전중심을 눈대중으로 잡으면 원자가 겹친다. H...O=C 를 1.75 A 로 두는
    # 위치를 역산해서 정한다: 받개 O 를 O-H 축 연장선에 놓고, 그 중점이 반전중심.
    acc = along_OH(1.75)
    co = ACOH[CO_O][1]
    cen = tuple((a + b) / 2 for a, b in zip(acc, co))
    second = [(e, tuple(2 * c - x for c, x in zip(cen, p))) for e, p in ACOH]
    out['AcOH_dimer'] = (list(ACOH) + second, 0, 1)
    return out


HEAD = "%pal nprocs 1 end\n%maxcore 6000\n"
OPT = """! r2SCAN-3c TightOpt CPCM(water)
%geom MaxIter 300 end
%cpcm smd true SMDsolvent "water" end
""" + HEAD + """* xyzfile {q} {m} {name}_geo.xyz
"""
FREQ = """! r2SCAN-3c NumFreq CPCM(water)
%cpcm smd true SMDsolvent "water" end
""" + HEAD + """* xyzfile {q} {m} {name}_opt.xyz
"""


def w(fn, txt):
    open(os.path.join(OUT, fn), 'w', encoding='utf-8', newline='\n').write(txt)


def check(name, at):
    n = len(at)
    worst = min(norm(sub(at[i][1], at[j][1])) for i in range(n) for j in range(i + 1, n))
    if worst < 0.85:
        sys.exit('%s: 원자 겹침 %.2f A' % (name, worst))
    return worst


def main():
    os.makedirs(os.path.join(OUT, 'logs'), exist_ok=True)
    names = []
    print('모형화합물 생성 (아세트산을 COOH 대리로):')
    for name, (at, q, m) in species().items():
        worst = check(name, at)
        print('  %-12s %2d원자  전하 %+d  최단 %.2f A' % (name, len(at), q, worst))
        w('%s_start.xyz' % name, '%d\n%s\n' % (len(at), name) +
          ''.join('%-3s %12.6f %12.6f %12.6f\n' % (e, *p) for e, p in at))
        f = dict(q=q, m=m, name=name)
        w('%s_opt.inp' % name, OPT.format(**f))
        w('%s_freq.inp' % name, FREQ.format(**f))
        names.append(name)
    w('species.txt', ''.join(n + '\n' for n in names))

    runner = textwrap.dedent("""
        #!/bin/bash
        #SBATCH -J cofreq
        #SBATCH -p long
        #SBATCH --array=1-NJ%NJ
        #SBATCH -c 1
        #SBATCH -t 0-06:00:00
        # --mem 금지 (RealMemory=1)
        #SBATCH -o logs/slurm-%A_%a.out
        #SBATCH --requeue

        cd "$SLURM_SUBMIT_DIR" || exit 1
        [ -f ~/software/activate_orca_mpi.sh ] && source ~/software/activate_orca_mpi.sh
        N=$(sed -n "${SLURM_ARRAY_TASK_ID}p" species.txt)
        [ -z "$N" ] && exit 0
        say(){ echo "[$(date '+%H:%M:%S')] $N $*" >> STATUS.txt; }
        command -v orca >/dev/null || { say "FATAL orca 없음"; exit 1; }

        run(){ local s=$1 t=$2
          [ -f "${N}.$s.done" ] && return 0
          local t0=$SECONDS
          timeout -k 30 "$t" orca "${N}_${s}.inp" > "${N}_${s}.out" 2>&1
          if grep -q "ORCA TERMINATED NORMALLY" "${N}_${s}.out"; then
            say "$s OK $((SECONDS-t0))s"; touch "${N}.$s.done"; rm -f "${N}_${s}"*.tmp*; return 0
          fi
          say "$s 실패 $((SECONDS-t0))s"; rm -f "${N}_${s}"*.tmp*; return 1
        }

        cp "${N}_start.xyz" "${N}_geo.xyz"
        run opt 5400 && [ -f "${N}_opt.xyz" ] || cp "${N}_geo.xyz" "${N}_opt.xyz"
        run freq 10800
        say DONE
    """).lstrip().replace('NJ', str(len(names)))
    w('run_cofreq.sh', runner)

    print('\nwrote %s/  (%d종)' % (OUT, len(names)))
    print('제출: sbatch cofreq/run_cofreq.sh')
    sh = open(os.path.join(OUT, 'run_cofreq.sh'), encoding='utf-8').read()
    assert not any('--mem' in l for l in sh.splitlines()
                   if l.strip().startswith('#SBATCH'))
    assert '\r' not in sh
    print('self-check OK: --mem 없음, LF 줄끝, 인풋 %d개' % (len(names) * 2))


if __name__ == '__main__':
    main()
