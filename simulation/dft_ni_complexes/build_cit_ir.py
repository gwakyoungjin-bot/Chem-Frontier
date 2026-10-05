"""시트르산 자체의 nu(C=O) 를 염화물 개수별로 — FT-IR 조성 시리즈의 in-silico 판.

왜 이게 결정적인가:
  실험은 조성을 바꾸며 시트르산의 C=O 를 본다. 조성을 바꾸는 것은 곧
  시트르산 1분자가 만나는 Cl- 개수를 바꾸는 것이다. 그러면 계산에서
  **Cl- 를 0, 1, 2 개 붙여가며 nu(C=O) 를 구하면 실험 시리즈를 그대로 재현**하는 셈이다.
  아세트산 대리가 아니라 **실제 시트르산**으로 하므로 실측과 직접 비교된다.

속도:
  전 원자 NumFreq 는 22원자에서 1시간이 넘는다. 우리는 C=O 신축만 필요하므로
  **부분 헤시안**으로 카르복실기 원자만 움직인다. 헤시안 차원이 수십분의 1 로 줄어
  몇 분이면 끝난다. 절대 파수는 조화근사라 실측보다 높게 나오지만, 우리가 쓰는 것은
  **Cl- 개수에 따른 이동량**이므로 문제되지 않는다.

기하:
  clsat 잡이 이미 최적화해 둔 H3Cit / H3Cit_Cl1 / H3Cit_Cl2 를 그대로 쓴다.
  따라서 opt 를 다시 돌리지 않는다.

실행 (서버에서):
    python build_cit_ir.py && sbatch cit_ir/run_cit_ir.sh
"""
import os
import sys
import textwrap

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'cit_ir')
CLSAT = os.path.join(HERE, 'clsat')
SPECIES = [('H3Cit', 0), ('H3Cit_Cl1', -1), ('H3Cit_Cl2', -2)]

HEAD = "%pal nprocs 1 end\n%maxcore 6000\n"
# 부분 헤시안: 카르복실 원자(C, O, O, H)만 움직인다. 인덱스는 생성 시 계산.
FREQ = """! r2SCAN-3c NumFreq CPCM(water)
%cpcm smd true SMDsolvent "water" end
%freq
  PARTIAL_Hess true
  Hess_Atoms {{{atoms}}} end
end
""" + HEAD + """* xyzfile {q} 1 {name}.xyz
"""


def read_xyz(p):
    L = [l.split() for l in open(p, encoding='utf-8') if l.strip()]
    n = int(L[0][0])
    return [(r[0], tuple(float(v) for v in r[1:4])) for r in L[2:2 + n]]


def dist(a, b):
    return sum((x - y) ** 2 for x, y in zip(a, b)) ** 0.5


def carboxyl_atoms(at):
    """카르복실기(C + 산소 2개 + 그 H)의 0-기반 인덱스."""
    idx = set()
    for i, (e, p) in enumerate(at):
        if e != 'C':
            continue
        os_ = [j for j, (e2, q) in enumerate(at) if e2 == 'O' and dist(p, q) < 1.60]
        if len(os_) < 2:
            continue                      # 카르복실 탄소가 아니다
        idx.add(i); idx.update(os_)
        for j in os_:
            for k, (e3, r) in enumerate(at):
                if e3 == 'H' and dist(at[j][1], r) < 1.15:
                    idx.add(k)
    return sorted(idx)


def main():
    os.makedirs(os.path.join(OUT, 'logs'), exist_ok=True)
    names = []
    print('부분 헤시안 대상 원자 (0-기반):')
    for name, q in SPECIES:
        src = None
        for cand in ('%s_opt.xyz' % name, '%s_sp.xyz' % name, '%s_start.xyz' % name):
            p = os.path.join(CLSAT, cand)
            if os.path.exists(p):
                src = p; break
        if src is None:
            print('  %-11s 기하 없음 — 건너뜀' % name)
            continue
        at = read_xyz(src)
        idx = carboxyl_atoms(at)
        if not idx:
            print('  %-11s 카르복실기 검출 실패' % name)
            continue
        print('  %-11s %2d원자 중 %2d개  (%s)  <- %s'
              % (name, len(at), len(idx), src.split(os.sep)[-1],
                 ' '.join(map(str, idx[:12]))))
        open(os.path.join(OUT, '%s.xyz' % name), 'w', encoding='utf-8',
             newline='\n').write('%d\n%s\n' % (len(at), name) +
                                 ''.join('%-3s %12.6f %12.6f %12.6f\n' % (e, *p)
                                         for e, p in at))
        open(os.path.join(OUT, '%s_pfreq.inp' % name), 'w', encoding='utf-8',
             newline='\n').write(FREQ.format(atoms=' '.join(map(str, idx)),
                                             q=q, name=name))
        names.append(name)

    open(os.path.join(OUT, 'species.txt'), 'w', encoding='utf-8',
         newline='\n').write(''.join(n + '\n' for n in names))

    runner = textwrap.dedent("""
        #!/bin/bash
        #SBATCH -J cit_ir
        #SBATCH -p long
        #SBATCH --array=1-NJ%NJ
        #SBATCH -c 1
        #SBATCH -t 0-02:00:00
        #SBATCH -o logs/slurm-%A_%a.out
        #SBATCH --requeue
        cd "$SLURM_SUBMIT_DIR" || exit 1
        [ -f ~/software/activate_orca_mpi.sh ] && source ~/software/activate_orca_mpi.sh
        N=$(sed -n "${SLURM_ARRAY_TASK_ID}p" species.txt); [ -z "$N" ] && exit 0
        say(){ echo "[$(date '+%H:%M:%S')] $N $*" >> STATUS.txt; }
        t0=$SECONDS
        timeout -k 30 4200 orca "${N}_pfreq.inp" > "${N}_pfreq.out" 2>&1
        grep -q "ORCA TERMINATED NORMALLY" "${N}_pfreq.out" \\
          && say "pfreq OK $((SECONDS-t0))s" || say "pfreq 실패 $((SECONDS-t0))s"
        rm -f "${N}_pfreq"*.tmp*
        say DONE
    """).lstrip().replace('NJ', str(max(1, len(names))))
    open(os.path.join(OUT, 'run_cit_ir.sh'), 'w', encoding='utf-8',
         newline='\n').write(runner)
    print('\nwrote %s/  (%d종)' % (OUT, len(names)))
    assert names, '생성된 종이 없다'


if __name__ == '__main__':
    main()
