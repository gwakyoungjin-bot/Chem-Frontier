"""시트르산 C=O 를 누가 움직이는가 — 환경 네 가지를 동시에 계산한다.

실측: CA 많은 쪽 1713.9 -> ChCl 많은 쪽 1729.2 (+15.3 cm-1)
       CA 가 많으면 낮고, ChCl 이 많으면 높다.

각 환경이 예측하는 방향:
  H3Cit_dimer   COOH...COOH 자기회합. CA 많은 쪽을 대표한다.
                C=O 가 수소결합을 받으므로 **가장 낮아야** 한다.
  H3Cit_chol    콜린 양이온의 -OH 가 C=O 에 수소결합.
                ACS Omega 2023 이 이 상호작용을 보고했다. ChCl 이 많을수록 늘어난다.
                이것도 C=O 를 **내린다**면 실측(ChCl 많을수록 높음)과 반대가 되어
                콜린 가설이 기각된다.
  H3Cit_w2      물 2개가 붙은 상태. 배경.
  H3Cit_Cl1/2   별도 잡(cit_ir)에서 계산 중.

즉 이 계산은 "무엇이 낮은 쪽 끝(1713.9)을 만드는가"를 가린다. 그게 이량체면
FT-IR 은 CA 자기회합을 보는 것이고, 콜린이면 ChCl 쪽이 낮아야 하므로 기각된다.

속도: 부분 헤시안으로 카르복실기만 푼다. opt 는 1시간 상한을 건다.

실행 (서버): python build_cit_env.py && sbatch cit_env/run_cit_env.sh
"""
import math
import os
import re
import sys
import textwrap

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, 'cit_env')
CIT = os.path.join(HERE, 'clsat', 'H3Cit_opt.xyz')
CHCL_COSMO = [os.path.join(ROOT, 'feature_extraction', 'orcacosmo',
                           'choline_chloride', 'choline_chloride_c000.orcacosmo'),
              os.path.join(os.path.expanduser('~'), 'software', 'cosmors_offline',
                           'openCOSMO-RS_conformer_pipeline', 'choline_chloride',
                           'COSMO_TZVPD', 'choline_chloride_c000.orcacosmo')]
BOHR = 0.529177210903


def read_xyz(p):
    L = [l.split() for l in open(p, encoding='utf-8') if l.strip()]
    return [(r[0], tuple(float(v) for v in r[1:4])) for r in L[2:2 + int(L[0][0])]]


def read_cosmo_xyz(p):
    """orcacosmo 의 #XYZ_FILE 절에서 좌표를 꺼낸다 (Angstrom)."""
    lines = open(p, encoding='utf-8', errors='replace').read().splitlines()
    i = next(k for k, l in enumerate(lines) if l.startswith('#XYZ_FILE'))
    n = int(lines[i + 1].split()[0])
    out = []
    for l in lines[i + 3:i + 3 + n]:
        q = l.split()
        out.append((q[0], tuple(float(v) for v in q[1:4])))
    return out


def d(a, b): return math.dist(a, b)
def sub(a, b): return tuple(x - y for x, y in zip(a, b))
def add(a, b): return tuple(x + y for x, y in zip(a, b))
def mul(a, s): return tuple(x * s for x in a)
def unit(a): n = math.sqrt(sum(x * x for x in a)); return mul(a, 1.0 / n)


def carboxyl(at):
    idx = set()
    for i, (e, p) in enumerate(at):
        if e != 'C':
            continue
        os_ = [j for j, (e2, q) in enumerate(at) if e2 == 'O' and d(p, q) < 1.60]
        if len(os_) < 2:
            continue
        idx.add(i); idx.update(os_)
        for j in os_:
            for k, (e3, r) in enumerate(at):
                if e3 == 'H' and d(at[j][1], r) < 1.15:
                    idx.add(k)
    return sorted(idx)


def acidic_h(at):
    """(H index, O index) 목록. 카르복실 O 에 붙은 H."""
    out = []
    for i, (e, p) in enumerate(at):
        if e != 'H':
            continue
        for j, (e2, q) in enumerate(at):
            if e2 == 'O' and d(p, q) < 1.15:
                for k, (e3, s) in enumerate(at):
                    if e3 == 'C' and d(q, s) < 1.6 and \
                       sum(1 for m, (e4, t) in enumerate(at)
                           if e4 == 'O' and d(s, t) < 1.6) >= 2:
                        out.append((i, j)); break
                break
    return out


def carbonyl_o(at):
    """C=O 산소 (H 가 안 붙은 카르복실 산소) 하나."""
    for i, (e, p) in enumerate(at):
        if e != 'O':
            continue
        has_h = any(e2 == 'H' and d(p, q) < 1.15 for e2, q in at)
        near_c = any(e2 == 'C' and d(p, q) < 1.6 for e2, q in at)
        if near_c and not has_h:
            return i
    return None


HEAD = "%pal nprocs 1 end\n%maxcore 6000\n"
OPT = """! r2SCAN-3c LooseOpt CPCM(water)
%geom MaxIter 120 end
%cpcm smd true SMDsolvent "water" end
""" + HEAD + """* xyzfile {q} 1 {name}_geo.xyz
"""
PF = """! r2SCAN-3c NumFreq CPCM(water)
%cpcm smd true SMDsolvent "water" end
%freq
  PARTIAL_Hess true
  Hess_Atoms {{{atoms}}} end
end
""" + HEAD + """* xyzfile {q} 1 {name}_fin.xyz
"""


def main():
    if not os.path.exists(CIT):
        sys.exit('H3Cit_opt.xyz 없음 (clsat 먼저)')
    cit = read_xyz(CIT)
    ah = acidic_h(cit)
    co = carbonyl_o(cit)
    sp = {}

    # 1) 이량체: 산성 H 와 상대 분자의 C=O 가 1.75 A 가 되도록 반전 배치
    hi, oi = ah[0]
    u = unit(sub(cit[hi][1], cit[oi][1]))
    acc = add(cit[hi][1], mul(u, 1.75))
    cen = tuple((a + b) / 2 for a, b in zip(acc, cit[co][1]))
    second = [(e, tuple(2 * c - x for c, x in zip(cen, p))) for e, p in cit]
    sp['H3Cit_dimer'] = (cit + second, 0)

    # 2) 물 2개를 서로 다른 산성 H 에
    w = list(cit)
    for k in (0, 1):
        h2, o2 = ah[k]
        v = unit(sub(cit[h2][1], cit[o2][1]))
        ow = add(cit[h2][1], mul(v, 1.75))
        w += [('O', ow), ('H', add(ow, (0.60, 0.76, 0.05))),
              ('H', add(ow, (0.55, -0.58, 0.58)))]
    sp['H3Cit_w2'] = (w, 0)

    # 3) 콜린 양이온: ChCl 이온쌍에서 Cl 을 빼고 C=O 쪽에 붙인다
    chp = next((p for p in CHCL_COSMO if os.path.exists(p)), None)
    if chp:
        pair = read_cosmo_xyz(chp)
        chol = [(e, p) for e, p in pair if e != 'Cl']
        # 콜린의 -OH 수소를 찾아 그것을 시트르산 C=O 쪽으로 1.75 A 에 둔다
        oh = [(i, j) for i, (e, p) in enumerate(chol) for j, (e2, q) in enumerate(chol)
              if e == 'H' and e2 == 'O' and d(p, q) < 1.15]
        if oh:
            hi2, oi2 = oh[0]
            target = add(cit[co][1], mul(unit(sub(cit[co][1], cit[co - 1][1])), 1.75))
            shift = sub(target, chol[hi2][1])
            chol = [(e, add(p, shift)) for e, p in chol]
            sp['H3Cit_chol'] = (cit + chol, 1)
    else:
        print('  (콜린 프로파일 없음 — H3Cit_chol 생략)')

    os.makedirs(os.path.join(OUT, 'logs'), exist_ok=True)
    names = []
    for name, (at, q) in sp.items():
        n = len(at)
        worst = min(d(at[i][1], at[j][1]) for i in range(n) for j in range(i + 1, n))
        if worst < 0.80:
            print('  %-14s 원자 겹침 %.2f A — 건너뜀' % (name, worst)); continue
        idx = carboxyl(at[:len(cit)])          # 첫 시트르산의 카르복실만
        print('  %-14s %2d원자 전하 %+d 최단 %.2f A  헤시안 %d개' % (name, n, q, worst, len(idx)))
        open(os.path.join(OUT, '%s_start.xyz' % name), 'w', encoding='utf-8',
             newline='\n').write('%d\n%s\n' % (n, name) +
                                 ''.join('%-3s %12.6f %12.6f %12.6f\n' % (e, *p) for e, p in at))
        f = dict(q=q, name=name, atoms=' '.join(map(str, idx)))
        open(os.path.join(OUT, '%s_opt.inp' % name), 'w', encoding='utf-8',
             newline='\n').write(OPT.format(**f))
        open(os.path.join(OUT, '%s_pfreq.inp' % name), 'w', encoding='utf-8',
             newline='\n').write(PF.format(**f))
        names.append(name)
    open(os.path.join(OUT, 'species.txt'), 'w', encoding='utf-8',
         newline='\n').write(''.join(n + '\n' for n in names))

    runner = textwrap.dedent("""
        #!/bin/bash
        #SBATCH -J cit_env
        #SBATCH -p long
        #SBATCH --array=1-NJ%NJ
        #SBATCH -c 1
        #SBATCH -t 0-03:00:00
        #SBATCH -o logs/slurm-%A_%a.out
        #SBATCH --requeue
        cd "$SLURM_SUBMIT_DIR" || exit 1
        [ -f ~/software/activate_orca_mpi.sh ] && source ~/software/activate_orca_mpi.sh
        N=$(sed -n "${SLURM_ARRAY_TASK_ID}p" species.txt); [ -z "$N" ] && exit 0
        say(){ echo "[$(date '+%H:%M:%S')] $N $*" >> STATUS.txt; }
        cp "${N}_start.xyz" "${N}_geo.xyz"; cp "${N}_geo.xyz" "${N}_fin.xyz"
        t0=$SECONDS
        timeout -k 30 3600 orca "${N}_opt.inp" > "${N}_opt.out" 2>&1
        if grep -q "ORCA TERMINATED NORMALLY" "${N}_opt.out"; then
          say "opt OK $((SECONDS-t0))s"; [ -f "${N}_opt.xyz" ] && cp "${N}_opt.xyz" "${N}_fin.xyz"
        else
          say "opt 미완 $((SECONDS-t0))s"; [ -f "${N}_opt.xyz" ] && cp "${N}_opt.xyz" "${N}_fin.xyz"
        fi
        rm -f "${N}_opt"*.tmp*
        t0=$SECONDS
        timeout -k 30 4200 orca "${N}_pfreq.inp" > "${N}_pfreq.out" 2>&1
        grep -q "ORCA TERMINATED NORMALLY" "${N}_pfreq.out" \\
          && say "pfreq OK $((SECONDS-t0))s" || say "pfreq 실패 $((SECONDS-t0))s"
        rm -f "${N}_pfreq"*.tmp*; say DONE
    """).lstrip().replace('NJ', str(max(1, len(names))))
    open(os.path.join(OUT, 'run_cit_env.sh'), 'w', encoding='utf-8',
         newline='\n').write(runner)
    print('\nwrote %s/  (%d종)' % (OUT, len(names)))
    assert names


if __name__ == '__main__':
    main()
