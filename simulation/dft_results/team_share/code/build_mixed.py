"""Mixed ligand [MCl2(H3Cit)(H2O)2] — 2시간 예산 경향성 계산 (2026-09-30).

왜:
  민이형/GPT 지적 — ChCl:CA 비율을 바꾸면 Cl- 와 citrate 비율이 **동시에** 바뀐다.
  순수 Cl 사다리와 순수 citrate 착물만 따로 계산하면 혼합 배위종을 놓친다.
  그리고 UV-Vis 345 nm 밴드는 d-d 가 아니라 LMCT 영역이라, citrate 가 붙은
  화학종의 지문일 수 있다. 그래서 에너지와 UV-Vis 를 같이 본다.

속도 설계 (요구: 2시간 내, 정확도보다 경향성):
  * 초기구조를 공짜로 얻는다 — 이미 DFT 최적화된 [Ni(H3Cit)(H2O)4] 에서
    물 2개를 Cl 로 치환. 시트르산 배위 모티프가 이미 수렴된 상태라 훨씬 빨리 붙는다.
  * freq 를 돌리지 않는다. 30원자 NumFreq 는 몇 시간이다. 따라서 dG 가 아니라 dE 다.
    ** 반응이 2분자 -> 3분자라 dS > 0 이므로 dG < dE. 즉 dE 는 citrate 결합을
       과소평가한다(보수적). dE 로도 citrate 가 이기면 결론이 견고하다. **
  * 단계마다 결과를 파일로 남긴다. 중간에 시간이 끊겨도 그때까지의 답이 있다.

비교 반응 (원자·전하 balanced):
  [MCl2(H2O)4] + H3Cit  ->  [MCl2(H3Cit)(H2O)2] + 2 H2O
  좌변 3종은 전부 기존 계산에 있다.

실행:
    python build_mixed.py
    sbatch mixed/run_mixed.sh
"""
import math
import os
import sys
import textwrap

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(os.path.dirname(HERE), 'dft_results', 'geometries',
                   'citrate', 'NiCit_aq4_opt.xyz')
OUT = os.path.join(HERE, 'mixed')

# 다중도는 기존 계산과 동일 (고스핀). r_cl 은 기존 수렴값.
METALS = {'Ni': dict(mult=3, r_cl=2.37), 'Co': dict(mult=4, r_cl=2.39),
          'Mn': dict(mult=6, r_cl=2.43)}
N_CL = 2                      # 물 2개를 Cl 로 치환 -> 전하 0


def read_xyz(p):
    L = [l.split() for l in open(p, encoding='utf-8') if l.strip()]
    n = int(L[0][0])
    return [(r[0], tuple(float(v) for v in r[1:4])) for r in L[2:2 + n]]


def dist(a, b):
    return math.dist(a, b)


def build(metal):
    at = read_xyz(SRC)
    mi = next(i for i, a in enumerate(at) if a[0] == 'Ni')
    m = at[mi][1]

    # 물 = H 2개를 달고 있는 O. 시트르산 O 와 구분해야 한다.
    waters = []
    for i, (e, p) in enumerate(at):
        if e != 'O':
            continue
        hs = [j for j, (e2, p2) in enumerate(at) if e2 == 'H' and dist(p, p2) < 1.15]
        if len(hs) == 2:
            waters.append((dist(m, p), i, hs))
    waters.sort(reverse=True)                 # 금속에서 먼 물부터 치환
    assert len(waters) >= N_CL, '배위된 물이 %d개뿐' % len(waters)

    drop, cl_dirs = set(), []
    for _, oi, hs in waters[:N_CL]:
        drop.add(oi); drop.update(hs)
        u = [(c - d) for c, d in zip(at[oi][1], m)]
        n = math.sqrt(sum(v * v for v in u))
        cl_dirs.append([v / n for v in u])    # 그 물이 있던 방향에 Cl 을 놓는다

    r = METALS[metal]['r_cl']
    new = [(metal if i == mi else e, p) for i, (e, p) in enumerate(at) if i not in drop]
    for u in cl_dirs:
        new.append(('Cl', tuple(c + v * r for c, v in zip(m, u))))
    return new


HEAD = "%pal nprocs 1 end\n%maxcore 6000\n"
XTB = """! XTB2 Opt ALPB(water)
%geom MaxIter 400 end
""" + HEAD + """* xyzfile 0 {m} {name}_geo.xyz
"""
# LooseOpt + MaxIter 제한: 정확도보다 시간. 경향성만 보면 되므로.
OPT = """! r2SCAN-3c LooseOpt CPCM(water)
%geom MaxIter 120 end
%cpcm smd true SMDsolvent "water" end
""" + HEAD + """* xyzfile 0 {m} {name}_geo.xyz
"""
SP = """! r2SCAN-3c CPCM(water)
%cpcm smd true SMDsolvent "water" end
""" + HEAD + """* xyzfile 0 {m} {name}_sp.xyz
"""
# UV-Vis: 345 nm 는 LMCT 영역이라 TD-DFT 가 쓸 만하다(d-d 와 달리).
# 하이브리드 PBE0 + 작은 기저 + TDA. 절대 파장이 아니라 밴드 유무만 본다.
TDA = """! PBE0 def2-SVP def2/J RIJCOSX CPCM(water)
""" + HEAD + """%tddft
  nroots 25
  maxdim 250
  tda true
end
* xyzfile 0 {m} {name}_sp.xyz
"""


def w(fn, txt):
    open(os.path.join(OUT, fn), 'w', encoding='utf-8', newline='\n').write(txt)


def main():
    os.makedirs(os.path.join(OUT, 'logs'), exist_ok=True)
    names = []
    print('혼합 배위종 생성 (기존 최적화 구조에서 물->Cl 치환):')
    for metal, M in METALS.items():
        at = build(metal)
        name = '%sCl2_H3Cit_aq2' % metal
        mi = next(i for i, a in enumerate(at) if a[0] == metal)
        d = sorted(dist(at[mi][1], p) for e, p in at
                   if e in ('O', 'Cl') and p != at[mi][1])
        ncoord = len([x for x in d if x < 2.7])
        print('  %-18s %2d원자  전하 0  mult %d  배위 %d개 (%s)'
              % (name, len(at), M['mult'], ncoord,
                 ', '.join('%.2f' % x for x in d[:6])))
        assert ncoord == 6, '%s: 6배위가 아님 (%d)' % (name, ncoord)
        w('%s_start.xyz' % name, '%d\n%s\n' % (len(at), name) +
          ''.join('%-3s %14.8f %14.8f %14.8f\n' % (e, *p) for e, p in at))
        f = dict(m=M['mult'], name=name)
        for tag, tpl in (('xtb', XTB), ('opt', OPT), ('sp', SP), ('tda', TDA)):
            w('%s_%s.inp' % (name, tag), tpl.format(**f))
        names.append(name)

    w('species.txt', ''.join(n + '\n' for n in names))

    runner = textwrap.dedent("""
        #!/bin/bash
        #SBATCH -J mixed
        #SBATCH -p long
        #SBATCH --array=1-NJ%NJ
        #SBATCH -c 1
        #SBATCH -t 0-04:00:00
        # --mem 금지 (RealMemory=1). 메모리는 %maxcore 로만.
        #SBATCH -o logs/slurm-%A_%a.out
        #SBATCH --requeue

        cd "$SLURM_SUBMIT_DIR" || exit 1
        [ -f ~/software/activate_orca_mpi.sh ] && source ~/software/activate_orca_mpi.sh
        N=$(sed -n "${SLURM_ARRAY_TASK_ID}p" species.txt)
        [ -z "$N" ] && exit 0
        say() { echo "[$(date '+%H:%M:%S')] $N $*" >> STATUS.txt; }
        command -v orca >/dev/null || { say "FATAL orca 없음"; exit 1; }

        run() {  # run <단계> <타임아웃초>
          local s=$1 t=$2
          [ -f "${N}.$s.done" ] && return 0
          local t0=$SECONDS
          timeout -k 30 "$t" orca "${N}_${s}.inp" > "${N}_${s}.out" 2>&1
          local dt=$((SECONDS-t0))
          if grep -q "ORCA TERMINATED NORMALLY" "${N}_${s}.out"; then
            say "$s OK ${dt}s"; touch "${N}.$s.done"; rm -f "${N}_${s}"*.tmp*; return 0
          fi
          say "$s 미완 ${dt}s (시간초과 또는 실패)"; rm -f "${N}_${s}"*.tmp*; return 1
        }

        cp "${N}_start.xyz" "${N}_geo.xyz"

        # 1) xTB 로 빠르게 다듬는다
        run xtb 600 && [ -f "${N}_xtb.xyz" ] && cp "${N}_xtb.xyz" "${N}_geo.xyz"

        # 2) xTB 기하로 먼저 SP — 여기서 이미 '답' 이 하나 확보된다.
        #    뒤 단계가 시간에 걸려 죽어도 빈손이 아니게 하는 장치.
        cp "${N}_geo.xyz" "${N}_sp.xyz"
        run sp 1200 && cp "${N}_sp.out" "${N}_sp_xtbgeom.out" && rm -f "${N}.sp.done"

        # 3) DFT 로 느슨하게 최적화 (1시간 상한)
        if run opt 3600; then
          [ -f "${N}_opt.xyz" ] && cp "${N}_opt.xyz" "${N}_sp.xyz"
        else
          # 미수렴이어도 마지막 기하가 xTB 보다는 낫다
          [ -f "${N}_opt.xyz" ] && cp "${N}_opt.xyz" "${N}_sp.xyz" && say "미수렴 기하 사용"
        fi

        # 4) 개선된 기하로 다시 SP
        run sp 1200

        # 5) UV-Vis (345 nm 확인). 실패해도 에너지 결론에는 영향 없음.
        run tda 2400 || say "tda 생략 -> UV-Vis 비교 불가"
        say DONE
    """).lstrip().replace('NJ', str(len(names)))
    w('run_mixed.sh', runner)

    print('\nwrote %s/  (%d종)' % (OUT, len(names)))
    print('최악 소요: xtb 10분 + sp 20분 + opt 60분 + sp 20분 + tda 40분 = 약 2.5시간')
    print('제출: sbatch mixed/run_mixed.sh')

    for n in names:
        for s in ('xtb', 'opt', 'sp', 'tda'):
            p = os.path.join(OUT, '%s_%s.inp' % (n, s))
            assert os.path.exists(p) and '{' not in open(p, encoding='utf-8').read()
    sh = open(os.path.join(OUT, 'run_mixed.sh'), encoding='utf-8').read()
    assert not any('--mem' in l for l in sh.splitlines()
                   if l.strip().startswith('#SBATCH'))
    assert '\r' not in sh
    print('self-check OK: 인풋 %d개, --mem 없음, LF 줄끝' % (len(names) * 4))


if __name__ == '__main__':
    main()
