"""Ni-시트르산 착물 ORCA 인풋 생성 (2026-09-24).

왜 필요한가:
  x(ChCl)=0.05 시료는 95%가 시트르산인데, 기존 화학종 세트(Cl/H2O)에 시트레이트가 없다.
  "산 과잉 쪽 우세종이 [Ni(H2O)6]2+ 가 맞느냐"는 정당한 반론이고, 심사에서 나올 질문이다.

무엇을 계산하나:
  [Ni(H3Cit)(H2O)4]2+  — 중성 시트르산이 중심 카르복실 C=O 와 중심 OH 로 이좌배위(5원환)
  H3Cit                — 자유 리간드 (반응식 양변 맞춤용)

  반응:  [Ni(H2O)6]2+ + H3Cit  ->  [Ni(H3Cit)(H2O)4]2+ + 2 H2O

  ** 중성 H3Cit 를 쓰는 이유 (중요) **
  시트레이트는 보통 탈양성자화해서 배위하지만, 그러면 반응식에 H+ 가 나타나고
  암시적 용매에서 양성자 기준상태는 악명 높게 부정확하다(수십 kJ/mol).
  우리 계는 산 과잉이라 중성 H3Cit 가 우세하므로, 양성자 이동 없는 배위로 두면
  **pKa/양성자 기준이 아예 필요 없다.** 정확도와 방어가능성 모두 이쪽이 낫다.
  (탈양성자화 배위는 더 강하므로, 이 계산은 시트르산 배위력의 **하한**이다 — 한계로 기재)

격리 이유:
  ~/dft_ni_complexes 의 잡배열이 실행 중이다. bash 는 스크립트를 실행 중에 조금씩 읽으므로
  run_array.sh 나 species.txt 를 건드리면 돌고 있는 잡이 깨진다.
  따라서 citrate/ 하위폴더에 독립된 species/러너를 만든다. 기존 파일은 절대 수정하지 않는다.

실행:
    cd ~/dft_ni_complexes
    source ~/venvs/cosmors/bin/activate     # rdkit 필요
    python build_citrate.py
    sbatch citrate/run_citrate.sh
"""
import math, os, sys, textwrap

R_NI_O = 2.06          # Ni-O(물) Angstrom
R_NI_O_LIG = 2.05      # Ni-O(시트르산)
R_OH, ANG_HOH = 0.97, 104.5
N_CONF = 300           # 킬레이트 가능한 배좌를 찾기 위한 후보 수
OO_MIN, OO_MAX = 2.3, 3.3   # 5원 킬레이트가 성립하는 O...O 거리 범위 (Angstrom)
CLEAR_MIN = 2.0        # Ni 와 비배위 원자의 최소 허용 거리

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'citrate')


# ------------------------------------------------------------ 벡터 보조
def sub(a, b): return tuple(x - y for x, y in zip(a, b))
def add(a, b): return tuple(x + y for x, y in zip(a, b))
def mul(a, s): return tuple(x * s for x in a)
def dot(a, b): return sum(x * y for x, y in zip(a, b))
def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
def norm(a): return math.sqrt(dot(a, a))
def unit(a):
    n = norm(a)
    if n < 1e-9:
        raise ValueError('영벡터 정규화')
    return mul(a, 1.0 / n)


def water(u, d):
    """Ni 에서 방향 u, 거리 d 에 물 배위 (O 가 Ni 를 향함)."""
    O = mul(u, d)
    p = (1., 0., 0.) if abs(u[0]) < 0.9 else (0., 1., 0.)
    v = unit(cross(p, u))
    h = math.radians(ANG_HOH / 2)
    out = [('O', O)]
    for s in (+1, -1):
        hv = add(mul(u, math.cos(h)), mul(v, s * math.sin(h)))
        out.append(('H', add(O, mul(hv, R_OH))))
    return out


# ------------------------------------------------- 시트르산 3D 구조 (RDKit)
def citric_acid():
    """RDKit 으로 시트르산 3D 좌표 + 배위에 쓸 두 산소의 인덱스를 찾는다.

    반환: (atoms, i_Ocarbonyl, i_Ohydroxyl)
      i_Ocarbonyl : 중심탄소에 붙은 COOH 의 카르보닐 산소 (C=O)
      i_Ohydroxyl : 중심탄소에 붙은 알코올 산소 (C-OH)
    이 둘이 Ni 와 5원 킬레이트 고리를 만든다 (알파-히드록시카르복실산의 전형적 배위).
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    smi = 'OC(=O)CC(O)(C(=O)O)CC(=O)O'          # citric acid
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    ps = AllChem.ETKDGv3()
    ps.randomSeed = 0xC17                        # 재현성
    # 배좌를 여럿 만든다. 자유산의 최저에너지 배좌는 카르복실이 엉뚱한 쪽을 보고 있어
    # 5원 킬레이트가 안 된다 — 실제로 Ni 를 놓을 자리가 없었음.
    # 킬레이트 가능한 배좌를 골라 쓰는 것이 맞다 (착물 형성 시 리간드는 재배열한다).
    if AllChem.EmbedMultipleConfs(m, numConfs=N_CONF, params=ps) == 0:
        raise SystemExit('RDKit 임베딩 실패')
    AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=2000)
    conf = m.GetConformer()

    # 주의 2가지 (둘 다 실제로 걸렸던 함정):
    #  1) AddHs 후에는 수소가 별도 원자라 GetTotalNumHs() 가 0 을 돌려준다.
    #     -> 수소 이웃을 직접 세야 한다.
    #  2) CH2 는 산소에 직접 붙어있지 않다. "산소를 가진 탄소 이웃 3개"로 찾으면 실패한다.
    nbrs = lambda a, sym: [n for n in a.GetNeighbors() if n.GetSymbol() == sym]

    # 중심탄소 = 수소 0개, 탄소이웃 3개(CH2 2 + COOH 1), 산소이웃 1개(알코올 OH)
    center = i_oh = None
    for a in m.GetAtoms():
        if a.GetSymbol() != 'C' or nbrs(a, 'H'):
            continue
        o, c = nbrs(a, 'O'), nbrs(a, 'C')
        if len(o) == 1 and len(c) >= 3 and nbrs(o[0], 'H'):
            center, i_oh = a, o[0].GetIdx()
            break
    if center is None:
        raise SystemExit('중심탄소를 찾지 못했습니다')

    # 중심탄소에 직접 붙은 카르복실기의 카르보닐 산소(C=O)
    i_oc = None
    for n in center.GetNeighbors():
        if n.GetSymbol() != 'C':
            continue
        for b in n.GetBonds():
            o = b.GetOtherAtom(n)
            if o.GetSymbol() == 'O' and b.GetBondType() == Chem.BondType.DOUBLE:
                i_oc = o.GetIdx()
        if i_oc is not None:
            break
    if i_oc is None:
        raise SystemExit('중심 카르복실의 C=O 산소를 찾지 못했습니다')

    confs = [[(a.GetSymbol(), tuple(c.GetAtomPosition(a.GetIdx()))) for a in m.GetAtoms()]
             for c in m.GetConformers()]
    return confs, i_oc, i_oh


def build_complex():
    """[Ni(H3Cit)(H2O)4]2+ 초기구조.

    배좌 N_CONF 개 x 원주 72지점을 전부 훑어, 배위산소 두 개에서 등거리이면서
    다른 원자와 가장 멀리 떨어지는 Ni 자리를 고른다. 거친 배치라도 XTB2 가 정리하지만,
    시작부터 원자가 겹쳐 있으면 XTB2 도 복구하지 못한다.
    """
    confs, i_oc, i_oh = citric_acid()
    others = [k for k in range(len(confs[0])) if k not in (i_oc, i_oh)]

    best = None                      # (clearance, 배좌, Ni좌표)
    n_ok = 0
    for lig in confs:
        O1, O2 = lig[i_oc][1], lig[i_oh][1]
        doo = norm(sub(O2, O1))
        if not (OO_MIN <= doo <= OO_MAX):        # 킬레이트 불가능한 배좌는 건너뜀
            continue
        n_ok += 1
        M = mul(add(O1, O2), 0.5)
        axis = unit(sub(O2, O1))
        half = doo / 2.0
        if half >= R_NI_O_LIG:
            continue
        h = math.sqrt(R_NI_O_LIG ** 2 - half ** 2)
        p = (1., 0., 0.) if abs(axis[0]) < 0.9 else (0., 1., 0.)
        b1 = unit(cross(axis, p))
        b2 = unit(cross(axis, b1))
        for deg in range(0, 360, 5):
            t = math.radians(deg)
            d = add(mul(b1, math.cos(t)), mul(b2, math.sin(t)))
            cand = add(M, mul(d, h))
            clear = min(norm(sub(cand, lig[k][1])) for k in others)
            if best is None or clear > best[0]:
                best = (clear, lig, cand)

    if best is None:
        raise SystemExit('킬레이트 가능한 배좌를 찾지 못했습니다 (O...O %.1f~%.1f A 조건)'
                         % (OO_MIN, OO_MAX))
    clear, lig, Ni = best
    print('  배좌 %d개 중 킬레이트 가능 %d개, 최선 배치의 여유거리 %.2f A' % (len(confs), n_ok, clear))
    if clear < CLEAR_MIN:
        raise SystemExit('Ni 를 놓을 여유가 부족합니다 (%.2f A < %.2f). 배위모드 재검토'
                         % (clear, CLEAR_MIN))
    O1, O2 = lig[i_oc][1], lig[i_oh][1]

    # Ni 기준 국소 팔면체 좌표계: e1~O1, e2~O2 평면, e3 = 법선
    u1 = unit(sub(O1, Ni))
    u2 = unit(sub(O2, Ni))
    e1 = u1
    e3 = unit(cross(u1, u2))
    e2 = unit(cross(e3, e1))
    free = [mul(e1, -1), mul(e2, -1), e3, mul(e3, -1)]           # 남은 4자리

    atoms = [('Ni', (0., 0., 0.))]
    atoms += [(sym, sub(pos, Ni)) for sym, pos in lig]           # Ni 를 원점으로 평행이동
    for u in free:
        atoms += water(u, R_NI_O)
    return atoms


def free_ligand():
    """자유 H3Cit. 착물과 같은 배좌를 쓸 이유가 없으므로 첫 번째(최저에너지급) 배좌를 쓴다."""
    confs, _, _ = citric_acid()
    lig = confs[0]
    cen = tuple(sum(p[1][k] for p in lig) / len(lig) for k in range(3))
    return [(sym, sub(pos, cen)) for sym, pos in lig]


# ------------------------------------------------------------ 건전성 점검
def sanity(name, atoms):
    """원자가 겹치거나 Ni 배위가 엉망이면 여기서 멈춘다. 잘못된 구조로 하루를 날리지 않기 위함."""
    n = len(atoms)
    worst = 9e9
    for i in range(n):
        for j in range(i + 1, n):
            d = norm(sub(atoms[i][1], atoms[j][1]))
            worst = min(worst, d)
            if d < 0.85:
                raise SystemExit('%s: 원자 %d(%s)-%d(%s) 거리 %.2f A — 구조 깨짐'
                                 % (name, i, atoms[i][0], j, atoms[j][0], d))
    ni = [i for i, a in enumerate(atoms) if a[0] == 'Ni']
    msg = '%-14s %2d atoms  최단거리 %.2f A' % (name, n, worst)
    if ni:
        d = sorted(norm(sub(atoms[ni[0]][1], a[1]))
                   for k, a in enumerate(atoms) if k != ni[0] and a[0] == 'O')
        coord = [x for x in d if x < 2.5]
        msg += '  Ni-O 배위 %d개 (%s)' % (len(coord), ', '.join('%.2f' % x for x in coord))
        if len(coord) != 6:
            raise SystemExit(msg + '  <- 팔면체 6배위가 아님')
    print('  ' + msg)
    return True


# ------------------------------------------------------------------- 인풋
HEAD = "%pal nprocs 1 end\n%maxcore 6000\n"

XTB = """! XTB2 Opt ALPB(water)
%geom MaxIter 800 end
""" + HEAD + """* xyzfile {q} {m} {name}_geo.xyz
"""
OPT = """! r2SCAN-3c TightOpt CPCM(water)
%geom MaxIter 600 end
%cpcm smd true SMDsolvent "water" end
""" + HEAD + """* xyzfile {q} {m} {name}_geo.xyz
"""
FREQ = """! r2SCAN-3c NumFreq CPCM(water)
%freq Temp 298.15, 333.15 end
%cpcm smd true SMDsolvent "water" end
""" + HEAD + """* xyzfile {q} {m} {name}_opt.xyz
"""

# (이름, 전하, 다중도, 원자, 단계)
def species():
    return [('NiCit_aq4', 2, 3, build_complex(), 'xtb,opt,freq'),
            ('H3Cit',     0, 1, free_ligand(),   'xtb,opt,freq')]


def w(fn, txt):
    open(os.path.join(OUT, fn), 'w', encoding='utf-8', newline='\n').write(txt)


def main():
    os.makedirs(os.path.join(OUT, 'logs'), exist_ok=True)
    rows = []
    print('구조 생성 및 점검:')
    for name, q, mult, at, steps in species():
        sanity(name, at)
        w('%s_start.xyz' % name,
          '%d\n%s\n' % (len(at), name) +
          ''.join('%-3s %14.8f %14.8f %14.8f\n' % (e, *xyz) for e, xyz in at))
        f = dict(q=q, m=mult, name=name)
        w(name + '_xtb.inp', XTB.format(**f))
        w(name + '_opt.inp', OPT.format(**f))
        w(name + '_freq.inp', FREQ.format(**f))
        rows.append('%s %s' % (name, steps))
    w('species.txt', '\n'.join(rows) + '\n')

    # 러너: ~/dft_ni_complexes/run_array.sh 와 같은 사다리 구조이나
    #  - 실행 중인 잡을 건드리지 않으려 파일을 분리했고
    #  - 분자가 3배 커서 시간상한(TMAX)과 walltime 을 크게 잡았다.
    w('run_citrate.sh', textwrap.dedent(r"""
        #!/bin/bash
        #SBATCH -J ni_cit
        #SBATCH -p long
        #SBATCH --array=1-2%2
        #SBATCH -c 1
        #SBATCH -t 4-00:00:00
        # --mem 금지: 이 노드는 RealMemory=1 이라 제출이 거부된다 (CLAUDE.md 참조).
        #SBATCH -o logs/slurm-%A_%a.out
        #SBATCH --requeue
        #
        # 시간예산 (최악): xtb 2x3600 + opt 3x86400 + freq 2x108000 = 484200s < 4일(345600s)?
        #   -> 아니므로 opt 3회가 아니라 2회, freq 1회 상한으로 조여 총 <4일이 되게 맞춘다:
        #      xtb 2x3600 + opt 3x72000 + freq 2x57600 = 338400s < 345600s. OK.
        # 정상 경로 추정: opt 7~15h + freq 15~25h = 하루 남짓.

        cd "$SLURM_SUBMIT_DIR" || exit 1
        [ -f ~/software/activate_orca_mpi.sh ] && source ~/software/activate_orca_mpi.sh

        LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" species.txt)
        N=${LINE%% *}
        [ -z "$N" ] && exit 0

        say() { echo "[$(date '+%m-%d %H:%M')] $N $*" >> STATUS.txt; }
        command -v orca >/dev/null || { say "FATAL orca not on PATH"; exit 1; }

        TMAX=0
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
            [ $? -eq 124 ] && say "$step TIMEOUT try=$try (상한 ${TMAX}s)"
            local dt=$((SECONDS-t0))
            if grep -q "ORCA TERMINATED NORMALLY" "${base}.out"; then
              say "$step OK try=$try ${dt}s ${extra:+[$extra]}"
              rm -f "${base}"*.tmp*
              return 0
            fi
            say "$step FAIL try=$try ${dt}s ${extra:+[$extra]}"
            [ -f "${base}_run.xyz" ] && cp "${base}_run.xyz" "${N}_geo.xyz"
            rm -f "${base}"*.tmp*
          done
          return 1
        }

        [ -f "${N}_geo.xyz" ] || cp "${N}_start.xyz" "${N}_geo.xyz"

        # 0단계 XTB2 사전최적화 — 33원자 유연한 분자라 이게 특히 중요하다
        TMAX=3600
        if [ ! -f "${N}.xtb.done" ]; then
          if attempt xtb "${N}_xtb" "" "SlowConv"; then
            [ -f "${N}_xtb_run.xyz" ] && cp "${N}_xtb_run.xyz" "${N}_geo.xyz"
            touch "${N}.xtb.done"
          else
            say "xtb 실패 -> 초기구조 그대로 진행"
          fi
        fi

        TMAX=72000
        if [ ! -f "${N}.opt.done" ]; then
          attempt opt "${N}_opt" "" "SlowConv" "SlowConv NormalOpt"
          if [ -f "${N}_opt_run.xyz" ]; then
            cp "${N}_opt_run.xyz" "${N}_opt.xyz"
            [ -f "${N}_opt_run.gbw" ] && cp "${N}_opt_run.gbw" "${N}_opt.gbw"
            grep -q "ORCA TERMINATED NORMALLY" "${N}_opt.out" && touch "${N}.opt.done" \
              || say "opt 미수렴 -> 마지막 기하 사용 (한계로 기록)"
          else
            cp "${N}_geo.xyz" "${N}_opt.xyz"; say "opt 결과 없음 -> XTB 기하로 대체"
          fi
        fi

        TMAX=57600
        if [ ! -f "${N}.freq.done" ]; then
          attempt freq "${N}_freq" "" "SlowConv" && touch "${N}.freq.done" \
            || say "freq 실패 -> 이 종은 dG 없음"
        fi

        say "DONE"
    """).lstrip())

    print('\nwrote %s/' % OUT)
    print('제출:  sbatch citrate/run_citrate.sh      확인:  cat citrate/STATUS.txt')


if __name__ == '__main__':
    main()
