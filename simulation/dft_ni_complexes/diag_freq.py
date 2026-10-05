"""허수진동수 원인 진단.

확인할 것:
 1. opt 가 실제로 수렴했는가 (HURRAY)
 2. freq 가 opt 와 같은 기하·같은 전자상태를 썼는가
    -> opt 최종 SCF 에너지 vs freq 첫 SCF 에너지가 같아야 한다. 다르면
       (a) 기하가 다르거나 (b) SCF 가 다른 해로 수렴한 것.
 3. 허수모드의 크기 분포
"""
import glob, os, re, sys

def last_scf(path):
    t = open(path, encoding='utf-8', errors='ignore').read()
    h = re.findall(r'FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)', t)
    return (float(h[0]) if h else None, float(h[-1]) if h else None, len(h))

def converged(path):
    t = open(path, encoding='utf-8', errors='ignore').read()
    return 'HURRAY' in t or 'THE OPTIMIZATION HAS CONVERGED' in t

def freqs(path):
    t = open(path, encoding='utf-8', errors='ignore').read()
    m = re.search(r'VIBRATIONAL FREQUENCIES(.*?)(NORMAL MODES|$)', t, re.S)
    if not m:
        return []
    return [float(x) for x in re.findall(r':\s+(-?\d+\.\d+)\s*cm\*\*-1', m.group(1))]

def geom_rmsd(a, b):
    """두 xyz 파일의 좌표 차이 (원자 순서 같다고 가정)."""
    def rd(p):
        L = [l.split() for l in open(p).read().splitlines()[2:] if l.strip()]
        return [(x[0], tuple(map(float, x[1:4]))) for x in L]
    A, B = rd(a), rd(b)
    if len(A) != len(B):
        return None
    s = sum(sum((p[1][k] - q[1][k]) ** 2 for k in range(3)) for p, q in zip(A, B))
    return (s / len(A)) ** 0.5

here = os.path.dirname(os.path.abspath(__file__))
targets = sorted(glob.glob(os.path.join(here, '*_freq.out'))) + \
          sorted(glob.glob(os.path.join(here, 'citrate', '*_freq.out')))

print('%-13s %-9s %-16s %-16s %-9s %s' %
      ('종', 'opt수렴', 'opt 최종 E', 'freq 첫 E', 'ΔE(Eh)', '허수(cm-1)'))
for f in targets:
    d = os.path.dirname(f)
    n = os.path.basename(f).replace('_freq.out', '')
    o = os.path.join(d, n + '_opt.out')
    if not os.path.exists(o):
        continue
    _, e_opt, _ = last_scf(o)
    e_frq, _, _ = last_scf(f)
    v = [x for x in freqs(f) if x < -1]
    de = (e_frq - e_opt) if (e_opt and e_frq) else None
    print('%-13s %-9s %-16.8f %-16.8f %-9s %s' % (
        n, 'YES' if converged(o) else 'NO', e_opt, e_frq,
        ('%+.6f' % de) if de is not None else '?',
        ', '.join('%.0f' % abs(x) for x in sorted(v)[:6]) or '없음'))

# freq 가 읽은 기하와 opt 결과 기하가 같은지
print()
for f in targets:
    d = os.path.dirname(f)
    n = os.path.basename(f).replace('_freq.out', '')
    a = os.path.join(d, n + '_opt.xyz')
    b = os.path.join(d, n + '_opt_run.xyz')
    if os.path.exists(a) and os.path.exists(b):
        r = geom_rmsd(a, b)
        print('  %-13s freq 입력기하 vs opt 산출기하 RMSD = %s' %
              (n, ('%.4f A' % r) if r is not None else '원자수 불일치!'))
