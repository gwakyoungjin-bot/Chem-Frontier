"""15종 착물을 한 화면에 3x5 격자로 배치한 xyz 하나를 만든다.

왜 다중프레임이 안 됐나:
  xyz 다중프레임은 **모든 프레임의 원자 수가 같을 때만** 트래젝토리로 읽힌다.
  우리 종은 19/17/15/7/5 로 제각각이라 Avogadro 가 첫 프레임만 읽고 멈춘다.
  그래서 '넘겨보기' 가 아니라 '한 번에 늘어놓기' 로 간다.

배치:  행 = 금속 (Ni / Co / Mn),  열 = Cl 개수 (0 -> 4)
  왼쪽에서 오른쪽으로 갈수록 Cl 이 하나씩 늘고, 3열(Cl2)과 4열(Cl3) 사이에서
  팔면체 -> 사면체로 꺾인다. 그 꺾임이 이 그림의 요점이다.

간격 SPACING 은 12 A. Avogadro 는 공유결합 반지름으로 결합을 자동 추론하므로
이웃 착물끼리 가까우면 없는 결합이 그려진다. 12 A 면 확실히 안 붙는다.

실행: python dft_ni_complexes/export_grid.py
"""
import os
import sys

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

HERE = os.path.dirname(os.path.abspath(__file__))
GEO = os.path.join(os.path.dirname(HERE), 'dft_results', 'geometries')
OPT = os.path.join(GEO, 'optimized')

ROWS = [('Ni', ['Ni_aq6', 'NiCl1_aq5', 'NiCl2_aq4', 'NiCl3_aq1', 'NiCl4']),
        ('Co', ['Co_aq6', 'CoCl1_aq5', 'CoCl2_aq4', 'CoCl3_aq1', 'CoCl4']),
        ('Mn', ['Mn_aq6', 'MnCl1_aq5', 'MnCl2_aq4', 'MnCl3_aq1', 'MnCl4'])]
SPACING = 12.0
METALS = ('Ni', 'Co', 'Mn')


def read_xyz(path):
    with open(path, encoding='utf-8') as f:
        lines = [l.rstrip('\n') for l in f if l.strip()]
    n = int(lines[0].split()[0])
    at = []
    for l in lines[2:2 + n]:
        p = l.split()
        at.append((p[0], float(p[1]), float(p[2]), float(p[3])))
    return at


def centre_on_metal(at):
    """금속 원자를 원점으로. 격자 칸마다 금속이 정렬돼야 비교가 쉽다."""
    for e, x, y, z in at:
        if e in METALS:
            return [(el, a - x, b - y, c - z) for el, a, b, c in at]
    cx = sum(a for _, a, _, _ in at) / len(at)      # 금속 없으면 무게중심
    cy = sum(b for _, _, b, _ in at) / len(at)
    cz = sum(c for _, _, _, c in at) / len(at)
    return [(el, a - cx, b - cy, c - cz) for el, a, b, c in at]


def main():
    atoms, placed, missing = [], [], []
    for r, (metal, names) in enumerate(ROWS):
        for c, name in enumerate(names):
            src = os.path.join(OPT, '%s_opt.xyz' % name)
            if not os.path.exists(src):
                missing.append(name)
                continue
            at = centre_on_metal(read_xyz(src))
            dx, dy = c * SPACING, -r * SPACING
            atoms += [(e, x + dx, y + dy, z) for e, x, y, z in at]
            placed.append((r, c, name, len(at)))

    out = os.path.join(GEO, 'grid_all15.xyz')
    with open(out, 'w', encoding='ascii', newline='\r\n') as f:
        f.write('%d\n' % len(atoms))
        f.write('15 complexes, rows Ni/Co/Mn, cols Cl0-Cl4, %.0fA grid, '
                'r2SCAN-3c/CPCM optimized\n' % SPACING)
        for e, x, y, z in atoms:
            f.write('%-3s %14.8f %14.8f %14.8f\n' % (e, x, y, z))

    for r, c, name, n in placed:
        print('  행%d 열%d  %-11s %2d atoms' % (r + 1, c + 1, name, n))
    if missing:
        print('  빠짐: %s' % ', '.join(missing))
    print('\nwrote %s' % out)
    print('원자 %d개, 착물 %d종' % (len(atoms), len(placed)))

    # 자체검증: 이웃 착물이 붙어 보이면 그림이 거짓말을 한다.
    # 서로 다른 칸의 원자 간 최소거리가 결합 판정 한계(~2.2 A)보다 커야 한다.
    import itertools
    cent = [(c * SPACING, -r * SPACING) for r, c, _, _ in placed]
    worst = min(((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) ** 0.5
                for a, b in itertools.combinations(cent, 2))
    assert worst >= 10.0, '격자 간격이 너무 좁다: %.1f A' % worst
    print('self-check OK: 칸 사이 최소 %.1f A (결합 오인 없음)' % worst)


if __name__ == '__main__':
    main()
