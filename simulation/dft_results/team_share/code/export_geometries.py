"""DFT 구조 17종을 Avogadro 로 한 번에 볼 수 있게 묶는다.

왜: 구조가 금속별 폴더에 흩어져 있어 "사다리 전체"를 한눈에 못 본다.
    다중프레임 xyz 로 묶으면 Avogadro 의 애니메이션 슬라이더로 넘겨볼 수 있다.

프레임 순서는 배위 사다리 순서다 (아쿠아 -> 클로로).  같은 금속 안에서 Cl 이
하나씩 늘어나고, aq6/Cl1/Cl2 는 팔면체, Cl3/Cl4 는 사면체로 넘어간다.

두 벌을 만든다:
  optimized/  서버에서 받은 **DFT 최적화 완료** 좌표 (r2SCAN-3c/CPCM(SMD water)).
              run_array.sh:99-101 이 ORCA 출력을 <종>_opt.xyz 로 복사한 것.
  (루트)      최적화 **전** 초기 구조. 최적화가 무엇을 바꿨는지 비교용.

프레임 제목에 허수진동수 합격 여부를 넣는다 — 불합격 종은 그 구조가 극소점이
아니라 안장점이므로 구조를 인용할 때 반드시 같이 말해야 한다(CO_REPORT.txt [1]).

실행: python dft_ni_complexes/export_geometries.py
"""
import os
import shutil
import sys

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

HERE = os.path.dirname(os.path.abspath(__file__))
GEO = os.path.join(os.path.dirname(HERE), 'dft_results', 'geometries')
OPT = os.path.join(GEO, 'optimized')

# (하위폴더, 종이름, 전하, 스핀다중도)
# 다중도: Ni(II) d8 -> 3,  Co(II) d7 -> 4,  Mn(II) d5 고스핀 -> 6
LADDER = [
    ('',          'Ni_aq6',    +2, 3), ('',          'NiCl1_aq5', +1, 3),
    ('',          'NiCl2_aq4',  0, 3), ('',          'NiCl3_aq1', -1, 3),
    ('',          'NiCl4',     -2, 3),
    ('cobalt',    'Co_aq6',    +2, 4), ('cobalt',    'CoCl1_aq5', +1, 4),
    ('cobalt',    'CoCl2_aq4',  0, 4), ('cobalt',    'CoCl3_aq1', -1, 4),
    ('cobalt',    'CoCl4',     -2, 4),
    ('manganese', 'Mn_aq6',    +2, 6), ('manganese', 'MnCl1_aq5', +1, 6),
    ('manganese', 'MnCl2_aq4',  0, 6), ('manganese', 'MnCl3_aq1', -1, 6),
    ('manganese', 'MnCl4',     -2, 6),
]
EXTRA = [('', 'H2O', 0, 1), ('', 'Cl_ion', -1, 1)]

# CO_REPORT.txt [1] 의 잔여 허수진동수 (cm-1). 합격 기준은 |값| <= 50.
IMAG = {'Co_aq6': 142, 'CoCl1_aq5': 125,
        'Mn_aq6': 111, 'MnCl2_aq4': 65, 'MnCl1_aq5': 136,
        'NiCl1_aq5': 28}

GEOM = 'Oh(팔면체)'
SHAPE = {'aq6': GEOM, 'Cl1_aq5': GEOM, 'Cl2_aq4': GEOM,
         'Cl3_aq1': 'Td(사면체)', 'Cl4': 'Td(사면체)'}


def shape_of(name):
    for suf, s in SHAPE.items():
        if name.endswith(suf):
            return s
    return '-'


def flag_of(name):
    v = IMAG.get(name)
    if v is None:
        return 'freq 합격'
    return '허수 %d cm-1 %s' % (v, '(50 이하: 물 회전, 무시 가능)' if v <= 50
                                else '** 안장점 — 구조 인용 시 명시 **')


def read_xyz(path):
    """(원자수, 좌표줄 리스트). 원본 주석줄은 버리고 새로 붙인다."""
    with open(path, encoding='utf-8') as f:
        lines = [l.rstrip('\n') for l in f if l.strip()]
    n = int(lines[0].split()[0])
    return n, lines[2:2 + n]


def build(pick, tag, out_name):
    """pick(sub, name) -> 파일경로 또는 None. 찾은 것만 묶는다."""
    frames, found = [], []
    for sub, name, q, mult in LADDER + EXTRA:
        src = pick(sub, name)
        if not src or not os.path.exists(src):
            continue
        n, body = read_xyz(src)
        title = '%s  charge=%+d  mult=%d  %s  |  %s  |  %s' % (
            name, q, mult, shape_of(name), flag_of(name), tag)
        frames.append('%d\n%s\n%s\n' % (n, title, '\n'.join(body)))
        found.append((name, n, q, mult))

    out = os.path.join(GEO, out_name)
    with open(out, 'w', encoding='utf-8', newline='\n') as f:
        f.write(''.join(frames))

    # 자체검증: 원자수 선언이 실제 줄수와 어긋나면 Avogadro 는 에러 없이
    # 조용히 마지막 구조만 보여준다. 그게 제일 위험하므로 여기서 잡는다.
    with open(out, encoding='utf-8') as f:
        lines = [l.rstrip('\n') for l in f]
    i = cnt = 0
    while i < len(lines) and lines[i].strip():
        n = int(lines[i].split()[0])
        assert len([l for l in lines[i + 2:i + 2 + n] if l.strip()]) == n, \
            '프레임 %d 원자수 불일치' % cnt
        i += 2 + n
        cnt += 1
    assert cnt == len(frames), '프레임 수 불일치'
    return out, found, cnt



def write_ms():
    """Materials Studio 용 개별 파일.  MS 의 XYZ 파서는 세 가지를 못 읽는다:
      1) 다중프레임 (한 파일에 구조 하나만)
      2) 주석줄의 비 ASCII 문자 -> 'cannot read number of atoms' 로 나온다
      3) LF 줄끝 (Windows 빌드라 CRLF 를 기대)
    그래서 ASCII 주석 + CRLF + 구조당 파일 하나로 따로 뽑는다.
    """
    d = os.path.join(GEO, 'materials_studio')
    os.makedirs(d, exist_ok=True)
    n_ok = 0
    for sub, name, q, mult in LADDER + EXTRA:
        src = None
        for fn in ('%s_opt.xyz' % name, '%s_start.xyz' % name):
            cand = os.path.join(OPT, fn)
            if os.path.exists(cand):
                src = cand
                break
        if not src:
            continue
        n, body = read_xyz(src)
        imag = IMAG.get(name)
        title = '%s charge=%+d mult=%d %s %s' % (
            name, q, mult,
            'Oh' if shape_of(name).startswith('Oh') else ('Td' if shape_of(name).startswith('Td') else 'ref'),
            'freq-OK' if imag is None else 'imag=%dcm-1' % imag)
        title = title.encode('ascii', 'ignore').decode('ascii')
        out = os.path.join(d, '%s.xyz' % name)
        # newline='\r\n' 이 핵심. MS 는 LF 만 있으면 첫 줄을 못 끊어 읽고
        # 'cannot read number of atoms' 를 낸다.
        with open(out, 'w', encoding='ascii', newline='\r\n') as f:
            f.write('%d\n%s\n%s\n' % (n, title, '\n'.join(body)))
        n_ok += 1
    print('\n[Materials Studio 용] %d 개 -> %s' % (n_ok, d))
    print('  ASCII 주석 + CRLF + 구조당 파일 1개')
    return d

def main():
    os.makedirs(GEO, exist_ok=True)

    # 1) 최적화 완료본 (서버에서 받은 것)
    def pick_opt(sub, name):
        for fn in ('%s_opt.xyz' % name, '%s_start.xyz' % name):   # Cl- 는 opt 없음
            p = os.path.join(OPT, fn)
            if os.path.exists(p):
                return p
        return None

    out1, found1, n1 = build(pick_opt, 'r2SCAN-3c/CPCM(SMD water) 최적화 완료',
                             'all_species_OPTIMIZED.xyz')
    print('[최적화 완료본]  %d 프레임' % n1)
    for name, n, q, mult in found1:
        print('  %-12s %2d atoms  q=%+d  mult=%d  %s  %s'
              % (name, n, q, mult, shape_of(name), flag_of(name)))

    # 2) 최적화 전 초기 구조 (비교용)
    out2, _, n2 = build(lambda sub, name: os.path.join(HERE, sub, '%s_start.xyz' % name),
                        '최적화 전 초기 구조', 'all_species_START.xyz')

    # 개별 파일도 한 폴더에 (Avogadro 에서 하나씩 열고 싶을 때)
    for sub, name, q, mult in LADDER + EXTRA:
        src = pick_opt(sub, name)
        if src:
            shutil.copyfile(src, os.path.join(GEO, '%s.xyz' % name))

    print('\nwrote %s  (%d 프레임)' % (out1, n1))
    print('wrote %s  (%d 프레임)' % (out2, n2))
    write_ms()
    print('self-check OK: 두 파일 모두 프레임·원자수 일치')


if __name__ == '__main__':
    main()
