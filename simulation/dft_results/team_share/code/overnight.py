"""밤새 혼자 돌면서 대칭 안장점 문제를 끝까지 해결하고 후처리까지 마치는 감독 스크립트.

배경 (2026-09-26):
  Ni_aq6(기준물질)와 NiCl2_aq4 가 고대칭 초기구조 탓에 안장점에 수렴해
  허수진동수가 각각 10개(430 cm-1), 8개(476 cm-1) 나왔다. ORCA 열역학은 허수모드를
  분배함수에서 빼므로 기준물질의 엔트로피가 틀어지고 그 오차가 모든 dG 에 실린다.
  build_redo.py 가 물 배향을 무작위로 틀어 대칭을 깨고 재계산한다.
  한 번에 될 보장이 없으므로 여기서 라운드를 돌린다.

하는 일:
  1) redo 잡이 끝날 때까지 기다린다
  2) 허수진동수를 검사한다 (|imag| <= IMAG_TOL 이면 합격)
  3) 불합격이면 다른 시드로 대칭을 다시 깨서 재제출 (최대 MAX_ROUNDS 회)
  4) 합격한 종만 본 폴더에 설치한다 (원본은 saddle_backup/ 으로 백업)
  5) speciation.py + extract_dd.py 실행
  6) MORNING_REPORT.txt 한 장으로 정리

설계 원칙 (무인 실행):
  - 어떤 단계가 실패해도 예외를 밖으로 던지지 않는다. 기록하고 다음으로 간다.
  - 무한 대기 금지. 전체 한도(DEADLINE_H)와 라운드별 한도를 둔다.
  - 개선이 없으면 원본을 유지한다. 나쁜 결과로 좋은 결과를 덮지 않는다.
  - 마지막에 무슨 일이 있어도 MORNING_REPORT.txt 를 쓴다.

실행:
    sbatch run_overnight.sh          # 권장 (세션 끊겨도 살아남음)
    python overnight.py              # 직접 실행도 가능
"""
import os, re, shutil, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
REDO = os.path.join(HERE, 'redo')
BACKUP = os.path.join(HERE, 'saddle_backup')
LOG = os.path.join(HERE, 'overnight.log')
REPORT = os.path.join(HERE, 'MORNING_REPORT.txt')

SPECIES = ['Ni_aq6', 'NiCl2_aq4']
MAX_ROUNDS = 3
IMAG_TOL = 50.0        # cm-1. 이보다 작은 허수는 물 회전 잔재라 무시 가능
POLL = 120             # 초
ROUND_LIMIT_H = 4.0    # 한 라운드가 이보다 길어지면 포기하고 다음으로
DEADLINE_H = 9.0       # 전체 한도. 일요일 아침 전에 반드시 끝나도록
JOBNAME = 'ni_redo'

T0 = time.time()


def log(msg):
    line = '[%s] %s' % (time.strftime('%m-%d %H:%M:%S'), msg)
    print(line, flush=True)
    try:
        with open(LOG, 'a', encoding='utf-8') as f:
            f.write(line + '\n')
    except Exception:
        pass


def sh(cmd, timeout=600):
    """셸 명령. 절대 예외를 던지지 않는다. (rc, 출력) 반환."""
    try:
        p = subprocess.run(cmd, shell=True, cwd=HERE, timeout=timeout,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        return p.returncode, p.stdout.decode('utf-8', 'replace')
    except subprocess.TimeoutExpired:
        return 124, '(timeout %ds)' % timeout
    except Exception as e:
        return 1, '(예외: %s)' % e


def elapsed_h():
    return (time.time() - T0) / 3600.0


def imag_of(path):
    """(허수 개수, 최대 크기 cm-1). 파일이 없거나 못 읽으면 (None, 0)."""
    try:
        txt = open(path, encoding='utf-8', errors='ignore').read()
    except Exception:
        return None, 0.0
    m = re.search(r'VIBRATIONAL FREQUENCIES(.*?)(NORMAL MODES|$)', txt, re.S)
    if not m:
        return None, 0.0
    v = [float(x) for x in re.findall(r':\s+(-?\d+\.\d+)\s*cm\*\*-1', m.group(1))]
    neg = [x for x in v if x < -1.0]
    return len(neg), (abs(min(neg)) if neg else 0.0)


def finished_ok(path):
    try:
        return 'ORCA TERMINATED NORMALLY' in open(path, encoding='utf-8',
                                                  errors='ignore').read()
    except Exception:
        return False


def jobs_alive():
    rc, out = sh('squeue -u "$USER" -h -n %s -o "%%i"' % JOBNAME, timeout=120)
    if rc != 0:
        log('squeue 실패(rc=%d) -> 잡이 살아있다고 가정하고 계속 대기' % rc)
        return True
    return bool(out.strip())


def wait_for_jobs(limit_h):
    """redo 잡이 끝날 때까지. 한도를 넘으면 False."""
    t = time.time()
    while jobs_alive():
        if (time.time() - t) / 3600.0 > limit_h:
            log('라운드 한도 %.1fh 초과 -> 대기 중단' % limit_h)
            return False
        if elapsed_h() > DEADLINE_H:
            log('전체 한도 %.1fh 초과 -> 대기 중단' % DEADLINE_H)
            return False
        time.sleep(POLL)
    return True


def check():
    """종별 (합격여부, 개수, 최대크기). freq 가 없거나 비정상 종료면 불합격."""
    out = {}
    for n in SPECIES:
        f = os.path.join(REDO, '%s_freq.out' % n)
        if not finished_ok(f):
            out[n] = (False, None, 0.0)
            continue
        c, mx = imag_of(f)
        out[n] = (c is not None and mx <= IMAG_TOL, c, mx)
    return out


def resubmit(seed):
    """다른 시드로 대칭을 다시 깨고 재제출."""
    # 이전 라운드 흔적 제거: .done 플래그와 _geo.xyz 가 남아 있으면 새 초기구조가 무시된다
    for n in SPECIES:
        for pat in ('%s.xtb.done', '%s.opt.done', '%s.freq.done', '%s_geo.xyz'):
            try:
                os.remove(os.path.join(REDO, pat % n))
            except OSError:
                pass
    rc, out = sh('python build_redo.py --seed %d' % seed)
    log('build_redo(seed=%d) rc=%d\n%s' % (seed, rc, out.strip()))
    if rc != 0:
        return False
    rc, out = sh('cd redo && sbatch run_redo.sh')
    log('sbatch rc=%d %s' % (rc, out.strip()))
    return rc == 0


def install(good):
    """합격한 종을 본 폴더에 설치. 원본은 백업. 나쁜 결과로 덮지 않는다."""
    os.makedirs(BACKUP, exist_ok=True)
    done = []
    for n, (ok, c, mx) in good.items():
        if not ok:
            log('%s 불합격 -> 원본 유지 (허수 %s개, 최대 %.0f)' % (n, c, mx))
            continue
        moved = 0
        for ext in ('_freq.out', '_opt.xyz', '_opt.out'):
            src = os.path.join(REDO, n + ext)
            dst = os.path.join(HERE, n + ext)
            if not os.path.exists(src):
                continue
            try:
                if os.path.exists(dst):
                    shutil.copy2(dst, os.path.join(BACKUP, n + ext))
                shutil.copy2(src, dst)
                moved += 1
            except Exception as e:
                log('%s%s 설치 실패: %s' % (n, ext, e))
        log('%s 설치 완료 (파일 %d개, 허수 %s개/최대 %.0f cm-1)' % (n, moved, c, mx))
        done.append(n)
    return done


def write_report(rounds, final, installed, post):
    lines = []
    A = lines.append
    A('=' * 68)
    A('일요일 아침 보고서 — %s' % time.strftime('%Y-%m-%d %H:%M'))
    A('=' * 68)
    A('')
    A('[1] 대칭 안장점 재계산')
    A('  라운드 %d회 수행, 총 소요 %.1f시간' % (rounds, elapsed_h()))
    A('  합격 기준: 허수진동수 최대 크기 <= %.0f cm-1' % IMAG_TOL)
    A('')
    A('  %-12s %-8s %-8s %-8s' % ('종', '합격', '허수개수', '최대cm-1'))
    for n, (ok, c, mx) in final.items():
        A('  %-12s %-8s %-8s %-8.0f' % (n, 'YES' if ok else 'NO',
                                        '?' if c is None else c, mx))
    A('')
    A('  본 폴더에 설치된 종: %s' % (', '.join(installed) if installed else '없음'))
    if len(installed) < len(SPECIES):
        A('  ** 설치 안 된 종이 있습니다. 그 종의 dG 는 여전히 안장점 기반이므로')
        A('     보고서에 인용하지 말고, 해당 종만 다시 돌리세요.')
        A('     원본은 saddle_backup/ 에 보관되어 있습니다.')
    A('')
    A('[2] 후처리')
    for name, (rc, out) in post.items():
        A('  --- %s (rc=%d) ---' % (name, rc))
        for l in out.strip().splitlines():
            A('    ' + l)
        A('')
    A('[3] 다음에 볼 것')
    A('  - speciation.csv / guideline_map.png : 조성-온도 지도')
    A('  - dd_bands.csv / dd_spectra.png      : d-d 밴드와 eps')
    A('  - overnight.log                      : 밤새 진행 로그 전체')
    A('  - saddle_backup/                     : 교체 전 원본(안장점) 결과')
    A('')
    A('  확인 포인트: shift 절대값이 50 kJ/mol 아래로 내려왔는지.')
    A('  아직 크면 dG 자체에 남은 문제가 있다는 뜻입니다.')
    A('')
    A('[4] 이 실행이 건드리지 않은 것 (알고 있어야 할 한계)')
    A('  - d-d 계산(dd_bands.csv)은 교체 전 기하로 돌린 결과 그대로입니다.')
    A('    d-d 에너지는 Ni-리간드 거리가 지배하고 물 배향에는 둔감하므로 큰 영향은')
    A('    없지만, 엄밀히 하려면 Ni_aq6 의 dd 를 새 기하로 다시 돌려야 합니다.')
    A('  - NiCl4 계산 파장(479-532 nm)은 문헌(~660-700 nm)보다 에너지가 약 35% 높습니다.')
    A('    CAS(8,5) 의 double-shell 누락 탓입니다. 절대 파장은 인용하지 마세요.')
    A('  - 계산 eps(Td) 합계는 약 29 M-1cm-1 로 문헌값(150-200)보다 작습니다.')
    A('    CASSCF 가 d-d 진동자세기를 과소평가하는 알려진 문제입니다.')
    A('    검출한계 논증에 어느 값을 쓸지 정하고, 그 선택을 보고서에 밝히세요.')
    txt = '\n'.join(lines) + '\n'
    try:
        open(REPORT, 'w', encoding='utf-8', newline='\n').write(txt)
        log('보고서 작성: %s' % REPORT)
    except Exception as e:
        log('보고서 작성 실패: %s' % e)
    print(txt)


def main():
    log('=' * 60)
    log('감독 시작. 대상=%s, 최대 %d라운드, 전체한도 %.1fh'
        % (','.join(SPECIES), MAX_ROUNDS, DEADLINE_H))

    final = {n: (False, None, 0.0) for n in SPECIES}
    rounds = 0
    for rnd in range(1, MAX_ROUNDS + 1):
        rounds = rnd
        log('--- 라운드 %d: 잡 대기 ---' % rnd)
        wait_for_jobs(ROUND_LIMIT_H)
        final = check()
        for n, (ok, c, mx) in final.items():
            log('  %s: %s (허수 %s개, 최대 %.0f cm-1)'
                % (n, '합격' if ok else '불합격', c, mx))
        if all(v[0] for v in final.values()):
            log('전 종 합격. 라운드 종료.')
            break
        if rnd == MAX_ROUNDS:
            log('최대 라운드 도달. 합격한 종만 설치하고 진행.')
            break
        if elapsed_h() > DEADLINE_H - ROUND_LIMIT_H:
            log('남은 시간이 부족해 추가 라운드를 시작하지 않음.')
            break
        seed = 20260926 + rnd * 977
        log('--- 라운드 %d 재제출 (seed=%d) ---' % (rnd + 1, seed))
        if not resubmit(seed):
            log('재제출 실패 -> 라운드 중단')
            break

    installed = install(final)

    log('--- 후처리 ---')
    post = {}
    for name, cmd in (('speciation.py', 'python speciation.py'),
                      ('extract_dd.py', 'python extract_dd.py')):
        rc, out = sh(cmd, timeout=1800)
        post[name] = (rc, out)
        log('%s rc=%d (%d줄)' % (name, rc, len(out.splitlines())))

    write_report(rounds, final, installed, post)
    log('감독 종료. 총 %.1f시간.' % elapsed_h())


if __name__ == '__main__':
    try:
        main()
    except Exception as e:               # 무슨 일이 있어도 보고서는 남긴다
        import traceback
        log('치명적 예외: %s' % e)
        log(traceback.format_exc())
        try:
            open(REPORT, 'a', encoding='utf-8').write(
                '\n[감독 스크립트가 예외로 종료됨]\n%s\n'
                'overnight.log 를 확인하세요.\n' % e)
        except Exception:
            pass
        sys.exit(0)                      # Slurm 에 실패로 남기지 않는다
