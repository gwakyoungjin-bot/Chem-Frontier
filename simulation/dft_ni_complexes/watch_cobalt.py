"""Co·Mn 계산이 끝나면 개입 없이 검증 -> 금속별 재시도 -> 후처리까지 이어간다.

왜 금속별로 독립인가:
  Co 가 실패해도 Mn 은 쓸 수 있어야 하고, 그 반대도 마찬가지다.
  한 금속이 막혔다고 다른 금속 결과까지 버리면 아침에 볼 게 없어진다.
  => 검사·재시도를 금속마다 따로 돌리고, 후처리는 '성공한 금속들' 로 진행한다.
     cobalt_window.py 가 없는 금속을 알아서 건너뛰므로 이게 성립한다.

흐름:
  1) co_dft 잡(Co+Mn 공용 이름)이 전부 끝날 때까지 대기
  2) 금속별로 허수진동수 검사 (|imag| <= 50 cm-1 합격)
  3) 불합격 금속만 다른 시드로 다시 만들어 재제출 (최대 MAX_ROUNDS)
  4) 후처리: speciation.py --water, cobalt_window.py --water
  5) CO_REPORT.txt 한 장

안전장치는 overnight.py 와 동일하며 그 헬퍼를 그대로 재사용한다:
  예외 비전파 / 무한대기 금지 / 나쁜 결과로 덮지 않음 / 무조건 보고서 작성.

실행:
    sbatch run_watch_cobalt.sh
"""
import os, sys, time

import overnight as O          # sh / log / imag_of / finished_ok / wait_for_jobs 재사용

HERE = os.path.dirname(os.path.abspath(__file__))

# 금속 -> (결과폴더, 종 목록, 재빌드 인자, 제출 명령)
METALS = {
    'Co': dict(dir='cobalt', build='--metal Co',
               run='cd cobalt && sbatch run_cobalt.sh',
               species=['CoCl4', 'Co_aq6', 'CoCl2_aq4', 'CoCl3_aq1', 'CoCl1_aq5']),
    'Mn': dict(dir='manganese', build='--metal Mn',
               run='cd manganese && sbatch run_manganese.sh',
               species=['MnCl4', 'Mn_aq6', 'MnCl2_aq4', 'MnCl3_aq1', 'MnCl1_aq5']),
}

MAX_ROUNDS = 3
ROUND_LIMIT_H = 8.0
DEADLINE_H = 20.0

O.JOBNAME = 'co_dft'                 # Co·Mn 둘 다 이 이름으로 제출돼 있다
O.LOG = os.path.join(HERE, 'watch_cobalt.log')
O.DEADLINE_H = DEADLINE_H
O.T0 = time.time()

REPORT = os.path.join(HERE, 'CO_REPORT.txt')


def check(metal):
    """{종: (합격, 허수개수, 최대크기)}"""
    m = METALS[metal]
    out = {}
    for n in m['species']:
        f = os.path.join(HERE, m['dir'], '%s_freq.out' % n)
        if not O.finished_ok(f):
            out[n] = (False, None, 0.0)
            continue
        c, mx = O.imag_of(f)
        out[n] = (c is not None and mx <= O.IMAG_TOL, c, mx)
    return out


def ok(status):
    return all(v[0] for v in status.values())


def resubmit(metal, seed):
    """그 금속만 다른 시드로 다시 만들어 제출."""
    m = METALS[metal]
    d = os.path.join(HERE, m['dir'])
    for n in m['species']:
        for pat in ('%s.xtb.done', '%s.opt.done', '%s.freq.done', '%s.dd.done',
                    '%s_geo.xyz'):
            try:
                os.remove(os.path.join(d, pat % n))
            except OSError:
                pass
    rc, out = O.sh('python build_cobalt.py %s --seed %d' % (m['build'], seed))
    O.log('[%s] build rc=%d\n%s' % (metal, rc, out.strip()))
    if rc != 0:
        return False
    rc, out = O.sh(m['run'])
    O.log('[%s] sbatch rc=%d %s' % (metal, rc, out.strip()))
    return rc == 0


def write_report(final, rounds, post, usable):
    A = []
    A.append('=' * 70)
    A.append('Co · Mn 계산 보고서 — %s' % time.strftime('%Y-%m-%d %H:%M'))
    A.append('=' * 70)
    A.append('')
    A.append('[0] 한 줄 요약')
    A.append('  바로 쓸 수 있는 금속: %s' % (', '.join(usable) if usable else '없음'))
    A.append('  그림: window_map.png (분리 영역도) / guideline_map.png (Ni 단독)')
    A.append('  표  : window.csv (금속별 전환선) / speciation.csv (Ni 화학종 분율)')
    A.append('')
    A.append('[1] 계산 상태 (합격 = 허수진동수 최대 |값| <= %.0f cm-1)' % O.IMAG_TOL)
    for metal in METALS:
        st = final.get(metal, {})
        A.append('')
        A.append('  == %s ==  라운드 %d회' % (metal, rounds.get(metal, 0)))
        A.append('  %-14s %-8s %-9s %-8s' % ('종', '합격', '허수개수', '최대cm-1'))
        for n in METALS[metal]['species']:
            o, c, mx = st.get(n, (False, None, 0.0))
            A.append('  %-14s %-8s %-9s %-8.0f'
                     % (n, 'YES' if o else 'NO', '?' if c is None else c, mx))
        bad = [n for n in METALS[metal]['species'] if not st.get(n, (False,))[0]]
        if bad:
            A.append('  ** 불합격: %s' % ', '.join(bad))
            A.append('     -> 이 금속의 dG 는 인용 전에 확인하세요.')
    A.append('')
    A.append('[2] 후처리 출력')
    for name, (rc, out) in post.items():
        A.append('  --- %s (rc=%d) ---' % (name, rc))
        for l in out.strip().splitlines():
            A.append('    ' + l)
        A.append('')
    A.append('[3] 해석 시 지켜야 할 선 (2026-09-27 확정 주장 범위)')
    A.append('  - 조성비 -> 활동도 -> 화학종 까지만 주장한다.')
    A.append('  - 화학종 -> 침출효율 은 제언(future work). 효율 수치를 말하지 않는다.')
    A.append('')
    A.append('[4] 반드시 같이 말할 한계')
    A.append('  - Ni 앵커로 구한 보정값(shift, alpha)을 Co·Mn 에 전용했다.')
    A.append('    용매 오차가 금속에 무관하다는 가정이다.')
    A.append('  - a_Cl 은 중성 이온쌍 근사의 대리변수다. 절대값 인용 금지.')
    A.append('  - Mn(II) 는 d5 고스핀이라 d-d 가 전부 스핀금지 -> 분광 검증 불가, dG 만 유효.')
    A.append('  - Ni<Mn<Co 같은 착물 형성 순서를 문헌으로 확정하지 못했다.')
    A.append('    즉 이 결과는 "아는 답의 확인" 이 아니라 **예측**이다. 원문 대조가 남아 있다.')
    txt = '\n'.join(A) + '\n'
    try:
        open(REPORT, 'w', encoding='utf-8', newline='\n').write(txt)
        O.log('보고서: %s' % REPORT)
    except Exception as e:
        O.log('보고서 작성 실패: %s' % e)
    print(txt)


def main():
    O.log('=' * 60)
    O.log('Co·Mn 감독 시작 (금속별 독립 재시도)')

    final = {m: {} for m in METALS}
    rounds = {m: 0 for m in METALS}
    pending = list(METALS)

    for rnd in range(1, MAX_ROUNDS + 1):
        if not pending:
            break
        O.log('--- 라운드 %d: co_dft 대기 (대상 %s) ---' % (rnd, ', '.join(pending)))
        O.wait_for_jobs(ROUND_LIMIT_H)

        still = []
        for metal in pending:
            rounds[metal] = rnd
            st = check(metal)
            final[metal] = st
            for n, (o, c, mx) in st.items():
                O.log('  [%s] %s: %s (허수 %s개, 최대 %.0f)'
                      % (metal, n, '합격' if o else '불합격', c, mx))
            if ok(st):
                O.log('  [%s] 전 종 합격.' % metal)
            else:
                still.append(metal)
        pending = still
        if not pending or rnd == MAX_ROUNDS:
            break
        if O.elapsed_h() > DEADLINE_H - ROUND_LIMIT_H:
            O.log('남은 시간 부족 -> 추가 라운드 없음.')
            break
        for metal in list(pending):
            seed = 20260927 + rnd * 613 + (0 if metal == 'Co' else 101)
            O.log('--- [%s] 라운드 %d 재제출 (seed=%d) ---' % (metal, rnd + 1, seed))
            if not resubmit(metal, seed):
                O.log('  [%s] 재제출 실패 -> 이 금속은 포기' % metal)
                pending.remove(metal)

    usable = [m for m in METALS if final.get(m) and ok(final[m])]
    O.log('사용 가능 금속: %s' % (', '.join(usable) or '없음'))

    O.log('--- 후처리 ---')
    post = {}
    for name, cmd in (('speciation.py --water', 'python speciation.py --water'),
                      ('cobalt_window.py --water', 'python cobalt_window.py --water')):
        rc, out = O.sh(cmd, timeout=1800)
        post[name] = (rc, out)
        O.log('%s rc=%d (%d줄)' % (name, rc, len(out.splitlines())))

    write_report(final, rounds, post, usable)
    O.log('감독 종료. 총 %.1f시간.' % O.elapsed_h())


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        import traceback
        O.log('치명적 예외: %s' % e)
        O.log(traceback.format_exc())
        try:
            open(REPORT, 'a', encoding='utf-8').write(
                '\n[감독이 예외로 종료됨]\n%s\nwatch_cobalt.log 확인\n' % e)
        except Exception:
            pass
        sys.exit(0)
