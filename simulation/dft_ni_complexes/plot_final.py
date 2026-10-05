"""포스터용 최종 그림 — 예측(지도)과 실측(바이알 사진)을 한 장에.  한/영 + PNG/SVG.

설계 의도:
  이전 그림들은 예쁘지만 **읽어야** 알 수 있었다. 발표 중 내 파트가 1분이면 치명적이다.
  그래서 (1) 제목을 주제가 아니라 결론으로, (2) 전문용어를 평이한 말로,
  (3) 바이알 사진을 조성 축 제자리에 붙여 예측·실측을 같은 x 에서 맞춰 보이게,
  (4) 관행 몰비를 굵은 빨간 선으로 — 그게 우리 표적이므로.

라벨 위치는 눈대중이 아니라 각 영역의 실제 폭을 계산해 잡았다
(teal 은 x>=0.70 에서 약 33 °C 폭 일정, yellow 는 우상단에서 최대 57 °C).

SVG 출력:
  svg.fonttype='none' 로 저장해 **텍스트가 텍스트로 남는다** -> PowerPoint 에서
  SVG 를 넣고 [도형으로 변환] 하면 글자·도형을 따로 편집할 수 있다.
  단, 그 PC 에 같은 한글 폰트가 있어야 하므로 Windows 기본 탑재인 Malgun Gothic 을 쓴다.
  사진은 SVG 안에 base64 PNG 로 박히므로 그대로 보인다.

실행:
    python dft_ni_complexes/plot_final.py          # 한글 + 영문, PNG + SVG 전부
입력: dft_results/window.csv (서버 계산), photos/vials_260922.png
"""
import csv, os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from PIL import Image

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV = os.path.join(ROOT, 'dft_results', 'window.csv')
PHOTO = os.path.join(ROOT, 'photos', 'vials_260922.png')
OUTDIR = os.path.join(ROOT, 'dft_results', 'figure_styles')

VIALS = [('1:19', 0.05, 1900), ('3:17', 0.15, 1050), ('1:4', 0.20, 1500),
         ('3:1', 0.75, 2350), ('19:1', 0.95, 2800)]
CROP_Y, CROP_W = (1900, 2120), 95        # 좁게 잘라야 0.05 간격에서 안 겹친다

INK, INK2 = '#0b0b0b', '#55544f'
RED, GREEN = '#ff4d4d', '#0f7a5a'
VIR = plt.get_cmap('viridis')
FILL = [VIR(0.06), VIR(0.52), VIR(0.90)]

TXT = {
    'ko': dict(
        font='Malgun Gothic',
        # 주장 범위(9/27): 화학종까지만. '분리·선택성·아무 일도 없다' 는 성능 언어라
        # 청중이 침출 성능으로 듣는다. 전부 '화학종' 표현으로 바꿨다 (2026-09-30).
        title='조성비가 어느 금속을 음이온으로 만들지 결정한다\n'
              '관행 몰비 1:2 는 세 금속이 모두 같은 화학종으로 남는 자리에 있다',
        sub='DFT(r2SCAN-3c) + COSMO-RS, 실험과 동일한 물 함량에서 계산.   '
            'Ni·Co 는 각각 자신의 문헌 실측값에 고정(anchor).',
        dead='세 금속 모두 양이온으로 남음\n화학종으로는 서로 구분되지 않음',
        # mathtext 의 아래/위 첨자가 닫는 괄호와 겹쳐 지저분했다. 발표에서도
        # "사면체 코발트" 라고 읽으므로 화학식 대신 말로 쓴다.
        co='Co → 음이온 (사면체 CoCl₄²⁻)\nNi → 여전히 양이온',
        mn='Mn 도 합류\n(Mn 미보정)',
        conv='관행 몰비 1:2 는\n세 금속이 같은 화학종인\n영역 한복판',
        exp='우리 실험 · 60 °C',
        # 50 % 문턱으로 보면 0.75 는 24 % 라 '음이온이 된다' 가 아니다.
        # 실제로 맞힌 것은 '사면체가 생기기 시작하는 문턱의 위치' 이므로 그렇게 쓴다.
        hit='사면체 Co 를 예측한 두 조성이\n실제로 색이 난 두 조성과 일치\n(계산 24 % · 96 %, 나머지 3개는 0)',
        xlab='x(ChCl)          ←  시트르산 많음          염(ChCl) 많음  →',
        ylab='온도  /  °C'),
    'en': dict(
        font='DejaVu Sans',
        title='Composition decides which metal becomes an anion —\n'
              'and the conventional 1:2 ratio sits where all three stay the same',
        sub='DFT (r2SCAN-3c) + COSMO-RS at the experimental water content.   '
            'Ni and Co each anchored to their own published measurement.',
        dead='all three metals stay CATIONIC\nspeciation cannot tell them apart',
        co='Co $\\rightarrow$ ANION (CoCl$_4^{2-}$)\nNi $\\rightarrow$ still a cation',
        mn='Mn joins Co\n(Mn uncalibrated)',
        conv='the conventional 1:2 ratio sits\nwhere all three share one species',
        exp='our experiments · 60 °C',
        hit='the two we predicted to be anionic\nare exactly the two that turned blue',
        xlab='x(ChCl)          ←  more citric acid          more salt (ChCl)  →',
        ylab='temperature  /  °C'),
}


def load():
    xs, co, mn = [], [], []
    for r in csv.DictReader(open(CSV, encoding='utf-8')):
        xs.append(float(r['x_ChCl']))
        for key, dst in (('T_Co_C', co), ('T_Mn_C', mn)):
            v = r[key]
            dst.append(np.nan if v.startswith('>') else
                       (-99.0 if v.startswith('<') else float(v)))
    return np.array(xs), np.array(co), np.array(mn)


def draw(lang):
    t = TXT[lang]
    xs, co, mn = load()
    ymin, ymax = 25.0, 125.0
    co = np.where(np.isnan(co), ymax, np.where(co < 0, ymin, co))
    mn = np.maximum(np.where(np.isnan(mn), ymax, np.where(mn < 0, ymin, mn)), co)

    plt.rcParams.update({
        'font.family': t['font'], 'font.size': 11,
        'axes.unicode_minus': False,          # 한글 폰트에서 마이너스가 깨지는 것 방지
        'text.color': INK, 'axes.labelcolor': INK,
        'xtick.color': INK2, 'ytick.color': INK2,
        'svg.fonttype': 'none',               # SVG 에서 글자를 글자로 유지 (PPT 편집용)
    })
    fig = plt.figure(figsize=(10.6, 8.0))
    ax = fig.add_axes([0.085, 0.315, 0.885, 0.495])
    axp = fig.add_axes([0.085, 0.125, 0.885, 0.155], sharex=ax)

    pad = 0.045
    x0, x1 = xs.min() - pad, xs.max() + pad

    # lw=0 : 폭이 0인 구간에서도 가장자리 선이 그려져 상단에 얇은 띠가 생겼다
    ax.fill_between(xs, ymin, co, color=FILL[0], lw=0)
    ax.fill_between(xs, co, mn, color=FILL[1], lw=0)
    ax.fill_between(xs, mn, ymax, color=FILL[2], lw=0)
    ax.axvspan(x0, xs[0], color=FILL[0], lw=0)
    for lo, hi, c in ((ymin, co[-1], FILL[0]), (co[-1], mn[-1], FILL[1]),
                      (mn[-1], ymax, FILL[2])):
        ax.fill_between([xs[-1], x1], lo, hi, color=c, lw=0)

    ok = co < ymax - 0.1
    ax.plot(xs[ok], co[ok], '-', lw=3.0, color='w', solid_capstyle='round')
    okm = (mn < ymax - 0.1) & (mn > co + 0.1)
    ax.plot(xs[okm], mn[okm], ls=(0, (6, 3)), lw=2.6, color='w',
            dash_capstyle='round')

    ax.text(0.075, 100, t['dead'], color='w', fontsize=13, fontweight='bold',
            ha='left', va='center', linespacing=1.5)
    # 띠 폭(약 33 °C)과 대각 방향을 감안해 2줄 — 3줄은 노란 영역까지 넘쳤다
    ax.text(0.858, 70, t['co'], color='w', fontsize=11, fontweight='bold',
            ha='center', va='center', linespacing=1.4)
    ax.text(0.875, 106, t['mn'], color='#14203a', fontsize=10, fontweight='bold',
            ha='center', va='center', linespacing=1.35)

    ax.axvline(0.333, color=RED, lw=3.2, zorder=5)
    ax.annotate('', xy=(0.345, 42), xytext=(0.60, 42), zorder=6,
                arrowprops=dict(arrowstyle='<|-', color=RED, lw=2.4,
                                mutation_scale=20))
    ax.text(0.615, 42, t['conv'], color=RED, fontsize=12, fontweight='bold',
            va='center', zorder=6, linespacing=1.35)

    ax.axhline(60, color='w', ls=(0, (1, 2.6)), lw=1.6)
    ax.text(x0 + 0.012, 63, t['exp'], color='w', fontsize=10)

    ax.set_xlim(x0, x1); ax.set_ylim(ymin, ymax)
    ax.set_ylabel(t['ylab'], labelpad=8)
    ax.tick_params(labelbottom=False, length=3.5)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)

    # --- 아래 띠: 실제 바이알 사진을 같은 조성 위치에 ---------------------
    src = Image.open(PHOTO).convert('RGB')
    for lab, x, cx in VIALS:
        im = src.crop((cx - CROP_W, CROP_Y[0], cx + CROP_W, CROP_Y[1]))
        im = im.resize((76, 88), Image.LANCZOS)
        hit = x in (0.75, 0.95)          # 모델이 '음이온' 이라 예측한 두 개
        axp.add_artist(AnnotationBbox(
            OffsetImage(np.asarray(im), zoom=0.95), (x, 0.62), frameon=True,
            pad=0.12, bboxprops=dict(edgecolor=GREEN if hit else '#c9c8c3',
                                     lw=3.2 if hit else 1.0)))
        axp.text(x, 0.015, lab, ha='center', va='bottom', fontsize=10.5,
                 fontweight='bold')
    axp.set_ylim(0, 1); axp.set_yticks([])
    axp.set_xlabel(t['xlab'], labelpad=10, fontsize=12)
    for s in ('top', 'right', 'left'):
        axp.spines[s].set_visible(False)
    axp.axvline(0.333, color=RED, lw=3.2, alpha=0.5, zorder=0)
    # 검증 문구는 사진 띠의 빈 구간(x 0.25~0.70)에
    axp.text(0.475, 0.62, t['hit'], ha='center', va='center', fontsize=11.5,
             fontweight='bold', color=GREEN, linespacing=1.45, zorder=6)

    fig.text(0.085, 0.975, t['title'], fontsize=17.5, fontweight='bold',
             ha='left', va='top', linespacing=1.32)
    fig.text(0.085, 0.862, t['sub'], fontsize=10, color=INK2, ha='left', va='top')

    for ext, dpi in (('png', 200), ('svg', None)):
        p = os.path.join(OUTDIR, 'D_poster_%s.%s' % (lang, ext))
        fig.savefig(p, facecolor='white', bbox_inches='tight',
                    **({'dpi': dpi} if dpi else {}))
        print('wrote', p)
    plt.close(fig)


if __name__ == '__main__':
    os.makedirs(OUTDIR, exist_ok=True)
    for lg in ('ko', 'en'):
        draw(lg)
