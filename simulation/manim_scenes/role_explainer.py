"""오늘 재정의된 시뮬레이션 파트의 역할·방법·결과를 설명하는 애니메이션.

5장 구성:
  1 문제     관행 몰비는 '언제 녹느냐'의 답이지 '무엇이 생기느냐'의 답이 아니다
  2 방법     DFT -> COSMO-RS -> 문헌 앵커 -> 지도, 4단계
  3 검증1    바이알 10개와 계산 f_Td 가 10/10 일치
  4 검증2    FT-IR C=O 계단도 같은 자리를 가리킨다
  5 결론     세 증거가 x(ChCl) ~ 0.5 에서 만난다

수치는 전부 실제 산출물에서 가져왔다:
  f_Td   dft_results/ftd_60C.csv 및 ftd_scan (60 C, 앵커 보정)
  C=O    ftir_trend.py (13조성 실측)
  보정값 dft_results/shifts.csv

렌더:
  .venv-manim\\Scripts\\python.exe -m manim -qm manim_scenes/role_explainer.py RoleExplainer
"""
from manim import *

KO = "Malgun Gothic"

BG = "#0e1413"
INK = "#e4edeb"
MUTED = "#93a4a1"
TEAL = "#4fc4b8"
WARM = "#e0a458"
RED = "#e06c5a"
DIM = "#2a3a37"

config.background_color = BG

# --- 실측·계산 데이터 (하드코딩하지 않고 여기 한 곳에만) ---------------
COMP = [  # (라벨, x(ChCl), f_Td(Co) at 60C, 관찰색)
    ("1:19", 0.05, 0.0000, "무색"), ("1:9", 0.10, 0.0000, "무색"),
    ("3:17", 0.15, 0.0000, "무색"), ("1:4", 0.20, 0.0000, "무색"),
    ("1:3", 0.25, 0.0000, "무색"), ("2:3", 0.40, 0.0000, "무색"),
    ("3:2", 0.60, 0.0053, "청록"), ("3:1", 0.75, 0.2424, "청록"),
    ("17:3", 0.85, 0.7668, "청록"), ("19:1", 0.95, 0.9614, "청록"),
]
CO_IR = [(0.05, 1713.9), (0.10, 1714.0), (0.15, 1714.4), (0.20, 1714.7),
         (0.25, 1714.8), (0.333, 1715.3), (0.40, 1715.6), (0.50, 1716.0),
         (0.60, 1722.3), (0.667, 1722.9), (0.75, 1723.5), (0.85, 1726.8),
         (0.95, 1729.2)]


def ko(t, size=32, color=INK, weight=NORMAL):
    return Text(t, font=KO, font_size=size, color=color, weight=weight)


def ticks(ax, xs, ys, xfmt="%.1f", yfmt="%.1f", size=16):
    """축 눈금을 Text 로 붙인다. manim 기본 include_numbers 는 LaTeX 를 요구하는데
    이 PC 에 LaTeX 가 없어서 렌더가 죽는다."""
    g = VGroup()
    for v in xs:
        g.add(Text(xfmt % v, font=KO, font_size=size, color=MUTED)
              .next_to(ax.c2p(v, ax.y_range[0]), DOWN, buff=0.14))
    for v in ys:
        g.add(Text(yfmt % v, font=KO, font_size=size, color=MUTED)
              .next_to(ax.c2p(ax.x_range[0], v), LEFT, buff=0.14))
    return g


def title_card(n, head, sub=None):
    g = VGroup()
    num = Text(f"{n}", font=KO, font_size=20, color=TEAL, weight=BOLD)
    h = ko(head, 40, INK, BOLD)
    g.add(num, h)
    if sub:
        g.add(ko(sub, 24, MUTED))
    g.arrange(DOWN, buff=0.28, aligned_edge=LEFT)
    return g


class RoleExplainer(Scene):
    def construct(self):
        self.s1_problem()
        self.s2_method()
        self.s3_vials()
        self.s4_ftir()
        self.s5_conclusion()

    # ---------------------------------------------------------------- 1
    def s1_problem(self):
        t = title_card("01", "관행 몰비는 무엇의 답인가",
                       "공융점은 '언제 녹느냐'를 정한 값이다")
        t.to_edge(UP, buff=0.8).to_edge(LEFT, buff=1.0)
        self.play(FadeIn(t, shift=UP * 0.3), run_time=1.2)

        q1 = ko("공융점  =  가장 낮은 온도에서 액체가 되는 조성", 30, MUTED)
        q2 = ko("화학종  =  녹은 금속이 어떤 모양이 되는가", 30, TEAL)
        vs = ko("서로 다른 질문이다", 34, WARM, BOLD)
        grp = VGroup(q1, q2, vs).arrange(DOWN, buff=0.55)
        grp.next_to(t, DOWN, buff=1.0)

        self.play(Write(q1), run_time=1.0)
        self.play(Write(q2), run_time=1.0)
        self.play(FadeIn(vs, scale=1.1), run_time=0.8)
        self.wait(1.4)
        self.play(FadeOut(VGroup(t, grp)), run_time=0.7)

    # ---------------------------------------------------------------- 2
    def s2_method(self):
        t = title_card("02", "무엇을 계산했나", "네 단계")
        t.to_edge(UP, buff=0.7).to_edge(LEFT, buff=1.0)
        self.play(FadeIn(t, shift=UP * 0.3), run_time=1.0)

        steps = [
            ("DFT", "착물 15종의 자유에너지", "물 6개 ↔ 염화물 4개"),
            ("COSMO-RS", "조성별 염화물 활동도", "농도가 아니라 활동도"),
            ("문헌 앵커", "용매모형 오차 보정", "Ni +65.1  Co +67.7 kJ/mol"),
            ("지도", "사면체 분율 지도", "60 °C 에서의 사면체 분율"),
        ]
        boxes = VGroup()
        for name, what, note in steps:
            box = RoundedRectangle(width=3.0, height=1.65, corner_radius=0.1,
                                   stroke_color=DIM, stroke_width=2,
                                   fill_color="#161e1d", fill_opacity=1)
            lbl = ko(name, 24, TEAL, BOLD)
            w = ko(what, 17, INK)
            n = ko(note, 13, MUTED)
            inner = VGroup(lbl, w, n).arrange(DOWN, buff=0.16)
            boxes.add(VGroup(box, inner.move_to(box.get_center())))
        boxes.arrange(RIGHT, buff=0.42).scale(0.88).next_to(t, DOWN, buff=0.9)
        boxes.set_x(0)

        arrows = VGroup(*[
            Arrow(boxes[i].get_right(), boxes[i + 1].get_left(),
                  buff=0.08, stroke_width=3, color=MUTED,
                  max_tip_length_to_length_ratio=0.28)
            for i in range(3)])

        for i, b in enumerate(boxes):
            self.play(FadeIn(b, shift=UP * 0.2), run_time=0.55)
            if i < 3:
                self.play(GrowArrow(arrows[i]), run_time=0.3)
        self.wait(1.6)
        self.play(FadeOut(VGroup(t, boxes, arrows)), run_time=0.7)

    # ---------------------------------------------------------------- 3
    def s3_vials(self):
        t = title_card("03", "예측과 실물을 맞춰봤다", "바이알 10개, 계산은 사진을 안 보고 나온 값")
        t.to_edge(UP, buff=0.55).to_edge(LEFT, buff=1.0)
        self.play(FadeIn(t, shift=UP * 0.3), run_time=1.0)

        ax = Axes(x_range=[0, 1.0, 0.2], y_range=[0, 1.0, 0.25],
                  x_length=9.0, y_length=2.6,
                  axis_config=dict(color=MUTED, stroke_width=2,
                                   include_tip=False,
                                   font_size=20),
                  ).next_to(t, DOWN, buff=0.55)
        ax.set_x(0.35)
        xl = ko("x(ChCl)   ←  시트르산 많음        염 많음  →", 19, MUTED)
        xl.next_to(ax, DOWN, buff=0.72)
        yl = ko("사면체 Co 분율", 18, MUTED).rotate(PI / 2).next_to(ax, LEFT, buff=0.95)
        tk = ticks(ax, [], [0, 0.5, 1.0])   # x 는 조성 라벨이 대신한다
        self.play(Create(ax), FadeIn(xl), FadeIn(yl), FadeIn(tk), run_time=1.2)

        # 계산 곡선
        curve = ax.plot_line_graph(
            [c[1] for c in COMP], [c[2] for c in COMP],
            line_color=TEAL, stroke_width=4, add_vertex_dots=False)
        self.play(Create(curve), run_time=1.6)
        lab = ko("계산 곡선", 20, TEAL, BOLD).move_to(ax.c2p(0.27, 0.78))
        self.play(FadeIn(lab), run_time=0.4)

        # 바이알 = 관찰색 원
        dots, texts = VGroup(), VGroup()
        for name, x, f, col in COMP:
            c = "#7fd8cf" if col == "청록" else "#f2f2ee"
            d = Circle(radius=0.135, fill_color=c, fill_opacity=1,
                       stroke_color=INK, stroke_width=1.6)
            d.move_to(ax.c2p(x, f))
            dots.add(d)
            texts.add(ko(name, 13, MUTED).next_to(d, DOWN if f < 0.4 else UP, buff=0.14))
        self.play(LaggedStart(*[FadeIn(d, scale=0.6) for d in dots],
                              lag_ratio=0.12), run_time=2.0)
        self.play(FadeIn(texts), run_time=0.6)

        # 전이 구간
        band = Rectangle(width=ax.c2p(0.60, 0)[0] - ax.c2p(0.40, 0)[0],
                         height=ax.y_length, fill_color=WARM, fill_opacity=0.14,
                         stroke_width=0)
        band.move_to([(ax.c2p(0.40, 0)[0] + ax.c2p(0.60, 0)[0]) / 2,
                      ax.c2p(0, 0.5)[1], 0])
        bl = ko("전이", 18, WARM, BOLD).next_to(band, UP, buff=0.08)
        self.play(FadeIn(band), FadeIn(bl), run_time=0.8)

        hit = ko("무색 6개 → 계산 전부 0          청록 4개 → 계산 전부 0 초과", 24, WARM, BOLD)
        hit.next_to(xl, DOWN, buff=0.34)
        hit.set_x(0)
        self.play(Write(hit), run_time=1.4)
        self.wait(2.0)
        self.play(FadeOut(VGroup(t, ax, tk, xl, yl, curve, lab, dots, texts,
                                 band, bl, hit)), run_time=0.7)

    # ---------------------------------------------------------------- 4
    def s4_ftir(self):
        t = title_card("04", "용매 구조도 같은 자리에서 바뀐다",
                       "FT-IR 카르보닐 C=O 파수, 13조성 실측")
        t.to_edge(UP, buff=0.55).to_edge(LEFT, buff=1.0)
        self.play(FadeIn(t, shift=UP * 0.3), run_time=1.0)

        ax = Axes(x_range=[0, 1.0, 0.2], y_range=[1712, 1731, 5],
                  x_length=9.2, y_length=3.1,
                  axis_config=dict(color=MUTED, stroke_width=2,
                                   include_tip=False, font_size=20),
                  ).next_to(t, DOWN, buff=0.8)
        ax.set_x(0.45)
        xl = ko("x(ChCl)", 19, MUTED).next_to(ax, DOWN, buff=0.58)
        yl = ko("C=O  /  cm⁻¹", 19, MUTED).rotate(PI / 2).next_to(ax, LEFT, buff=0.62)
        tk = ticks(ax, [0, 0.2, 0.4, 0.6, 0.8, 1.0],
                   [1715, 1720, 1725, 1730], yfmt="%.0f")
        self.play(Create(ax), FadeIn(xl), FadeIn(yl), FadeIn(tk), run_time=1.1)

        g = ax.plot_line_graph([p[0] for p in CO_IR], [p[1] for p in CO_IR],
                               line_color=RED, stroke_width=4,
                               vertex_dot_style=dict(fill_color=RED,
                                                     stroke_width=0),
                               vertex_dot_radius=0.055)
        self.play(Create(g), run_time=2.0)

        band = Rectangle(width=ax.c2p(0.60, 1712)[0] - ax.c2p(0.50, 1712)[0],
                         height=ax.y_length, fill_color=WARM, fill_opacity=0.16,
                         stroke_width=0)
        band.move_to([(ax.c2p(0.50, 0)[0] + ax.c2p(0.60, 0)[0]) / 2,
                      ax.c2p(0, 1721.5)[1], 0])
        note = ko("x > 0.5 에서 기울기 급증", 21, WARM, BOLD).next_to(band, UP, buff=0.1)
        self.play(FadeIn(band), FadeIn(note), run_time=0.9)

        msg = ko("금속이 아니라 용매를 본 값 — 독립적인 제3의 증거", 24, WARM, BOLD)
        msg.next_to(xl, DOWN, buff=0.42)
        self.play(Write(msg), run_time=1.2)
        self.wait(2.0)
        self.play(FadeOut(VGroup(t, ax, tk, xl, yl, g, band, note, msg)), run_time=0.7)

    # ---------------------------------------------------------------- 5
    def s5_conclusion(self):
        head = ko("서로 다른 두 물리량이 같은 자리에서 변한다", 38, INK, BOLD).to_edge(UP, buff=1.1)
        self.play(FadeIn(head, shift=UP * 0.3), run_time=1.0)

        rows = [("육안 · 금속 화학종", "0.40 ↔ 0.60"),
                ("계산 · 금속 화학종", "0.45 ↔ 0.50"),
                ("FT-IR · 용매 구조", "0.50 ↔ 0.60")]
        g = VGroup()
        for a, b in rows:
            left = ko(a, 26, MUTED)
            right = ko(b, 26, TEAL, BOLD)
            line = VGroup(left, right).arrange(RIGHT, buff=0.9)
            g.add(line)
        g.arrange(DOWN, buff=0.42, aligned_edge=LEFT).next_to(head, DOWN, buff=0.9)
        g.set_x(0)
        for line in g:
            self.play(FadeIn(line, shift=RIGHT * 0.3), run_time=0.7)

        box = RoundedRectangle(width=10.4, height=1.5, corner_radius=0.12,
                               stroke_color=TEAL, stroke_width=2.5,
                               fill_color="#13201e", fill_opacity=1)
        concl = ko("조성비가 금속 화학종을 바꾸는 전이점은  x(ChCl) ≈ 0.5", 30, INK, BOLD)
        box.next_to(g, DOWN, buff=0.85)
        concl.move_to(box.get_center())
        self.play(FadeIn(box), Write(concl), run_time=1.5)

        limit = ko("침출효율·최적 조성은 주장하지 않는다 — 후속 과제", 20, MUTED)
        limit.next_to(box, DOWN, buff=0.4)
        self.play(FadeIn(limit), run_time=0.8)
        self.wait(2.6)
        self.play(FadeOut(VGroup(head, g, box, concl, limit)), run_time=1.0)
