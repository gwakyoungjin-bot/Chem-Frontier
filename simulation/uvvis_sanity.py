"""2차 UV-Vis 5개 몰비 데이터 건전성 점검. python uvvis_sanity.py"""
import glob, os, statistics as st, math
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

def load(f):
    x, on = {}, False
    for ln in open(f, encoding='utf-8', errors='replace'):
        ln = ln.strip()
        if ln == 'XYDATA': on = True; continue
        if not on: continue
        p = ln.split(',')
        if len(p) != 2: break
        try: x[float(p[0])] = float(p[1])
        except ValueError: break
    return x

D = {os.path.basename(f)[:-4].replace('_', ':'): load(f)
     for f in sorted(glob.glob('1;4.csv_260920/*.csv'))}
mean = lambda r, a, b: st.mean(v for w, v in r.items() if a <= w <= b)

base = {k: mean(r, 600, 700) for k, r in D.items()}   # 흡수체 없음 → 0이어야 함
peak = {k: mean(r, 335, 345) for k, r in D.items()}
assert all(abs(v) > 0.03 for v in base.values()), "베이스라인이 0이면 이 점검은 불필요"

r = (lambda a, b: (lambda ma, mb: sum((i-ma)*(j-mb) for i, j in zip(a, b)) /
     math.sqrt(sum((i-ma)**2 for i in a)*sum((j-mb)**2 for j in b)))(st.mean(a), st.mean(b))
     )(list(base.values()), list(peak.values()))

print(f"{'몰비':8}{'raw340':>9}{'base650':>9}{'보정340':>9}")
for k in D: print(f"{k:8}{peak[k]:+9.3f}{base[k]:+9.3f}{peak[k]-base[k]:+9.3f}")
print(f"\nPearson r(베이스라인 오프셋, 340nm 피크) = {r:.3f}")
print("raw 순위 :", ' > '.join(sorted(D, key=peak.get, reverse=True)))
print("보정 순위:", ' > '.join(sorted(D, key=lambda k: peak[k]-base[k], reverse=True)))

# --- 엑셀 요약표(218-221nm 피크)가 단일점 노이즈 스파이크인지 판정 ---
XL_PEAK = {'1:4': 221, '1:19': 221, '3:1': 219, '3:17': 221, '19:1': 218}
print("\n엑셀 보고피크 vs 이웃 2-3nm 평균 (실제 밴드면 차이 ~0):")
for k, p in XL_PEAK.items():
    nb = st.mean(D[k][p+d] for d in (-3, -2, 2, 3))
    print(f"  {k:6} λ={p} A={D[k][p]:.3f} 이웃={nb:.3f} 차이={D[k][p]-nb:+.3f}")
print("\n5점 평활 후 230-500nm λmax (실제 밴드면 살아남음):")
for k, r in D.items():
    sm = {w: st.mean(r[w+d] for d in (-2, -1, 0, 1, 2)) for w in range(232, 499)}
    mx = max(sm, key=sm.get)
    print(f"  {k:6} λmax={mx} A={sm[mx]-base[k]:+.3f}")
print("\n검출기 포화(A=7) 상한 -> 이 아래 전부 무효:")
for k, r in D.items():
    print(f"  {k:6} {max(w for w in r if r[w] >= 6.9):.0f} nm")

fig, ax = plt.subplots(1, 2, figsize=(11, 4))
for k, s in D.items():
    w = sorted(s, reverse=True)
    ax[0].plot(w, [s[i] for i in w], lw=.7, label=k)
    ax[1].plot(w, [s[i]-base[k] for i in w], lw=.7, label=k)
for a, t in zip(ax, ['raw', 'baseline-corrected (600-700nm 기준)']):
    a.axvline(340, color='k', ls=':', lw=.8); a.set_xlim(250, 800); a.set_ylim(-.4, 1)
    a.set_title(t); a.set_xlabel('nm'); a.set_ylabel('A'); a.legend(fontsize=7)
ax[0].annotate('lamp change\n@340nm', (340, .85), fontsize=7)
plt.tight_layout(); plt.savefig('uvvis_sanity.png', dpi=130)
print("\n-> uvvis_sanity.png")
