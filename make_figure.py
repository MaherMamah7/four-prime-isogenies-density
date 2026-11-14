#!/usr/bin/env python3
"""
Figure: log10 Pr[R(E) > 0] against the bit length of p, one panel per lambda,
with the observed fraction of successes and the lower bound of Proposition 1.

  python3 make_figure.py                      # uses data/results.csv, data/lower_bound.csv
  python3 make_figure.py --results results.csv --bound lower_bound.csv --lambdas 8 10

Inputs are the outputs of `density.py merge --out` and `density.py bound --out`.
Requires matplotlib.
"""
import argparse
import csv
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter

ap = argparse.ArgumentParser()
ap.add_argument("--results", default="data/results.csv")
ap.add_argument("--bound", default="data/lower_bound.csv")
ap.add_argument("--lambdas", type=float, nargs="+", default=[8, 10])
ap.add_argument("--out", default="density_figure")
a = ap.parse_args()

obs = [r for r in csv.DictReader(open(a.results))]
low = [r for r in csv.DictReader(open(a.bound))]

INK, INK2, GRID, BLUE, GREEN = "#0b0b0b", "#52514e", "#e4e3df", "#2a78d6", "#1baf7a"
plt.rcParams.update({"font.size": 9, "font.family": "serif", "axes.edgecolor": INK2,
                     "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2})
fig, axes = plt.subplots(1, len(a.lambdas), figsize=(3.5 * len(a.lambdas), 3.5), sharey=True, squeeze=False)
xmax = max(int(r["bits"]) for r in obs)
for ax, lam in zip(axes[0], a.lambdas):
    lb = sorted((int(r["bits"]), float(r["lower_bound"])) for r in low
                if float(r["lambda"]) == lam and int(r["bits"]) <= xmax + 4)
    pts = sorted((int(r["bits"]), float(r["pr"]), float(r["se"])) for r in obs if float(r["lambda"]) == lam)
    if lb:
        ax.plot([b for b, _ in lb], [v for _, v in lb], color=GREEN, lw=2)
        ax.text(lb[-1][0] + 2.5, lb[-1][1], "proven\nlower bound", color=INK, va="center", fontsize=7.5)
    if pts:
        x = [b for b, _, _ in pts]; y = [v for _, v, _ in pts]
        ax.plot(x, y, color=BLUE, lw=1.6, ls=":")
        ax.errorbar(x, y, yerr=[2 * s for _, _, s in pts], fmt="o", color=BLUE, ms=6, mec="white",
                    mew=1.2, elinewidth=1.2, capsize=3, zorder=4)
        ax.text(x[-1] + 2.5, y[-1], "observed", color=INK, va="center", fontsize=7.5)
    ax.set_yscale("log"); ax.set_ylim(1e-12, 3)
    ax.yaxis.set_major_locator(FixedLocator([10.0 ** k for k in range(-12, 1, 2)]))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{round(math.log10(v))}"))
    ax.yaxis.set_minor_locator(FixedLocator([]))
    lo = min([b for b, _ in lb] + [b for b, _, _ in pts])
    ax.set_xlim(lo - 5, xmax + 15)
    ax.set_xticks(list(range(20 * math.ceil((lo - 5) / 20), xmax + 1, 20)))
    ax.set_xlabel(r"bit length of $p$")
    ax.set_title(rf"$\lambda={lam:g}$", fontsize=9, color=INK, loc="left")
    ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
axes[0][0].set_ylabel(r"$\log_{10}\Pr_E[R(E)>0]$")
fig.tight_layout()
fig.savefig(a.out + ".pdf"); fig.savefig(a.out + ".png", dpi=200)
print(f"wrote {a.out}.pdf and {a.out}.png")
