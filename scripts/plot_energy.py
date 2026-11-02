# Plot the energy curves written by scripts/energy.jl --csv with matplotlib, for a LaTeX manuscript: the plot of
# energy.jl --plot. Does not run the solver.
#
#   python3 scripts/plot_energy.py <energy>.csv [out.png] [--legend best|...] [--width inches] [--no-usetex]
#
# E(t) - E(0), or (E(0) - E(t))/E(0) for energy.jl --relative, against t. From a case study: one panel per scheme
# (FOPD, FODD, SODD, named in its legend title), one curve per λ. From history files: one panel, one curve per run
# (FOPD, FODD, SODD, ...), and the energy of the exact solution, if the case has one, as the black reference.
# --legend places the legend in each panel (the positions of plot.py; default best). Style, --width and --no-usetex
# as in plot.py. Default output: the file with .png.
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import style  # noqa: E402
from plot import LEGENDS, read_data  # noqa: E402
from style import plt  # noqa: E402

YLABELS = {"relative": r"$(E(0) - E(t))/E(0)$", "difference": r"$E(t) - E(0)$"}


def plot_energy(path, out=None, legend="best", width=None):
    meta, names, d = read_data(path)
    if names != ["curve", "lambda", "t", "value"]:
        sys.exit(f"{path} is not a curve file of energy.jl --csv")
    ks = list(dict.fromkeys(int(k) for k in d[:, 0]))
    # A case study has one panel per scheme; history files one panel.
    panels = list(dict.fromkeys(meta[f"panel{k}"] for k in ks if f"panel{k}" in meta)) or [None]
    n = len(panels)
    nc = 2 if n == 4 else min(n, 3)
    nr = -(-n // nc)
    width = width or (3.4 if nc == 1 else 6.5)
    fig, axes = plt.subplots(nr, nc, figsize=(width, 0.9 * width / nc * nr), layout="constrained", squeeze=False)
    for ax in axes.flat[n:]:
        ax.set_visible(False)
    for ax, p in zip(axes.flat, panels):
        ax.set_xlabel("$t$")
        ax.set_ylabel(YLABELS[meta.get("measure", "difference")])
        j = 0
        for k in ks:
            if p is not None and meta.get(f"panel{k}") != p:
                continue
            rows = d[d[:, 0] == k]
            if k == 0:                                            # the exact solution
                style.reference(ax, rows[:, 2], rows[:, 3], label="Exact")
                continue
            j += 1
            lab = f"$\\lambda = {style.number(rows[0, 1])[1:]}" if p is not None else style.text(meta[f"curve{k}"])
            style.curve(ax, rows[:, 2], rows[:, 3], j, lab)
        loc, anchor = LEGENDS[legend]
        if loc is not None:
            ax.legend(loc=loc, bbox_to_anchor=anchor, title=style.text(p) if p is not None else None)
    out = out or re.sub(r"\.csv$", ".png", path)
    fig.savefig(out)
    print("wrote", out)


FLAGS = ("--no-usetex",)
USAGE = "usage: python3 scripts/plot_energy.py <energy>.csv [out.png] [--legend best|...] [--width inches] [--no-usetex]"


def main(argv):
    argv, vals = list(argv), {}
    for k in ("--legend", "--width"):
        if k not in argv:
            continue
        i = argv.index(k)
        if i == len(argv) - 1 or argv[i + 1].startswith("--"):
            sys.exit(f"{k} needs a value\n{USAGE}")
        vals[k] = argv[i + 1]
        del argv[i:i + 2]
    legend = vals.get("--legend", "best")
    if legend not in LEGENDS:
        sys.exit(f"--legend is one of {', '.join(LEGENDS)}\n{USAGE}")
    width = float(vals["--width"]) if "--width" in vals else None
    flags = [a for a in argv if a.startswith("--")]
    args = [a for a in argv if not a.startswith("--")]
    for f in flags:
        if f not in FLAGS:
            sys.exit(f"unknown option {f}\n{USAGE}")
    isout = re.compile(r"\.(png|pdf|svg|pgf|eps)$").search
    inputs, outs = [a for a in args if not isout(a)], [a for a in args if isout(a)]
    if len(inputs) != 1 or len(outs) > 1:
        sys.exit(USAGE)
    style.apply(usetex=False if "--no-usetex" in flags else None)
    plot_energy(inputs[0], out=outs[0] if outs else None, legend=legend, width=width)


if __name__ == "__main__":
    main(sys.argv[1:])
