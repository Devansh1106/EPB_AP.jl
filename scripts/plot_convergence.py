# Plot the errors written by scripts/convergence.jl --csv with matplotlib, for a LaTeX manuscript: the plot of
# convergence.jl --plot. Does not run the solver.
#
#   python3 scripts/plot_convergence.py <errors>.csv [more.csv ...] [out.png] [--var rho,u,phi] [--norm L1,L2,Linf]
#           [--lambda λ[,λ2,...]] [--slope 1,2] [--legend best|below|...] [--width inches] [--no-usetex]
#
# Log–log errors against N, one panel per variable (--var: only these, in this order; default all of the first file),
# one curve per (file, λ, norm) in each panel, each with its own colour and marker (scripts/style.py). --norm keeps
# only these norms (default all three), --lambda only these λ (default all of each file). Curves are labelled by what
# differs between them: the scheme (FOPD, FODD, SODD; with the limiter or η_c if two files differ in it), λ, the norm.
# Black reference lines of the slopes in --slope (default: the orders of the files' schemes, 1 or 2), solid for the
# first slope and dashed, dotted for the next ones, start at half the lowest error on the coarsest grid. The panel of
# ap (the AP residual e^φ - ρ, which has no order) has instead a dotted line at λ² for each λ.
# Several files (e.g. FOPD and SODD runs of one case) are drawn in the same panels. --legend places the legend in each
# panel (the positions of plot.py; default best), or `below`: one legend for all panels, under them (for many curves). Style, --width and --no-usetex as in plot.py. Default output: the first file
# with .png (`_compare.png` for several), with _<norms> before it when --norm leaves out a norm (e.g. _L1).
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402

import style  # noqa: E402
from plot import LEGENDS, read_data  # noqa: E402
from style import plt  # noqa: E402

NORMS = ("L1", "L2", "Linf")
NORM_LABELS = {"L1": "$L^1$", "L2": "$L^2$", "Linf": r"$L^\infty$"}
LINESTYLES = ("-", "--", ":", "-.")


def variables(names):
    return [n[:-3] for n in names if n.endswith("_L1")]


def scheme_tag(meta, limiter, eta):
    s = style.scheme_label(meta.get("scheme", "?"))
    if limiter and "limiter" in meta:
        return f"{s} ({style.text(meta['limiter'])})"
    return rf"{s} ($\eta_c = {meta['eta_c']}$)" if eta and "eta_c" in meta else s


def plot_errors(paths, out=None, vars=None, norms=NORMS, lambdas=None, slopes=None, legend="best", width=None):
    files = [read_data(p) for p in paths]
    for p, (_, names, _) in zip(paths, files):
        if names[:2] != ["lambda", "N"]:
            sys.exit(f"{p} is not an error file of convergence.jl --csv")
    vars = vars or variables(files[0][1])
    for p, (_, names, _) in zip(paths, files):
        bad = [v for v in vars if v not in variables(names)]
        if bad:
            sys.exit(f"unknown variable(s) {', '.join(bad)} in {p}; available: {', '.join(variables(names))}")
    # One curve per (file, λ, norm): (file index, λ, norm, rows of that λ).
    curves = []
    for i, (_, _, d) in enumerate(files):
        for lam in dict.fromkeys(d[:, 0]):
            if lambdas is None or np.any(np.isclose(lam, lambdas, rtol=1e-9, atol=0)):
                curves += [(i, lam, n, d[d[:, 0] == lam]) for n in norms]
    if not curves:
        sys.exit("no data of the --lambda values")
    limiter = len({m["limiter"] for m, _, _ in files if "limiter" in m}) > 1
    eta = len({float(m["eta_c"]) for m, _, _ in files if "eta_c" in m}) > 1
    many = (len({c[0] for c in curves}) > 1, len({c[1] for c in curves}) > 1, len(norms) > 1)

    def label(i, lam, n):
        parts = [scheme_tag(files[i][0], limiter, eta) if many[0] else None,
                 f"$\\lambda = {style.number(lam)[1:]}" if many[1] else None,
                 NORM_LABELS[n] if many[2] or not any(many) else None]
        return ", ".join(p for p in parts if p)

    if slopes is None:
        slopes = sorted({2.0 if m.get("scheme") == "second" else 1.0 for m, _, _ in files})
    Ns = sorted({int(N) for c in curves for N in c[3][:, 1]})
    n = len(vars)
    nc = 2 if n == 4 else min(n, 3)
    nr = -(-n // nc)
    width = width or (3.4 if nc == 1 else 6.5)
    fig, axes = plt.subplots(nr, nc, figsize=(width, 0.9 * width / nc * nr), layout="constrained", squeeze=False)
    for ax in axes.flat[n:]:
        ax.set_visible(False)
    for ax, v in zip(axes.flat, vars):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("$N$")
        ax.set_ylabel(style.symbol(v) if v == "ap" else f"error in {style.symbol(v)}")
        cols = {i: files[i][1] for i in range(len(files))}
        data = [(i, lam, nm, rows[:, 1], rows[:, cols[i].index(f"{v}_{nm}")]) for i, lam, nm, rows in curves]
        x = np.array(Ns, dtype=float)
        if v == "ap":                                             # no order: the level λ² instead
            lams = list(dict.fromkeys(c[1] for c in curves))
            for lam in lams:
                lab = r"$\lambda^2$" if len(lams) == 1 else rf"$\lambda^2$, $\lambda = {style.number(lam)[1:]}"
                style.reference(ax, x, np.full(len(x), lam**2), label=lab, ls=":")
        else:
            y0 = 0.5 * min(e[0] for *_, e in data)
            for j, q in enumerate(slopes):
                style.reference(ax, x, y0 * (x[0] / x) ** q, label=f"slope {q:g}", ls=LINESTYLES[j % len(LINESTYLES)])
        for k, (i, lam, nm, N, e) in enumerate(data, start=1):
            style.curve(ax, N, e, k, label(i, lam, nm))
        ax.set_xticks(Ns, [str(N) for N in Ns])
        ax.xaxis.set_minor_locator(plt.NullLocator())
        loc, anchor = LEGENDS.get(legend, (None, None))
        if loc is not None:
            ax.legend(loc=loc, bbox_to_anchor=anchor)
    if legend == "below":                                         # one legend for all panels, under them
        entries = {}                                              # every panel's entries, each once
        for ax in axes.flat[:n]:
            for h, lab in zip(*ax.get_legend_handles_labels()):
                entries.setdefault(lab, h)
        fig.legend(entries.values(), entries.keys(), loc="outside lower center",
                   ncol=min(len(entries), 4 if nc > 1 else 2))
    if out is None:
        nrm = "_" + "_".join(norms) if len(norms) < len(NORMS) else ""
        out = re.sub(r"\.csv$", ("" if len(paths) == 1 else "_compare") + nrm + ".png", paths[0])
    fig.savefig(out)
    print("wrote", out)


FLAGS = ("--no-usetex",)
USAGE = ("usage: python3 scripts/plot_convergence.py <errors>.csv [more.csv ...] [out.png] [--var rho,u,phi] "
         "[--norm L1,L2,Linf] [--lambda λ[,λ2,...]] [--slope 1,2] [--legend best|below|...] [--width inches] [--no-usetex]")


def main(argv):
    argv, vals = list(argv), {}
    for k in ("--var", "--norm", "--lambda", "--slope", "--legend", "--width"):
        if k not in argv:
            continue
        i = argv.index(k)
        if i == len(argv) - 1 or argv[i + 1].startswith("--"):
            sys.exit(f"{k} needs a value\n{USAGE}")
        vals[k] = argv[i + 1]
        del argv[i:i + 2]

    def items(k):
        return [v.strip() for v in vals[k].split(",") if v.strip()] if k in vals else None

    vars = items("--var")
    norms = items("--norm") or list(NORMS)
    if not all(n in NORMS for n in norms):
        sys.exit(f"--norm takes a list of {', '.join(NORMS)}\n{USAGE}")
    lambdas = [float(v) for v in items("--lambda")] if "--lambda" in vals else None
    slopes = [float(v) for v in items("--slope")] if "--slope" in vals else None
    legend = vals.get("--legend", "best")
    if legend not in LEGENDS and legend != "below":
        sys.exit(f"--legend is one of {', '.join(LEGENDS)}, below\n{USAGE}")
    width = float(vals["--width"]) if "--width" in vals else None
    flags = [a for a in argv if a.startswith("--")]
    args = [a for a in argv if not a.startswith("--")]
    for f in flags:
        if f not in FLAGS:
            sys.exit(f"unknown option {f}\n{USAGE}")
    isout = re.compile(r"\.(png|pdf|svg|pgf|eps)$").search
    inputs, outs = [a for a in args if not isout(a)], [a for a in args if isout(a)]
    if not inputs or len(outs) > 1:
        sys.exit(USAGE)
    style.apply(usetex=False if "--no-usetex" in flags else None)
    plot_errors(inputs, out=outs[0] if outs else None, vars=vars, norms=norms, lambdas=lambdas, slopes=slopes,
                legend=legend, width=width)


if __name__ == "__main__":
    main(sys.argv[1:])
