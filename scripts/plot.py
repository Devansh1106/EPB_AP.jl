# Plot CSV files written by EPB_AP (solution, snapshot or history files) with matplotlib, for a LaTeX manuscript; the
# same plots as scripts/plot.jl. Does not run the solver.
#
#   python3 scripts/plot.py data/<file>.csv [more.csv ...] [out.png] [--var rho,u,phi] [--zoom x1,x2[,y1,y2]]
#           [--inset topright|top|bottom|center|...] [--legend bottomleft|...] [--no-initial] [--no-exact]
#           [--width inches] [--no-usetex]
#
# The first column is the abscissa; every other column gets a panel (--var: only these columns, in this order), its
# axis labelled by the variable's symbol (ρ, u, ϕ, ...). The references are black lines of width 1: the exact solution
# (column `<name>_exact`) solid, the initial data (`<name>_init`) solid as well when the panel has no exact solution,
# dashed otherwise. The numerical solutions each get their own colour and marker (scripts/style.py), labelled by
# scheme (FOPD, FODD, SODD); --no-initial (or --no-init), --no-exact omit those columns. --legend places the legend
# (default best, bottomleft with --zoom).
# --zoom x1,x2 adds to each panel an inset magnifying [x1, x2] (y-range fitted; x1,x2,y1,y2 sets it, the same in all
# panels), the region marked by a grey box joined to the inset's corners; --inset places it (top right by default).
# Several files are drawn in the same panels, each curve labelled by its scheme and the metadata that differs: case, N,
# λ, t; if that is not enough, any other header entry that differs (limiter, θ, timestep, ...); else the file name.
# There is no title. Initial data are drawn once per case, exact solutions once per case and time.
# Style (scripts/style.py): Computer Modern, 9 pt labels, 8 pt ticks and legend, rendered by LaTeX (--no-usetex:
# matplotlib's own). --width is the figure width in inches (default 3.4 for one panel, 6.5 = \textwidth for more), so
# that the figure is included at its size and its text matches the manuscript's. Default output: the first file with
# .png (`_compare.png` for several), at 300 dpi; .pdf, .svg, .pgf, .eps give vector files.
# The exact solution of riemann (and riemann_speed) is the λ → 0 one: plotting it for λ > 1e-4 prints a warning.
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402

import style  # noqa: E402
from style import plt  # noqa: E402


def read_data(path):
    meta, names, rows = {}, [], []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                k, _, v = line[1:].partition("=")
                meta[k.strip()] = v.strip()
            elif not names:
                names = line.split(",")
            elif line:
                rows.append([float(v) for v in line.split(",")])
    return meta, names, np.array(rows)


KEYS = ("case", "scheme", "N", "lambda", "t")
MATH = {"N": "N", "lambda": r"\lambda", "t": "t"}


def entry(meta, k):
    v = meta.get(k, "?")
    if k == "scheme":
        return style.scheme_label(v)
    if k in MATH:
        n = style.number(v)
        return f"${MATH[k]} = {n[1:]}" if n.startswith("$") else f"${MATH[k]}$ = {n}"
    return style.text(f"{k} = {v}")


def tag(meta, keys):
    return ", ".join(entry(meta, k) for k in keys)


def file_labels(paths, files):
    """Legend label of each file: its scheme (FOPD, FODD, SODD) and the other KEYS that differ between the files; while
    two labels are the same, the other header entries that differ too (but the run's step count and status); if that is
    not enough, the file name."""
    def differs(k):
        return len({m.get(k, "?") for m, _, _ in files}) > 1

    def label(m, keys):
        return ", ".join(s for s in (entry(m, "scheme") if "scheme" in m else "", tag(m, keys)) if s)

    keys = [k for k in KEYS if differs(k) and k != "scheme"]
    labs = [label(m, keys) for m, _, _ in files]
    if len(set(labs)) < len(labs):
        more = sorted({k for m, _, _ in files for k in m if differs(k) and k not in KEYS and k not in ("steps", "status")})
        labs = [label(m, keys + more) for m, _, _ in files]
    if len(set(labs)) == len(labs) and all(labs):
        return labs
    return [style.text(os.path.basename(p)) for p in paths]


# Cases whose "exact" solution is the quasineutral (λ → 0) one, and the largest λ at which it is still a fair reference.
LIMIT_EXACT = ("riemann", "riemann_speed")
LIMIT_LAMBDA = 1e-4


def warn_limit_exact(paths, files):
    for path, (m, names, _) in zip(paths, files):
        if m.get("case", "") not in LIMIT_EXACT or not any(n.endswith("_exact") for n in names):
            continue
        try:
            lam = float(m.get("lambda", ""))
        except ValueError:
            continue
        if lam > LIMIT_LAMBDA:
            print(f"Warning: {os.path.basename(path)}: the exact solution of {m['case']} is the quasineutral (λ → 0) "
                  f"one, [1, (7.8)–(7.10)]; at λ = {lam:g} > {LIMIT_LAMBDA:g} it is not the solution of this problem "
                  "(--no-exact omits it)", file=sys.stderr)


def draw(ax, c, files, labels, init, exact, window=None):
    """Draws the panel of column c: the numerical solutions, then on top the references (initial data, dashed if an
    exact solution follows; exact). With window = (x1, x2), only the points in [x1, x2] (and one on each side).
    Returns the values drawn, for the axis limits."""
    drawn, ys = [], []

    def sel(x):
        if window is None:
            return slice(None)
        i = np.flatnonzero((x >= window[0]) & (x <= window[1]))
        return slice(0, 0) if len(i) == 0 else slice(max(i[0] - 1, 0), min(i[-1] + 1, len(x) - 1) + 1)

    def key(suffix, m):
        return (suffix, m.get("case")) if suffix == "_init" else (suffix, m.get("case"), m.get("t", ""))

    has_exact = exact and any(c + "_exact" in names for _, names, _ in files)

    def others(suffix, name):
        for i, (m, names, data) in enumerate(files):
            if c + suffix not in names or key(suffix, m) in drawn:
                continue
            drawn.append(key(suffix, m))
            j = names.index(c + suffix)
            lab = name if sum(d[0] == suffix for d in drawn) == 1 else f"{name}, {labels[i]}"
            r = sel(data[:, 0])
            style.reference(ax, data[r, 0], data[r, j], label=lab, ls="--" if suffix == "_init" and has_exact else "-")
            if suffix == "_init" or window is not None:
                ys.extend(data[r, j])

    k = 0
    for i, (m, names, data) in enumerate(files):
        if c not in names:
            continue
        j = names.index(c)
        r = sel(data[:, 0])
        ys.extend(data[r, j])
        k += 1
        style.curve(ax, data[r, 0], data[r, j], k, label=labels[i])
    if init:
        others("_init", "Initial")
    if exact:
        others("_exact", "Exact")
    return np.array(ys)


# Inset positions (x, y, width, height in the panel's axes coordinates) for --inset.
W, H, GAP = 0.4, 0.32, 0.06
INSETS = {"bottom": ((1 - W) / 2, GAP), "bottomleft": (GAP, GAP), "bottomright": (1 - GAP - W, GAP),
          "top": ((1 - W) / 2, 1 - GAP - H), "topleft": (GAP, 1 - GAP - H), "topright": (1 - GAP - W, 1 - GAP - H),
          "center": ((1 - W) / 2, (1 - H) / 2)}

# --legend positions: matplotlib's loc, and bbox_to_anchor for those outside the panel.
LEGENDS = {"best": ("best", None), "topleft": ("upper left", None), "top": ("upper center", None),
           "topright": ("upper right", None), "left": ("center left", None), "right": ("center right", None),
           "bottomleft": ("lower left", None), "bottom": ("lower center", None), "bottomright": ("lower right", None),
           "outertopright": ("upper left", (1.02, 1.0)), "outerright": ("center left", (1.02, 0.5)),
           "outerbottom": ("upper center", (0.5, -0.18)), "none": (None, None)}


def plot_files(paths, out=None, vars=None, init=True, exact=True, zoom=None, inset="topright", legend=None,
               width=None):
    out = out or re.sub(r"\.csv$", ".png" if len(paths) == 1 else "_compare.png", paths[0])
    legend = legend or ("best" if zoom is None else "bottomleft")
    files = [read_data(p) for p in paths]
    labels = file_labels(paths, files)
    if exact:
        warn_limit_exact(paths, files)
    names1 = files[0][1]
    cols = [n for n in names1 if n != names1[0] and not n.endswith("_exact") and not n.endswith("_init")]
    if vars is not None:
        bad = [v for v in vars if v not in cols]
        if bad:
            sys.exit(f"unknown variable(s) {', '.join(bad)} in --var; available: {', '.join(cols)}")
        cols = vars
    nc = min(len(cols), 3)
    nr = -(-len(cols) // nc)                                      # near-square panels, at most 3 per row
    width = width or (3.4 if nc == 1 else 6.5)
    fig, axes = plt.subplots(nr, nc, figsize=(width, 0.9 * width / nc * nr), layout="constrained", squeeze=False)
    for ax in axes.flat[len(cols):]:
        ax.set_visible(False)
    for ax, c in zip(axes.flat, cols):
        ax.set_xlabel(style.symbol(names1[0]))
        ax.set_ylabel(style.symbol(c))
        ys = draw(ax, c, files, labels, init, exact)
        ys = ys[np.isfinite(ys)]
        if c in ("dt", "ap_residual") and len(ys) and np.all(ys > 0) and ys.max() >= 10 * ys.min():  # a decade or more
            ax.set_yscale("log")
        if len(ys):
            lo, hi = ys.min(), ys.max()
            w = max(1.0, abs(hi))
            if hi - lo < 1e-9 * w:                                # constant up to round-off (e.g. mass)
                ax.set_ylim(lo - 1e-6 * w, hi + 1e-6 * w)
        ax.set_xlim(min(d[:, 0].min() for _, _, d in files), max(d[:, 0].max() for _, _, d in files))
        if zoom is not None:
            # [x1, x2] magnified, its y-range [y1, y2] or fitted to the curves there; a grey box marks the region,
            # joined to the inset's corners.
            axins = ax.inset_axes([*INSETS[inset], W, H])
            ys = draw(axins, c, files, labels, init, exact, window=zoom[:2])
            ys = ys[np.isfinite(ys)]
            if len(ys) == 0:
                sys.exit(f"no data of {c} in --zoom {zoom[0]:g},{zoom[1]:g}")
            lo, hi = zoom[2:] if len(zoom) == 4 else (ys.min(), ys.max())
            pad = 0.0 if len(zoom) == 4 else max(0.05 * (hi - lo), 1e-12 * max(1.0, abs(hi)))
            axins.set_xlim(*zoom[:2])
            axins.set_ylim(lo - pad, hi + pad)
            axins.tick_params(labelsize=6, length=2)
            ind = ax.indicate_inset_zoom(axins, edgecolor="0.45", alpha=1.0, linewidth=0.6)
            for cn in ind.connectors:
                cn.set(color="0.45", linewidth=0.5, alpha=1.0)
        loc, anchor = LEGENDS[legend]
        if loc is not None:
            ax.legend(loc=loc, bbox_to_anchor=anchor)
    fig.savefig(out)
    print("wrote", out)


FLAGS = ("--no-init", "--no-initial", "--no-exact", "--no-usetex")
USAGE = ("usage: python3 scripts/plot.py data/<file>.csv [more.csv ...] [out.png] [--var rho,u,phi] "
         "[--zoom x1,x2[,y1,y2]] [--inset topright|top|bottom|center|...] [--legend bottomleft|...] "
         "[--no-init(ial)] [--no-exact] [--width inches] [--no-usetex]")


def main(argv):
    argv, vals = list(argv), {}
    for k in ("--var", "--zoom", "--inset", "--legend", "--width"):
        if k not in argv:
            continue
        i = argv.index(k)
        if i == len(argv) - 1 or argv[i + 1].startswith("--"):
            sys.exit(f"{k} needs a value\n{USAGE}")
        vals[k] = argv[i + 1]
        del argv[i:i + 2]
    vars = [v.strip() for v in vals["--var"].split(",") if v.strip()] if "--var" in vals else None
    zoom = tuple(float(v) for v in vals["--zoom"].split(",")) if "--zoom" in vals else None
    if zoom is not None and not (len(zoom) in (2, 4) and zoom[0] < zoom[1] and (len(zoom) == 2 or zoom[2] < zoom[3])):
        sys.exit(f"--zoom takes x1,x2 (x1 < x2) or x1,x2,y1,y2 (y1 < y2: one y-range for all panels)\n{USAGE}")
    inset = vals.get("--inset", "topright")
    if inset not in INSETS:
        sys.exit(f"--inset is one of {', '.join(sorted(INSETS))}\n{USAGE}")
    legend = vals.get("--legend")
    if legend is not None and legend not in LEGENDS:
        sys.exit(f"--legend is one of {', '.join(LEGENDS)}\n{USAGE}")
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
    plot_files(inputs, out=outs[0] if outs else None, vars=vars, zoom=zoom, inset=inset, legend=legend, width=width,
               init="--no-init" not in flags and "--no-initial" not in flags, exact="--no-exact" not in flags)


if __name__ == "__main__":
    main(sys.argv[1:])
