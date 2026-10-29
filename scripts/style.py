# Plot style of the Python scripts, as scripts/style.jl, sized for a LaTeX manuscript: Computer Modern text (rendered
# by LaTeX, or by matplotlib's mathtext with --no-usetex), 9 pt labels, 8 pt ticks and legend, full frame with ticks
# inside, no grid, no titles. A reference (exact solution, initial data) is a black line of width 1; every other curve
# has its own colour (the palette of the Julia plots) and marker and width 1.2.
import shutil

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

COLORS = ("#009AFA", "#E36F47", "#3EA44E", "#C371D2", "#AC8E18", "#00AAAE", "#ED5E93", "#C68225", "#00A98D",
          "#8E971D")
MARKERS = ("o", "s", "^", "D", "v", "*", "p", "h", "X", "P")
REF_LW = 1.0
LW = 1.2

# Axis symbol of a solution variable (CSV column), in LaTeX.
SYMBOLS = {"x": r"$x$", "t": r"$t$", "rho": r"$\rho$", "u": r"$u$", "phi": r"$\phi$", "m": r"$\rho u$",
           "energy": r"$E$", "min_rho": r"$\min \rho$", "mass": "mass", "ap_residual": r"$\Vert e^\phi - \rho\Vert_\infty$",
           "dt": r"$\Delta t$", "cfl_ratio": "CFL ratio"}

# Legend name of a scheme (`scheme` header entry): shift (§3.3), dissipation (§3.4), second order (§4).
SCHEME_LABELS = {"shift": "FOPD", "dissipation": "FODD", "second": "SODD"}

_usetex = False


def apply(usetex=None):
    """Sets the rcParams; usetex: LaTeX text (default: if `latex` is installed)."""
    global _usetex
    _usetex = shutil.which("latex") is not None if usetex is None else usetex
    plt.rcParams.update({
        "text.usetex": _usetex,
        "text.latex.preamble": r"\usepackage{amsmath}",
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman", "cmr10"],
        "mathtext.fontset": "cm",
        "axes.formatter.use_mathtext": True,
        "axes.formatter.useoffset": False,
        "axes.unicode_minus": True,
        "font.size": 9,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "axes.linewidth": 0.6,
        "axes.grid": False,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "legend.frameon": True,
        "legend.fancybox": False,
        "legend.framealpha": 1.0,
        "legend.edgecolor": "black",
        "legend.borderpad": 0.4,
        "legend.handlelength": 2.2,
        "patch.linewidth": 0.5,
        "lines.markersize": 3.5,
        "savefig.dpi": 300,
    })


def text(s):
    """Plain text for a label: LaTeX's special characters escaped when LaTeX renders it."""
    if not _usetex:
        return s.replace("$", r"\$")
    for a, b in (("\\", r"\textbackslash{}"), ("_", r"\_"), ("%", r"\%"), ("&", r"\&"), ("#", r"\#"),
                 ("$", r"\$"), ("{", r"\{"), ("}", r"\}"), ("^", r"\^{}"), ("~", r"\~{}")):
        s = s.replace(a, b)
    return s


def symbol(name):
    return SYMBOLS.get(name, text(name))


def scheme_label(name):
    return SCHEME_LABELS.get(name, text(name))


def number(v):
    """A header value as LaTeX math: 1e-4 is 10^{-4}, 2.5e-6 is 2.5 \\times 10^{-6}; other text is escaped."""
    try:
        x = float(v)
    except ValueError:
        return text(v)
    if x == 0 or not np.isfinite(x) or 1e-2 <= abs(x) < 1e4:
        return f"${x:g}$"
    e = int(np.floor(np.log10(abs(x))))
    m = x / 10.0**e
    return f"$10^{{{e}}}$" if np.isclose(m, 1) else f"${m:g} \\times 10^{{{e}}}$"


def reference(ax, x, y, label="reference", ls="-", **kw):
    return ax.plot(x, y, color="black", lw=REF_LW, ls=ls, label=label, **kw)


def curve(ax, x, y, k, label, nmarks=15, **kw):
    """Curve number k (1, 2, ...): its colour, the line through all points, its marker on at most nmarks of them
    (staggered with k, so that curves on the same grid do not hide each other's markers); one legend entry."""
    x, y = np.asarray(x), np.asarray(y)
    n = len(x)
    if n <= nmarks:
        idx = np.arange(n)
    else:
        step = (n - 1) / (nmarks - 1)
        idx = np.unique(np.clip(np.round(np.mod((k - 1) * 0.3 * step, step) + step * np.arange(nmarks)), 0, n - 1)
                        .astype(int))
    idx = [i for i in idx if np.isfinite(y[i])]
    return ax.plot(x, y, color=COLORS[(k - 1) % len(COLORS)], lw=LW, marker=MARKERS[(k - 1) % len(MARKERS)],
                   markevery=idx, markeredgecolor="black", markeredgewidth=0.4, label=label, **kw)
