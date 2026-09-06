"""Shared matplotlib style for all figures in the JFM manuscript
(manuscript/main_v3.tex, class JFM-FLM_Au).

House rule (Sparsh, 2026-09-06): every plot in the paper must match the
paper's own typography — same font family (CUP Times) and the same
printed size for ticks, labels and any annotated text.  The way to
guarantee the printed size is to BUILD each figure at its true printed
width and include it at natural size:

    \\includegraphics[width=\\textwidth]{...}   ->  figsize width = FULL
    width=0.5\\textwidth                        ->  HALF

Never build wide and let LaTeX scale down — that shrinks the fonts.

Fonts: Times New Roman (the closest system match to CUP Times) with the
STIX math set, which is Times-compatible.  Base size 9 pt = the printed
size of JFM figure lettering (between \\small captions and 10 pt body).
No in-figure titles or suptitles: the caption carries that text in JFM.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: F401  (re-export for callers)

PT = 1.0 / 72.27          # TeX point in inches
FULL = 32 * 12 * PT       # \textwidth = 32 pc = 5.313 in
HALF = 0.5 * FULL

matplotlib.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "STIXGeneral", "Nimbus Roman",
                   "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "legend.frameon": False,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.minor.width": 0.45,
    "ytick.minor.width": 0.45,
    "lines.linewidth": 1.1,
    "lines.markersize": 3.5,
    "grid.color": "0.92",
    "grid.linewidth": 0.5,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,   # embed TrueType so Times survives in the PDF
})


def strip(ax, grid_axis="y"):
    """JFM axes hygiene: no top/right spines, light horizontal grid."""
    ax.spines[["top", "right"]].set_visible(False)
    if grid_axis:
        ax.grid(axis=grid_axis, zorder=0)


def save(fig, path_noext):
    """Write .pdf (manuscript) and .png (email/preview) at natural size."""
    for ext in ("pdf", "png"):
        fig.savefig(path_noext + "." + ext, bbox_inches="tight")
    print("saved", path_noext + ".pdf/.png")
