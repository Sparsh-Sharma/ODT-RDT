"""Shared matplotlib style for all figures in the JFM manuscript
(manuscript/main_v3.tex, class JFM-FLM_Au).

CANON (Sparsh's SoundPower/LE-noise papers, module jfm_rapids.py; house
rule 2026-09-06: use this style whenever plotting for JFM):
  - Times-like serif (STIX, bundled with matplotlib, no LaTeX needed)
    for text AND maths; base 8 pt, axis labels 9 pt, never below 7 pt.
  - Thin BOXED axes (0.6 pt), ticks INWARD on all four sides, NO grid.
  - Panel labels italic (a), (b), ... via panel().
  - Build at the TRUE printed width and include at 100 %:
    \\includegraphics[width=\\textwidth] -> figsize width FULL
    (JFM text width = 32 pc = 5.31 in); never build wide and let LaTeX
    scale down — that shrinks the fonts.
  - save(): vector PDF (600 dpi raster parts) + PNG proof.

COLOUR (consistent across the whole paper — same entity, same colour):
  systems   : DNS black, exact linear RDT purple, standard ODT blue.
  components: phi_22/upwash blue, phi_11 orange, phi_33 green
              (within-system component plots).
  variants  : ODT intervention family ramp (VARIANTS list), baseline
              black.
  contrasts : plane strain orange vs axisymmetric blue; slow blue vs
              rapid orange.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: F401  (re-export for callers)

PT = 1.0 / 72.27          # TeX point in inches
FULL = 32 * 12 * PT       # \textwidth = 32 pc = 5.31 in
HALF = 0.5 * FULL

# --- the paper-wide palette ---------------------------------------------
BLUE = "#2a78d6"
ORANGE = "#eb6834"
GREEN = "#1baf7a"
PURPLE = "#8e44ad"
GOLD = "#d29a00"
RED = "#c22f2f"
BLACK = "k"

COL = {
    # systems
    "dns": BLACK, "rdt": PURPLE, "odt": BLUE,
    # components (within one system)
    "phi22": BLUE, "phi11": ORANGE, "phi33": GREEN,
    # strain-geometry contrast
    "plane": ORANGE, "axi": BLUE, "b22": GREEN,
    # strain-rapidity contrast
    "slow": BLUE, "rapid": ORANGE,
}
# ODT intervention family (baseline first, then the ramp)
VARIANTS = [BLACK, BLUE, GREEN, ORANGE, PURPLE, GOLD]


def set_jfm():
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["STIXGeneral", "Times New Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 8,
        "axes.titlesize": 8.5,
        "axes.labelsize": 9,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7,
        "axes.linewidth": 0.6,
        "lines.linewidth": 1.1,
        "lines.markersize": 3.2,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "xtick.major.size": 2.8, "ytick.major.size": 2.8,
        "xtick.minor.size": 1.5, "ytick.minor.size": 1.5,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.grid": False,
        "legend.frameon": False, "legend.handlelength": 1.5,
        "legend.borderpad": 0.2, "legend.labelspacing": 0.25,
        "savefig.dpi": 600, "savefig.bbox": "tight",
        "savefig.pad_inches": 0.01,
        "figure.dpi": 150,
        "pdf.fonttype": 42,
    })


set_jfm()


def panel(ax, lb, x=0.02, y=0.97, color="k"):
    """Italic JFM panel label (a), (b), ... in the axes corner."""
    ax.text(x, y, f"({lb})", transform=ax.transAxes, fontstyle="italic",
            fontsize=8.5, va="top", ha="left", color=color,
            bbox=dict(boxstyle="square,pad=0.12", fc="white", ec="none",
                      alpha=0.75))


def save(fig, path_noext):
    """Write .pdf (manuscript, 600 dpi raster parts) + .png proof."""
    fig.savefig(path_noext + ".pdf")
    fig.savefig(path_noext + ".png", dpi=200)
    print("saved", path_noext + ".pdf/.png")
