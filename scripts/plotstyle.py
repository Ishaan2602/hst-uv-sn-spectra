# shared figure style for every product plot: the paper-grade look (serif, inward major+minor ticks on
# all four sides, thick spines). call apply() once right after importing pyplot so the per-SN diagnostic
# pngs match the published-figure standard instead of shipping as a default matplotlib dump.
import matplotlib as mpl

_STYLE = {
    "font.family": "serif",
    "font.size": 12,
    "axes.linewidth": 1.2,
    "axes.grid": False,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "xtick.major.size": 6, "ytick.major.size": 6,
    "xtick.minor.size": 3, "ytick.minor.size": 3,
    "legend.frameon": False,
    "savefig.bbox": "tight",
}


def apply():
    mpl.rcParams.update(_STYLE)
