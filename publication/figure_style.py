"""Shared visual specification for the manuscript's quantitative figures."""

import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap

INK = '#1E2833'
MUTED = '#5D6974'
GRID = '#E6EBEF'
BLUE = '#245792'
TEAL = '#007F73'
ROSE = '#B94F42'
VIOLET = '#705799'
GOLD = '#A8711C'
GRAY = '#78838C'
ICE = '#EBF3FB'
MINT = '#E4F4F1'

METHOD_COLORS = {
    'tmm': BLUE, 'marker_matched': TEAL,
    'marker_denominator_only': GOLD, 'oracle': INK,
    'tmm_exclude_factor': VIOLET, 'tmm_exclude_both': '#4E5862',
    'Observed total': GRAY, 'Fixed total': ROSE,
    'TMM': BLUE, 'DESeq2 poscounts': TEAL,
}
REFERENCE_COLORS = {'LIBD': BLUE, 'SEAAD': ROSE}
CLASS_COLORS = [BLUE, TEAL, GOLD, VIOLET, ROSE, '#668F69', '#A78264']
SEQUENTIAL = LinearSegmentedColormap.from_list(
    'paper3_absolute', ['#FAF8F2', '#F3DBA3', '#E9A270', '#C96251', '#77364B'])
DIVERGING = LinearSegmentedColormap.from_list(
    'paper3_difference', ['#254E89', '#90B7DB', '#FAFAF9', '#E9AE99', '#A13D38'])
DENSITY = LinearSegmentedColormap.from_list(
    'paper3_density', ['#EAF1F8', '#C5D4E4', '#829EBC', '#3F6588', '#24394F'])


def apply():
    mpl.rcParams.update({
        'font.family':'sans-serif', 'font.sans-serif':['Arial','DejaVu Sans'],
        'font.size':9, 'axes.labelsize':9, 'axes.titlesize':9,
        'axes.titleweight':'normal', 'axes.titlepad':10,
        'xtick.labelsize':8, 'ytick.labelsize':8, 'legend.fontsize':8,
        'text.color':INK, 'axes.labelcolor':INK,
        'xtick.color':INK, 'ytick.color':INK,
        'axes.spines.top':False, 'axes.spines.right':False,
        'axes.linewidth':.65, 'axes.edgecolor':'#444B52',
        'xtick.major.width':.6, 'ytick.major.width':.6,
        'xtick.major.size':3, 'ytick.major.size':3,
        'xtick.direction':'out', 'ytick.direction':'out',
        'lines.linewidth':1.45, 'lines.markersize':4.5,
        'legend.frameon':False, 'legend.handlelength':1.8,
        'legend.columnspacing':1.6, 'legend.handletextpad':.6,
        'pdf.fonttype':42, 'ps.fonttype':42, 'svg.fonttype':'none', 'svg.hashsalt':'NormCouple-Tx',
        'savefig.facecolor':'white', 'figure.facecolor':'white',
        'axes.facecolor':'white', 'axes.unicode_minus':True,
    })


def clean_axes(fig):
    for ax in fig.axes:
        if not ax.axison:
            continue
        ax.set_axisbelow(True)
        if ax.images:
            ax.tick_params(length=0)
            for spine in ax.spines.values():
                spine.set_visible(False)
        elif ax.get_label() == '<colorbar>':
            ax.tick_params(length=2, width=.5, labelsize=8)
            ax.spines['outline'].set_visible(False)
        else:
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.grid(axis='y', color=GRID, linewidth=.45)


def panel(ax, label):
    position = ax.get_position()
    left = max(-.15, (.025-position.x0)/position.width)
    ax.text(left, 1.06, label, transform=ax.transAxes,
            fontsize=11, fontweight='bold', va='bottom', color=INK)


def route_title(ax, ref, target):
    names = {'LIBD':'LIBD', 'SEAAD':'SEA-AD', 'GSE174367':'external'}
    ax.set_title(names[ref]+' → '+names[target], pad=9)


def check_text_bounds(fig):
    """Reject clipped panel labels, axis labels, annotations, and legends."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    labels = list(fig.texts)
    legends = list(fig.legends)
    for ax in fig.axes:
        labels.extend(ax.texts)
        labels.extend([ax.title, ax.xaxis.label, ax.yaxis.label])
        if ax.get_legend() is not None:
            legends.append(ax.get_legend())
    for legend in legends:
        labels.extend(legend.get_texts())
        labels.append(legend.get_title())
    outside = []
    for label in labels:
        if not label.get_visible() or not label.get_text():
            continue
        box = label.get_window_extent(renderer)
        if (box.x0 < -1 or box.y0 < -1 or
                box.x1 > fig.bbox.width+1 or box.y1 > fig.bbox.height+1):
            outside.append(label.get_text())
    if outside:
        raise ValueError('Text outside figure canvas: '+repr(outside))
