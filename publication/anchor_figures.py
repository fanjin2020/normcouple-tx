"""Figure for the held-out anchor experiment and observed assay pairs."""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from . import figure_style as style
from pipeline.cbc_anchors import LABELS, METHODS


def draw(summary, pairs, output, formats, formal=True):
    style.apply()
    output = Path(output)
    output.mkdir(parents=True,exist_ok=True)
    fig = plt.figure(figsize=(7.1,5.65))
    grid = fig.add_gridspec(2,1,left=.15,right=.95,bottom=.10,top=.91,
                            height_ratios=[1,1.18],hspace=.67)
    ax = fig.add_subplot(grid[0])
    p = pairs.assign(change=100*pairs.total_relative_change).sort_values('change')
    y = np.arange(len(p))
    ax.hlines(y,0,p.change,color=style.GRAY,lw=1.15,zorder=2)
    ax.scatter(p.change,y,c=np.where(p.change>=0,style.BLUE,style.TEAL),s=25,zorder=3)
    ax.axvline(0,ls=':',color=style.INK,lw=.8)
    blocks = p.block.str.replace('_',' ',regex=False)
    ax.set_yticks(y,[f'{b} / r{r}' for b,r in zip(blocks,p['round'])],fontsize=7)
    ax.set(xlabel='Counterfactual library-total change after equal-depth substitution (%)',
           title=f'Observed test-count contrast: {len(p)} protocol pairs in {p.donor.nunique()} donors')
    style.panel(ax,'a')
    sub=grid[1].subgridspec(1,3,wspace=.26,width_ratios=[1,1,1])
    palette=[style.BLUE,style.GRAY,style.TEAL,style.GOLD,style.ROSE,style.INK]
    for j,(metric,outcome,title,threshold) in enumerate([
        ('fdr','null','Global-null FDR',5),('power','sparse','Power',None),
        ('coverage','sparse','95% interval coverage',95)]):
        ax=fig.add_subplot(sub[j])
        focus=summary.loc[(summary.probe=='marker') & (summary.outcome==outcome)]
        for i,(method,color) in enumerate(zip(METHODS,palette)):
            for direction,mk,off in [('Astro','o',.165),('Micro','s',-.165)]:
                r=focus.loc[(focus.direction==direction) & (focus.method==method)].iloc[0]
                m,lo,hi=100*r[metric],100*r[metric+'_lo'],100*r[metric+'_hi']
                ax.errorbar(m,len(METHODS)-1-i+off,xerr=[[max(0,m-lo)],[max(0,hi-m)]],
                            color=color,fmt=mk,lw=1.0,ms=4.3,
                            mfc='white' if mk=='s' else color,capsize=1.8)
        ax.set_yticks(range(len(METHODS)),LABELS[::-1] if j==0 else [])
        ax.tick_params(axis='y',length=0)
        if j==0:
            ax.set_xlim(left=0)
        if threshold is not None:
            ax.axvline(threshold,c=style.MUTED,lw=.8,ls=':')
        ax.set(xlabel='%',title=title,ylim=(-.47,len(METHODS)-.53))
        style.panel(ax,'bcd'[j])
    fig.text(.58,.024,'Filled circle: astrocyte shift   |   Open square: microglial shift',
             ha='center',fontsize=8,color=style.MUTED)
    if not formal:
        fig.text(.51,.984,'SOFTWARE ACCEPTANCE: NOT RESEARCH RESULTS',ha='center',va='top',fontsize=8)
    style.clean_axes(fig)
    style.check_text_bounds(fig)
    for fmt in tuple(dict.fromkeys(formats)):
        kwargs={'dpi':300,'facecolor':'white'}
        if fmt=='pdf':
            kwargs['metadata']={'CreationDate':None}
        fig.savefig(output/f'FigS12.{fmt}',**kwargs)
    plt.close(fig)
