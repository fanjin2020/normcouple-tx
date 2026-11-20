"""Empirical-anchor quantitative grids, using the established matplotlib workflow.

Contract: describe the size/direction of fixed-marker disease associations,
not normalization superiority. Panel evidence = coefficient distribution,
point-estimate dose fractions, and donor-level marker library uncertainty.
PDF/SVG keep editable text; all numbers are in sibling CSV source files.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from anchor_core import CLASSES


def save(fig, final, name):
    for ext in ('pdf','svg','png'):
        fig.savefig(final/(name+'.'+ext),dpi=300)
    plt.close(fig)


def make_figures(summary, genes, libraries, final):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.spines.top':False,
                         'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    primary = genes[(genes.model=='primary')&(genes.normalizer=='TMM_all')&(genes.scope=='owner')]
    fig,axes = plt.subplots(1,2,figsize=(10.0,5.2),layout='constrained',sharex=True)
    colors = {'LIBD':'#007C83','SEAAD':'#AD4A73'}
    for ax,(ref,color) in zip(axes,colors.items()):
        for j,ct in enumerate(CLASSES):
            data = primary[(primary.reference==ref)&(primary.cell_type==ct)]
            val = pd.to_numeric(data.log2FC,errors='coerce').dropna().to_numpy() if 'log2FC' in data else np.array([])
            if len(val):
                jitter = np.linspace(-.17,.17,len(val))
                ax.scatter(val,j+jitter,s=9,color=color,alpha=.45,linewidths=0,rasterized=True)
                ax.plot(np.median(val),j,'|',color='black',ms=12,mew=1.3)
            else:
                ax.text(0,j,'Unavailable',ha='center',va='center',fontsize=7,color='#666666')
        for value in (-1,-np.log2(1.5),0,np.log2(1.5),1):
            ax.axvline(value,color='#777777' if value==0 else '#BBBBBB',lw=.7,zorder=0,
                       ls='-' if value==0 else ':')
        ax.set_yticks(range(len(CLASSES)),CLASSES)
        ax.invert_yaxis()
        ax.set_title(ref+' reference markers',loc='left',fontweight='bold')
        ax.set_xlabel('AD vs control log2 fold-change')
    save(fig,final,'Fig_anchor_gene_effects')

    s = summary[(summary.model=='primary')&(summary.normalizer=='TMM_all')&(summary.scope=='owner')]
    labels = [ref+' / '+ct for ref in colors for ct in CLASSES]
    fig,axes = plt.subplots(1,2,figsize=(9.8,6.1),layout='constrained')
    for ax,col,title in zip(axes,['point_up_fraction_tested','point_down_fraction_tested'],['Higher in AD','Lower in AD']):
        grid = np.full((len(labels),3),np.nan)
        for j,label in enumerate(labels):
            ref,ct = label.split(' / ')
            for k,dose in enumerate((25,50,100)):
                f = s[(s.reference==ref)&(s.cell_type==ct)&(s.dose_percent==dose)]
                if len(f):
                    grid[j,k] = f[col].iloc[0]*100
        cmap = plt.get_cmap('YlGnBu').copy();cmap.set_bad('#E9E9E9')
        im = ax.imshow(np.ma.masked_invalid(grid),vmin=0,vmax=100,cmap=cmap,aspect='auto')
        for (j,k),v in np.ndenumerate(grid):
            ax.text(k,j,'NA' if not np.isfinite(v) else f'{v:.0f}',ha='center',va='center',
                    fontsize=7,color='white' if np.isfinite(v) and v>55 else 'black')
        ax.set_xticks(range(3),['1.25','1.50','2.00'])
        ax.set_yticks(range(len(labels)),labels)
        ax.set_title(title,loc='left',fontweight='bold')
        ax.set_xlabel('Fold-change magnitude threshold')
    fig.colorbar(im,ax=axes,label='Point estimates beyond threshold (% of tested markers)',shrink=.8)
    save(fig,final,'Fig_anchor_dose_context')

    d = libraries[(libraries.model=='primary')&(libraries.scope=='owner')&
                  (libraries.endpoint=='marker_fraction_log2')].copy()
    fig,ax = plt.subplots(figsize=(8.5,6.0),layout='constrained')
    for j,label in enumerate(labels):
        ref,ct=label.split(' / ')
        r=d[(d.reference==ref)&(d.cell_type==ct)&d.available]
        if len(r):
            v=r.iloc[0]
            ax.errorbar(v.beta,j,xerr=[[v.beta-v.lower],[v.upper-v.beta]],fmt='o',ms=4,
                        color=colors[ref],capsize=2,lw=1)
        else:
            ax.text(0,j,'Unavailable',ha='center',va='center',fontsize=7,color='#666666')
    ax.axvline(0,color='#777777',lw=.8,ls='--')
    ax.set_yticks(range(len(labels)),labels);ax.invert_yaxis()
    ax.set_xlabel('AD association with log2 marker UMI fraction (95% HC3 interval)')
    save(fig,final,'Fig_anchor_marker_library')
    (final/'FIGURE_LEGENDS.md').write_text('''# Empirical marker anchor figures

## Gene effects
Each dot is a fixed reference marker evaluated in its assigned broad cell type;
black ticks show the median among tested genes, not confidence intervals.
Dotted lines mark +/-log2(1.5) and +/-1. Primary age/sex/batch-adjusted donor
pseudobulk edgeR fits use TMM factors from expressed genes. Donor numbers,
expression-filter exclusions and untested planned markers are reported in CSV.

## Dose context
Cells show the percentage of tested fixed markers whose signed point estimates
exceed the indicated fold-change magnitude; these are NOT significant fractions
or prevalence estimates with known truth. NA denotes an unavailable design.
The corresponding threshold-test support, using edgeR glmTreat and BH over the
full planned panel family, is in marker_dose_anchor_summary.csv. No independence
of genes is assumed to construct binomial intervals.

## Marker library
Points are covariate-adjusted donor-level AD coefficients for log2((marker
UMI+0.5)/(total UMI+1)); bars are 95% HC3 t intervals. Both counts use exactly the
same author-annotated nuclei. This is a relative captured-library endpoint, not
an absolute RNA-yield measurement. No nuclei are treated as biological repeats.

None of these observational plots estimates true FDR, establishes the cause of
a marker change, or tests whether marker-restricted normalization is superior.
''',encoding='utf-8')
