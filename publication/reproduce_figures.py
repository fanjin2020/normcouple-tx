"""Reproduce all manuscript figures from the supplied numerical source tables."""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline.publication import figure_index, load_frame
from pipeline.io import separate

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from matplotlib.ticker import FormatStrFormatter, MaxNLocator, NullLocator
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm

import figure_panels as core
import figure_style as style
import figure_data as fd

HERE = Path(__file__).resolve().parent
CAPTIONS = json.loads((HERE/'figure_captions.json').read_text(encoding='utf-8'))
COL = style.METHOD_COLORS
METHODS = ['tmm','tmm_exclude_factor','tmm_exclude_both','marker_denominator_only','marker_matched']
LABELS = ['TMM','Factor-only','Complete\nexclusion','Denom.','MRN']
ROUTES, SHORT = core.ROUTES, core.SHORT
MANIFEST = {}
READS = set()
SOFTWARE_ACCEPTANCE = False


def read(name):
    READS.add(name)
    return load_frame(DATA, name)


def finish(fig, name, caption=None):
    entry = INDEX[name]
    if READS != set(entry['views']):
        raise ValueError(f"{name}: plotted inputs {sorted(READS)} differ from the figure index {sorted(entry['views'])}")
    if SOFTWARE_ACCEPTANCE:
        fig.text(.5,.99,'SOFTWARE ACCEPTANCE: NOT RESEARCH RESULTS',ha='center',va='top',fontsize=6,color='#666666')
    style.clean_axes(fig)
    style.check_text_bounds(fig)
    for fmt in FORMATS:
        kwargs = {'dpi':300 if fmt=='png' else 600, 'facecolor':'white'}
        if fmt=='pdf':
            kwargs['metadata']={'CreationDate':None}
        elif fmt=='svg':
            kwargs['metadata']={'Date':None}
        if fmt=='tiff':
            kwargs['pil_kwargs']={'compression':'tiff_lzw'}
        fig.savefig(OUT/(name+'.'+fmt), **kwargs)
        if fmt=='tiff':
            target=OUT/(name+'.'+fmt)
            with Image.open(target) as source:
                rgb=source.convert('RGB')
            rgb.save(target,compression='tiff_lzw',dpi=(600,600))
    MANIFEST[name]={'caption':caption or CAPTIONS[name], **entry,
                    'width_inches':fig.get_figwidth(), 'height_inches':fig.get_figheight()}
    plt.close(fig)


def select(frame, **conditions):
    for key,val in conditions.items():
        frame=frame[frame[key].eq(val)]
    if len(frame)!=1:
        raise ValueError((conditions,len(frame)))
    return frame.iloc[0]


def figure5():
    f=read('common_support_summary.csv')
    fig,axs=plt.subplots(2,2,figsize=(6.85,4.95))
    fig.subplots_adjust(left=.145,right=.98,bottom=.18,top=.91,hspace=.67,wspace=.57)
    for ax,route,tag,letter in zip(axs.flat,ROUTES,SHORT,'abcd'):
        for j,kind in enumerate(['legacy','matched','random']):
            for k,method in enumerate(['tmm','marker_matched']):
                y=2-j+(k-.5)*.25
                values=[100*select(f,reference_source=route[0],target_source=route[1],
                    subset_kind=kind,support_rule=rule,outcome='null',method=method).fdp
                   for rule in ['allow_zero','positive_only']]
                marker='o' if k==0 else 's'
                ax.plot(values,[y,y],color=COL[method],lw=1.4,zorder=2)
                ax.plot(values[0],y,color=COL[method],marker=marker,ms=5,mfc='white',mew=1.1,zorder=3)
                ax.plot(values[1],y,color=COL[method],marker=marker,ms=5,mfc=COL[method],mew=.8,zorder=4)
        ax.set_yticks([2,1,0],['Predefined','Matched','Random'])
        ax.set_ylim(-.48,2.48);ax.set_xlim(0,90);ax.set_xticks([0,30,60,90])
        ax.set_xlabel('Global-null FDR (%)')
        ax.axvline(5,color=style.MUTED,ls=':',lw=.8)
        style.route_title(ax,*route)
        core.letter(ax,letter)
    handles=[Line2D([],[],color=COL['tmm'],marker='o',label='TMM'),
             Line2D([],[],color=COL['marker_matched'],marker='s',label='MRN'),
             Line2D([],[],color=style.INK,marker='o',mfc='white',ls='none',label='Zero allowed'),
             Line2D([],[],color=style.INK,marker='o',ls='none',label='Positive only')]
    fig.legend(handles=handles,loc='lower center',ncol=4,bbox_to_anchor=(.5,.02))
    finish(fig,'Fig5')


def figure6():
    points=fd.subset_points(read('common_support_donor_subsets.csv'))
    fig=plt.figure(figsize=(6.85,5.32))
    comparison=fig.add_gridspec(1,2,left=.207,right=.962,top=.867,bottom=.515,wspace=.28,
                                width_ratios=[1,1.62])
    other=fig.add_gridspec(1,3,left=.125,right=.962,top=.353,bottom=.175,wspace=.30)
    ref,target=ROUTES[2]
    focus=points[points.reference_source.eq(ref)&points.target_source.eq(target)].set_index('base_subset')
    random=sorted(k for k in focus.index if k.startswith('random_') and k!='random_mean')
    keys=['legacy','random_mean',*random]
    labels=['Predefined','Random mean',*[f'Random {j}' for j in range(1,len(random)+1)]]
    marks=['*','D',*(['o']*len(random))]
    colors=[style.TEAL,style.ROSE,*([style.GRAY]*len(random))]
    fig.text(.204,.963,'LIBD reference · 15 external target donors',weight='bold',fontsize=10.5)
    full_fdr=float(100*points.fdp.max())
    for j,(column,ax,tag) in enumerate(zip(['proportion_rmse','fdp'],[fig.add_subplot(c) for c in comparison],'ab')):
        v=focus.loc[keys,column].to_numpy()*(100 if column=='fdp' else 1)
        y=np.arange(len(v))[::-1]
        ax.axhspan(y[0]-.42,y[0]+.42,color=style.MINT,alpha=.82,lw=0)
        ax.axhspan(y[1]-.42,y[1]+.42,color='#FBEDE8',alpha=.82,lw=0)
        for idx,(x,yi,marker,color) in enumerate(zip(v,y,marks,colors)):
            ax.plot(x,yi,marker=marker,ms=9 if idx==0 else 6,c=color,ls='none',
                    mfc='white' if idx==1 else color,mew=1.35 if idx==1 else .75)
            if column=='fdp' or idx<2:
                ax.annotate(f'{x:.2f}%' if column=='fdp' else f'{x:.4f}',(x,yi),
                            xytext=(7,0),textcoords='offset points',va='center',
                            fontsize=8.1,weight='bold' if idx<2 else 'normal',color=color)
        ax.set_ylim(-.55,len(y)-.45)
        ax.set_yticks(y,labels if j==0 else [])
        ax.tick_params(axis='y',length=0)
        ax.set_xlabel('Proportion RMSE' if j==0 else 'Global-null FDR (%)')
        if j==0:
            span=np.ptp(v)
            ax.set_xlim(v.min()-.14*span,v.max()+.55*span)
            ax.xaxis.set_major_locator(MaxNLocator(3))
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
        else:
            ax.set_xlim(-2,max(60,v.max()*1.22))
            ax.axvline(5,color=style.MUTED,lw=.9,ls=':')
        core.letter(ax,tag)
    styles=[('zero','Unperturbed',style.GRAY,'x'),('matched_','Matched',style.BLUE,'s'),
            ('random_','Random',style.GRAY,'o'),('legacy','Predefined',style.TEAL,'*'),
            ('random_mean','Random mean',style.ROSE,'D')]
    for i,(cell,route) in enumerate(zip(other,[ROUTES[0],ROUTES[1],ROUTES[3]])):
        ax=fig.add_subplot(cell)
        g=points[points.reference_source.eq(route[0])&points.target_source.eq(route[1])]
        for prefix,_,color,marker in styles:
            mask=g.base_subset.str.startswith(prefix) if prefix.endswith('_') else g.base_subset.eq(prefix)
            if prefix=='random_':
                mask &= ~g.base_subset.eq('random_mean')
            v=g[mask]
            ax.plot(v.proportion_rmse,100*v.fdp,marker=marker,ls='none',color=color,
                    ms=7 if marker=='*' else 4.7,mfc='white' if marker=='D' else color,mew=.95)
        ax.axhline(5,color=style.MUTED,lw=.8,ls=':')
        ax.set_ylim(-2,max(60,full_fdr*1.05))
        ax.yaxis.set_major_locator(MaxNLocator(3))
        ax.xaxis.set_major_locator(MaxNLocator(2,prune='both'))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
        ax.set_xlabel('Proportion RMSE')
        if i==0:
            ax.set_ylabel('FDR (%)')
        else:
            ax.tick_params(axis='y',labelleft=False)
        style.route_title(ax,*route)
        core.letter(ax,'cde'[i])
    handles=[Line2D([],[],label=label,marker=marker,color=color,ls='none',
                    mfc='white' if marker=='D' else color,ms=6) for _,label,color,marker in styles]
    fig.legend(handles=handles,ncol=5,loc='lower center',bbox_to_anchor=(.52,.012),
               columnspacing=1.05,handlelength=1.05)
    finish(fig,'Fig6')


def figure8():
    f=read('clinical_coefficient_pairs.csv.gz'); s=read('clinical_stability.csv')
    fig=plt.figure(figsize=(6.85,5.68))
    grid=fig.add_gridspec(2,2,left=.13,right=.875,bottom=.135,top=.91,
                          height_ratios=[1.45,1],wspace=.49,hspace=.63)
    axs=[fig.add_subplot(cell) for cell in grid]
    refs=['LIBD','SEAAD']; counts={}; artists=[]
    both=f[['tmm_beta','mrn_beta']].to_numpy()
    low,high=both.min(),both.max()
    pad=.048*(high-low); lim=[low-pad,high+pad]
    keys={'TMM only':(style.ROSE,'^'),'MRN only':(style.TEAL,'o')}
    for ax,ref in zip(axs[:2],refs):
        g=f[f.reference.eq(ref)]
        cats,counts[ref]=fd.clinical_status(g,select(s,reference=ref,method='marker_matched'))
        hexb=ax.hexbin(g.tmm_beta,g.mrn_beta,gridsize=49,extent=[*lim,*lim],mincnt=1,
                       cmap=style.DENSITY,linewidths=0,zorder=1)
        artists.append(hexb)
        ax.plot(lim,lim,c=style.INK,lw=.8,ls='--',alpha=.70,zorder=2)
        for key,(color,marker) in keys.items():
            v=g[cats==key]
            ax.scatter(v.tmm_beta,v.mrn_beta,s=7 if key=='TMM only' else 9,marker=marker,
                       c=color,edgecolors='none',alpha=.76,rasterized=True,zorder=3)
        ax.set(xlim=lim,ylim=lim,xlabel='TMM disease coefficient',ylabel='MRN disease coefficient',
               title=ref.replace('SEAAD','SEA-AD')+' reference',aspect='equal')
    maximum=max(a.get_array().max() for a in artists)
    for a in artists:
        a.set_norm(LogNorm(vmin=1,vmax=maximum))
    cax=fig.add_axes([.90,.552,.017,.309])
    cb=fig.colorbar(artists[0],cax=cax)
    cb.set_label('Genes per bin (log scale)',labelpad=5)
    ax=axs[2]
    y=np.arange(3)[::-1]
    for j,ref in enumerate(refs):
        for k,(status,color) in enumerate([('Both',style.GRAY),('TMM only',style.ROSE),('MRN only',style.TEAL)]):
            x=counts[ref][status]
            yy=y[k]+(.14 if j==0 else -.14)
            ax.barh(yy,x,height=.23,fc=color if j==0 else 'white',ec=color,lw=.85,zorder=3)
            ax.annotate(f'{x:,}',(x,yy),xytext=(4,0),textcoords='offset points',
                        va='center',color=style.INK,fontsize=7.1)
    ymax=max(counts[ref][k] for ref in refs for k in ('Both','TMM only','MRN only'))
    ax.set(xlim=(0,1.26*ymax),ylim=(-.48,2.47),xlabel='Significant genes',title='Discovery membership')
    ax.set_yticks(y,['Both','TMM only','MRN only'])
    ax.xaxis.set_major_locator(MaxNLocator(3))
    ax.tick_params(axis='y',length=0)
    ax=axs[3]
    for ref,c,ls in [('LIBD',style.BLUE,'-'),('SEAAD',style.TEAL,'--')]:
        g=f[f.reference.eq(ref)];ratio=np.sort((g.mrn_se/g.tmm_se).to_numpy())
        ax.plot(ratio,np.arange(1,len(ratio)+1)/len(ratio),color=c,ls=ls,label=ref.replace('SEAAD','SEA-AD'))
    ax.axvline(1,color=style.MUTED,ls=':',lw=.8)
    ax.set(xlabel='MRN / TMM standard error',ylabel='Fraction of genes',title='Precision')
    ax.set_xscale('log');ax.set_xticks([1,1.2,1.4,1.6,1.8]);ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.xaxis.set_minor_locator(NullLocator());ax.legend(loc='lower right')
    fig.text(.52,.03,'Filled bars: LIBD    |    Outlined bars: SEA-AD',fontsize=8,ha='center',color=style.MUTED)
    for ax,l in zip(axs,'abcd'):
        core.letter(ax,l)
    finish(fig,'Fig8')


def figure9():
    f=read('clinical_finite.csv')
    measured='clinical_probes.csv' in VIEWS
    fig,axs=plt.subplots(2,2,figsize=(6.85,6.15))
    fig.subplots_adjust(left=.14,right=.97,bottom=.155 if measured else .12,top=.90,hspace=.89,wspace=.48)
    dose=read('clinical_probes.csv') if measured else None
    if measured:
        dose=fd.fusiform_doses(f,dose)
    refs=['LIBD','SEAAD']
    arms=[('tmm','o','TMM','-'),('tmm_exclude_factor','^','Factor-only','--'),
          ('tmm_exclude_both','s','Three invariant arms',':')]
    for j,ref in enumerate(refs):
        for i,metric in enumerate(['threshold_crossing_fraction','coefficient_change_rms'][:1 if measured else 2]):
            ax=axs[i+int(measured),j]; sub=f[f.reference.eq(ref)]
            base=sub[sub.method.eq('tmm_exclude_both')].sort_values('log_count_step')
            for other in ('marker_denominator_only','marker_matched'):
                v=sub[sub.method.eq(other)].sort_values('log_count_step')
                if (base.empty or not np.array_equal(base.log_count_step,v.log_count_step) or
                        not (np.equal(base[metric].to_numpy(),0).all() and np.equal(v[metric].to_numpy(),0).all())):
                    raise ValueError('The pooled zero trace requires three equal, zero-valued paired contrasts')
            for m,mk,label,ls in arms:
                v=sub[sub.method.eq(m)].sort_values('log_count_step')
                x=100*np.expm1(v.log_count_step);y=v[metric]*(100 if i==0 else 1)
                ax.plot(x,y,marker=mk,ms=4.6,color=COL[m],ls=ls,label=label)
                if i==0 and m=='tmm':
                    ax.annotate(f'{y.iloc[-1]:.2f}%',(x.iloc[-1],y.iloc[-1]),xytext=(-3,7),
                                textcoords='offset points',ha='right',color=style.INK,weight='bold',fontsize=8)
            ax.set_xticks([0,5,10]);ax.set_xlabel('Tested-count increase (%)')
            ax.set_ylabel('BH crossings (%)' if i==0 else 'Coefficient-change RMS')
            if i==0:
                ax.set_ylim(-1,1.25*100*f.threshold_crossing_fraction.max()+1)
                ax.set_title(ref.replace('SEAAD','SEA-AD')+' reference')
    fig.legend(*axs[int(measured),0].get_legend_handles_labels(),ncol=3,loc='center',
               bbox_to_anchor=(.54,.038 if measured else .512),
               handlelength=1.8,columnspacing=1.1)
    if measured:
        ax=axs[0,1]
        for k,level in enumerate(('low','middle','high')):
            for i,ref in enumerate(refs):
                for direction,mk,offset in [('up','^',-.065),('down','v',.065)]:
                    row=select(dose,probe=f'single_{level}_{direction}')
                    value=row[ref+'_singleton_beta_change']*1e5
                    if not np.isfinite(value): raise ValueError('Missing single-gene coefficient change')
                    ax.scatter(k+(i-.5)*.36+offset,value,marker=mk,facecolors=style.BLUE if i==0 else 'white',
                               s=38,edgecolors=style.BLUE if i==0 else style.ROSE,linewidths=1,zorder=3)
        ax.axhline(0,color=style.MUTED,lw=.8);ax.set_xlim(-.55,2.55)
        shares=[select(dose,probe=f'single_{k}_up').AD_count_share_pct for k in ('low','middle','high')]
        ax.set_xticks(range(3),[f'{s:.3g}%' for s in shares]); ax.set_xlabel('Single gene: share of AD counts')
        ax.set_ylabel(r'Own coefficient change ($10^{-5}$ units)')
        own=dose[dose.probe.str.startswith('single_')][[ref+'_own_BH_crossings' for ref in refs]]
        if own.shape!=(6,2) or not own.isin([0,1]).all().all():
            raise ValueError('Expected six paired singleton interventions under both references')
        ax.set_title(f'Single-gene probe · {int(own.to_numpy().sum())}/{own.size} own calls crossed')
        ax.text(.5,-.46,'Triangles: ▲ increase, ▼ decrease\nFilled: LIBD   |   Open: SEA-AD',
                ha='center',va='top',transform=ax.transAxes,fontsize=7.9,color=style.MUTED)
        ax=axs[0,0]
        for j,(group,color,mk) in enumerate([('AD',style.BLUE,'o'),('control',style.ROSE,'s')]):
            for k,direction in enumerate(('down','up')):
                row=select(dose,probe='test_family_'+direction)
                q1,med,q3=[row[f'{group}_effective_{v}_pct'] for v in ('q25','median','q75')]
                ax.errorbar(k+(j-.5)*.18,med,yerr=[[med-q1],[q3-med]],color=color,fmt=mk,ms=5,
                            capsize=3,label=('AD' if j==0 else 'Control') if not k else None)
                if direction=='up':
                    ax.annotate(f'{med:+.2f}%',(k+(j-.5)*.18,q3),xytext=(-3 if j==0 else 4,3),
                                textcoords='offset points',ha='right' if j==0 else 'left',
                                color=color,fontsize=8,weight='bold')
        ticks=[select(dose,probe='test_family_'+t).AD_input_change_pct for t in ('down','up')]
        ax.set_xticks([0,1],[f'{x:+.2f}%' for x in ticks]);ax.set_xlim(-.55,1.55)
        yl=ax.get_ylim();ax.set_ylim(yl[0],yl[1]+.16*(yl[1]-yl[0]))
        ax.axhline(0,color=style.MUTED,lw=.8);ax.set_xlabel('Full-family count perturbation in AD')
        ax.set_ylabel('Effective denominator Δ (%)')
        share=select(dose,probe='test_family_up').AD_count_share_pct
        ax.set_title(f'Full family: {share:.2f}% of AD counts')
        ax.legend(loc='lower center',ncol=2,bbox_to_anchor=(.5,-.55),columnspacing=.8)
    for ax,l in zip(axs.flat,'abcd'):core.letter(ax,l)
    caption=None if measured else ('Fig. 9 Clinical dose response in computed runs without measured probe evidence. '
        'Panels a,b give BH threshold crossings; c,d give the RMS coefficient change at the same three doses '
        'under LIBD and SEA-AD references. Complete exclusion, denominator-only and MRN coincide at zero. '
        'Measured singleton and effective-denominator results require the paired clinical evidence directory')
    finish(fig,'Fig9',caption)


def figureS2():
    f=read('dose_retention.csv');c=read('dose_by_class.csv')
    fig,axs=plt.subplots(2,1,figsize=(6.85,3.6),gridspec_kw={'height_ratios':[.58,1.65]})
    fig.subplots_adjust(left=.12,right=.90,bottom=.14,top=.895,hspace=.59)
    vals=[100*select(f,reference_source=r,target_source=t).retained_fraction for r,t in ROUTES]
    ax=axs[0];ax.axis('off');ax.set(xlim=(0,1),ylim=(0,1))
    ax.text(0,.82,'Pooled counts retained',va='center',weight='bold',fontsize=9)
    ax.plot([0,.99],[.57,.57],c=style.INK,lw=.65)
    for x,tag,v in zip([.005,.26,.51,.765],SHORT,vals):
        ax.text(x,.19,tag,va='center',fontsize=8)
        ax.text(x+.226,.19,f'{v:.2f}%',ha='right',va='center',fontsize=9,weight='bold')
    m=np.array([[100*select(c,reference_source=r,target_source=t,cell_class=k).retained_fraction for k in range(1,8)] for r,t in ROUTES])
    im=axs[1].imshow(m,aspect='auto',cmap=style.SEQUENTIAL,vmin=0,vmax=100)
    axs[1].set_xticks(np.arange(-.5,7,1),minor=True);axs[1].set_yticks(np.arange(-.5,4,1),minor=True)
    axs[1].grid(which='minor',color='white',linewidth=.8);axs[1].tick_params(which='minor',length=0)
    axs[1].set_xticks(range(7),['Astro','EndoMural','Excit','Inhib','Micro','OPC','Oligo'],fontsize=8)
    axs[1].set_yticks(range(4),SHORT)
    for (i,j),v in np.ndenumerate(m):
        rgb=np.array(im.cmap(v/100)[:3]);lin=np.where(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)**2.4)
        lum=float(lin@np.array([.2126,.7152,.0722]));color='black' if (lum+.05)/.05>=1.05/(lum+.05) else 'white'
        axs[1].text(j,i,f'{v:.1f}',ha='center',va='center',color=color,fontsize=8)
    fig.colorbar(im,ax=axs[1],pad=.025,fraction=.035,label='Retained (%)')
    for ax,l in zip(axs,'ab'):core.letter(ax,l)
    finish(fig,'FigS2')


def figureS3():
    f=read('common_support_contrasts.csv')
    f=f[f.comparison.eq('marker_matched_minus_marker_denominator_only')&f.support_rule.eq('positive_only')]
    fig,axs=plt.subplots(1,3,figsize=(6.85,3.05))
    fig.subplots_adjust(left=.1,right=.98,bottom=.26,top=.88,wspace=.48)
    for ax,metric,outcome,factor,yl,letter in zip(axs,['fdp','power','proportion_rmse'],['null','up025','null'],
            [100,100,1],['Δ FDR (pp)','Δ power (pp)','Δ proportion RMSE'],'abc'):
        for j,kind,color,marker in [(0,'legacy',style.ROSE,'*'),(1,'matched',style.TEAL,'D'),(2,'random',style.VIOLET,'o')]:
            y=[factor*select(f,reference_source=r,target_source=t,subset_kind=kind,outcome=outcome)[metric] for r,t in ROUTES]
            ax.scatter(np.arange(4)+(j-1)*.16,y,c=color,marker=marker,s=30,label={'legacy':'Predefined','matched':'Matched','random':'Random'}[kind])
        ax.axhline(0,ls=':',color='#777777',lw=.8);core.routes_axes(ax);ax.set_ylabel(yl);core.letter(ax,letter)
    fig.legend(*axs[0].get_legend_handles_labels(),ncol=3,loc='lower center',frameon=False)
    finish(fig,'FigS3')


def figureS4():
    f=read('clinical_design.csv');fig,axs=plt.subplots(1,2,figsize=(6.85,3.25))
    fig.subplots_adjust(left=.11,right=.97,bottom=.32,top=.9,wspace=.4)
    for ax,metric,yl,l in zip(axs,['exposure_r2','exposure_residual_ss'],['Disease–nuisance R²','Residual exposure sum of squares'],'ab'):
        for ref,c,marker in [('LIBD',style.BLUE,'o'),('SEAAD',style.ROSE,'s')]:
            vals=[select(f,reference=ref,method=m)[metric] for m in METHODS]
            ax.plot(range(5),vals,color=c,marker=marker,ls='none',label=ref.replace('SEAAD','SEA-AD'))
        ax.set_xticks(range(5),LABELS,rotation=35,ha='right',fontsize=8);ax.set_ylabel(yl);core.letter(ax,l)
    fig.legend(*axs[0].get_legend_handles_labels(),loc='lower center',ncol=2,frameon=False)
    finish(fig,'FigS4')


def figureS7():
    q=read('diagnostic_qualification.csv');d=read('clinical_derivatives.csv')
    fig=plt.figure(figsize=(6.85,5.5));gs=fig.add_gridspec(2,2,height_ratios=[1.2,1],hspace=.78,wspace=.44,left=.12,right=.98,top=.91,bottom=.18)
    a=fig.add_subplot(gs[0,:]);b=fig.add_subplot(gs[1,0]);c=fig.add_subplot(gs[1,1])
    arms=METHODS+['nnls_linear','nnls_marker_matched']
    for j,col,color,label in [(0,'planned_records','#D7DDE2','Planned'),(1,'original_finite_scores',style.BLUE,'Finite'),(2,'qualified_scores',style.TEAL,'Qualified')]:
        vals=[select(q,method=m,outcome='null')[col] for m in arms]
        a.bar(np.arange(7)+(j-1)*.24,vals,width=.23,color=color,edgecolor='white',linewidth=.5,label=label)
    a.set_xticks(range(7),['TMM','Factor-only','Complete','Denom.','MRN','Linear','Linear-M'],fontsize=8)
    a.set_ylabel('Records per outcome');a.set_ylim(0,1600)
    a.legend(frameon=False,ncol=3,fontsize=8,loc='upper right',bbox_to_anchor=(1,1.19))
    for ax,ref in [(b,'LIBD'),(c,'SEAAD')]:
        for m,label,marker in [('tmm','TMM','o'),('tmm_exclude_factor','Factor-only','^')]:
            g=d[d.reference.eq(ref)&d.method.eq(m)].sort_values('derivative_step')
            ax.plot(g.derivative_step,g.derivative_agreement,marker=marker,c=COL[m],label=label)
        ax.axhline(.1,ls=':',color='#777777',lw=.8);ax.set_xscale('log');ax.set_yscale('log')
        ax.set_xlabel('Central-difference step');ax.set_ylabel('Relative disagreement')
        ax.text(.5,1.035,ref.replace('SEAAD','SEA-AD'),transform=ax.transAxes,va='bottom',ha='center')
    for ax,l in zip([a,b,c],'abc'):core.letter(ax,l)
    fig.legend(*b.get_legend_handles_labels(),loc='lower center',ncol=2,frameon=False)
    finish(fig,'FigS7')


def figureS8():
    f=read('cross_method.csv')
    fig,axs=plt.subplots(2,2,figsize=(6.85,4.9))
    fig.subplots_adjust(left=.13,right=.98,bottom=.20,top=.91,hspace=.58,wspace=.38)
    methods=[('dtangle_tmm','dtangle + TMM',style.BLUE,'o'),
             ('dtangle_marker','dtangle + MRN',style.TEAL,'s'),
             ('music_marker','MuSiC, fixed markers',style.ROSE,'^')]
    if f.method.eq('linear_svr_marker').any():
        methods.append(('linear_svr_marker','Linear SVR, fixed markers',style.VIOLET,'D'))
    for row,effect in enumerate([0.,.5]):
        for col,probe in enumerate(['test_counts','marker_counts']):
            ax=axs[row,col]
            for j,(method,label,color,marker) in enumerate(methods):
                g=[select(f,reference_source=r,target_source=t,effect=effect,method=method,probe=probe) for r,t in ROUTES]
                ax.errorbar(np.arange(4)+(j-(len(methods)-1)/2)*.16,[z.beta_rms for z in g],
                    yerr=[z.beta_rms_mc_se for z in g],fmt=marker,color=color,ms=4.5,
                    capsize=2,elinewidth=.8,label=label,linestyle='none')
            ax.axhline(0,color='#A7ADB2',lw=.6,zorder=0)
            ax.set_xticks(range(4),SHORT,fontsize=8)
            ax.set_ylabel('Coefficient-change RMS')
            ax.set_title(('Tested-count probe' if col==0 else 'Marker-count probe')+
                         (' | null' if row==0 else ' | effect 0.5'))
            ax.ticklabel_format(axis='y',style='sci',scilimits=(-3,3),useMathText=True)
    for ax,label in zip(axs.flat,'abcd'):core.letter(ax,label)
    fig.legend(*axs[0,0].get_legend_handles_labels(),loc='lower center',ncol=2,frameon=False,fontsize=8)
    finish(fig,'FigS8')


def figureS9():
    f=read('cross_precision.csv')
    fig,axs=plt.subplots(2,2,figsize=(6.85,5.2))
    fig.subplots_adjust(left=.14,right=.98,bottom=.20,top=.91,hspace=.55,wspace=.45)
    methods=[('dtangle_tmm','dtangle + TMM',style.BLUE,'o'),
             ('dtangle_marker','dtangle + MRN',style.TEAL,'s'),
             ('music_marker','MuSiC, fixed markers',style.ROSE,'^'),
             ('linear_svr_marker','Linear SVR, fixed markers',style.VIOLET,'D')]
    for row,(metric,effect,label) in enumerate([('fdp',0.,'FDR change (pp)'),('power',.5,'Power change (pp)')]):
        for col,probe in enumerate(['test_counts','marker_counts']):
            ax=axs[row,col]
            for j,(method,name,color,marker) in enumerate(methods):
                g=[select(f,reference_source=r,target_source=t,effect=effect,method=method,probe=probe,metric=metric) for r,t in ROUTES]
                ax.errorbar(np.arange(4)+(j-1.5)*.15,[100*z.donor_equal_change for z in g],
                    yerr=[100*z.half_width for z in g],fmt=marker,color=color,ms=4,capsize=2,
                    elinewidth=.8,linestyle='none',label=name)
            ax.axhline(0,color='#777777',ls=':',lw=.8)
            ax.set_xticks(range(4),SHORT,fontsize=8);ax.set_ylabel(label)
            if row==0:ax.set_title('Tested-count probe' if col==0 else 'Marker-count probe')
    for ax,label in zip(axs.flat,'abcd'):core.letter(ax,label)
    fig.legend(*axs[0,0].get_legend_handles_labels(),loc='lower center',ncol=2,frameon=False,fontsize=8)
    finish(fig,'FigS9')


def figureS10():
    f=read('clinical_estimator_contrasts.csv'); design=read('clinical_estimators.csv')
    methods=['dtangle_tmm','dtangle_marker','music_marker','music_broad','music_non_test','linear_svr_marker']
    labels=['dtangle + TMM','dtangle + MRN','MuSiC markers','MuSiC broad','MuSiC non-test','Linear SVR']
    fig,axs=plt.subplots(2,3,figsize=(6.85,5.8),gridspec_kw={'width_ratios':[1,1,.55]})
    fig.subplots_adjust(left=.205,right=.982,bottom=.16,top=.91,hspace=.68,wspace=.34)
    styles=[('LIBD','test_counts',style.INK,'o'),('LIBD','marker_counts',style.INK,'^'),
            ('SEAAD','test_counts',style.ROSE,'o'),('SEAAD','marker_counts',style.ROSE,'^')]
    cols=[('beta_rms','Coefficient-change RMS',1),('crossing_fraction','BH crossings (%)',100)]
    for row,cohort in enumerate(['GSE125583','GSE53697']):
        for col,(metric,label,factor) in enumerate(cols):
            ax=axs[row,col]
            for j,(ref,probe,color,marker) in enumerate(styles):
                for k,method in enumerate(methods):
                    r=select(f,cohort=cohort,reference=ref,method=method,probe=probe)
                    value=r[metric]
                    if pd.isna(value):
                        if str(r.available).lower()=='true':raise ValueError('Estimable clinical contrast has no value')
                        continue
                    if str(r.available).lower()!='true':raise ValueError('Unavailable clinical contrast has a value')
                    ax.scatter(factor*value,k+(j-1.5)*.16,color=color,marker=marker,s=18,zorder=3)
            ax.set_yticks(range(6),labels if col==0 else []);ax.set_ylim(5.6,-.6)
            all_values=f.loc[f.probe.isin(['marker_counts','test_counts']) & f[metric].notna(),metric]*factor
            upper=max(float(all_values.max()),1e-8);ax.set_xlim(-.045*upper,1.07*upper)
            ax.xaxis.set_major_locator(MaxNLocator(3));ax.set_xlabel(label)
            ax.axvline(0,color=style.MUTED,lw=.7)
            ax.set_title(cohort if col==0 else ('Test and marker probes' if row==0 else ''))
        ax=axs[row,2];ax.set_xlim(-.52,1.52);ax.set_ylim(5.6,-.6)
        ax.set_yticks([]);ax.set_xticks([0,1],['L','S']);ax.xaxis.tick_top();ax.tick_params(length=0)
        ax.set_title('Fit status',pad=10)
        for side in ax.spines.values():side.set_visible(False)
        for k,method in enumerate(methods):
            for j,ref in enumerate(('LIBD','SEAAD')):
                marks=[]
                for probe in ('test_counts','marker_counts'):
                    r=select(f,cohort=cohort,reference=ref,method=method,probe=probe)
                    dr=select(design,cohort=cohort,reference=ref,method=method,probe=probe)
                    reason='' if pd.isna(dr.reason) else str(dr.reason)
                    if str(r.available).lower()=='true':mark='OK'
                    elif 'MARKER_COVERAGE' in (reason+' '+str(r.reason)).upper():mark='M'
                    elif 'HC3' in (reason+' '+str(r.reason)).upper():mark='H'
                    else:mark='NA'
                    marks.append(mark)
                status=marks[0] if marks[0]==marks[1] else '/'.join(marks)
                ax.text(j,k,status,ha='center',va='center',fontsize=8,color=style.MUTED if status=='OK' else '#95422E',
                        fontweight='normal' if status=='OK' else 'bold')
    for ax,l in zip(axs.flat,'abcdef'):core.letter(ax,l)
    handles=[Line2D([],[],marker=mk,ls='none',color=c,ms=5,label=ref.replace('SEAAD','SEA-AD')+
                    (' / tested' if p=='test_counts' else ' / marker')) for ref,p,c,mk in styles]
    fig.legend(handles=handles,ncol=4,loc='lower center',bbox_to_anchor=(.53,.045),columnspacing=.9,handlelength=1)
    fig.text(.53,.025,'Fit status: L = LIBD, S = SEA-AD; M = marker coverage, H = HC3 non-estimable.',
             ha='center',fontsize=8)
    finish(fig,'FigS10')


def figureS11():
    data=fd.imaging_rows(read('imaging_comparison.csv'))
    classes=['Astro','EndoMural','Excit','Inhib','Micro','OligoOPC']
    if set(data.index)!=set(classes):raise ValueError('Incomplete imaging cell classes')
    data=data.loc[classes]
    labels=['Astrocytes','Endothelial / mural','Excitatory neurons','Inhibitory neurons','Microglia','Oligodendrocytes / OPC']
    fig=plt.figure(figsize=(6.85,3.85))
    gs=fig.add_gridspec(1,2,left=.255,right=.965,bottom=.23,top=.85,width_ratios=[3.45,1],wspace=.10)
    ax=fig.add_subplot(gs[0]);sizes=fig.add_subplot(gs[1])
    y=np.arange(len(data))
    for j,(ref,color,marker) in enumerate([('LIBD',style.INK,'o'),('SEAAD','#A45C6F','s')]):
        ax.plot(data[ref],y+(j-.5)*.31,marker=marker,ls='none',color=color,mfc=color if j==0 else 'white',
                ms=6,mew=1.3,label=ref.replace('SEAAD','SEA-AD'),zorder=4)
    bound=max(1, float(np.ceil(4*np.abs(data[['LIBD','SEAAD']].to_numpy()).max())/4))
    ax.set_xlim(-1.08*bound,1.08*bound);ax.set_ylim(5.55,-.55);ax.set_yticks(y,labels)
    ax.axvline(0,c=style.INK,lw=.9);ax.set_xlabel('Difference in mean absolute error (percentage points)')
    ax.text(.02,1.055,'Lower with MRN',transform=ax.transAxes,color=style.MUTED,fontsize=8)
    ax.text(.98,1.055,'Higher with MRN',ha='right',transform=ax.transAxes,color=style.MUTED,fontsize=8)
    sizes.set(xlim=(-.3,1.25),ylim=(5.55,-.55));sizes.axis('off')
    for x,word,col in [(0,'Donors','donors'),(1,'Blocks','blocks')]:
        sizes.text(x,1.055,word,ha='center',transform=sizes.get_xaxis_transform(),fontsize=8,weight='bold')
        for iy,n in enumerate(data[col]):
            sizes.text(x,iy,str(int(n)),ha='center',va='center',fontsize=9)
    fig.legend(*ax.get_legend_handles_labels(),loc='lower center',ncol=2,bbox_to_anchor=(.54,.01))
    finish(fig,'FigS11')


def main():
    global DATA,OUT,FORMATS,SOFTWARE_ACCEPTANCE,VIEWS,INDEX
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data',type=Path,default=HERE/'data')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--formats',nargs='+',choices=['pdf','png','eps','tiff','svg'],default=['pdf'])
    p.add_argument('--software-acceptance',action='store_true')
    p.add_argument('--estimators-only',action='store_true')
    args=p.parse_args();DATA=args.data;OUT=separate(args.output,DATA);FORMATS=tuple(dict.fromkeys(args.formats));OUT.mkdir(parents=True,exist_ok=True)
    SOFTWARE_ACCEPTANCE=args.software_acceptance
    core.finish=finish;core.read=read
    functions=[core.figure1,core.figure2,core.figure3,core.figure4,figure5,figure6,core.figure7,figure8,figure9,
               core.figureS1,figureS2,figureS3,figureS4,core.figureS5,core.figureS6,figureS7,figureS8]
    MANIFEST.clear()
    VIEWS = json.loads((DATA/'catalog.json').read_text(encoding='utf-8'))['views']
    INDEX = figure_index(VIEWS,args.estimators_only)
    if 'clinical_estimators.csv' in VIEWS:
        functions += [figureS9,figureS10]
    if 'imaging_comparison.csv' in VIEWS:
        functions.append(figureS11)
    if args.estimators_only:
        functions=[figureS8,figureS9,figureS10]
    for fn in functions:
        READS.clear()
        fn();print('Rendered',len(MANIFEST),'of',len(functions),flush=True)
    if list(MANIFEST)!=list(INDEX):
        raise ValueError('Unexpected figure set for the supplied numerical inputs')
    (OUT/'figure_sources.json').write_text(json.dumps(MANIFEST,ensure_ascii=False,indent=2),encoding='utf-8')


if __name__=='__main__':
    main()
