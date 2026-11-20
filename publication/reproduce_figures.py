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
    points = fd.subset_points(read('common_support_donor_subsets.csv'))
    # One large contrast, with the other three settings on the same FDR scale.
    fig = plt.figure(figsize=(6.85, 5.6))
    top = fig.add_gridspec(1, 2, left=.115, right=.965, bottom=.565, top=.91,
                           width_ratios=[1.28, 1], wspace=.29)
    bottom = fig.add_gridspec(1, 3, left=.115, right=.965, bottom=.175, top=.397, wspace=.35)
    hero, ledger = fig.add_subplot(top[0]), fig.add_subplot(top[1])
    route_order = [ROUTES[2], ROUTES[0], ROUTES[1], ROUTES[3]]
    shape = [('zero', 'Unperturbed', style.GRAY, 'x'), ('matched_', 'Matched', style.TEAL, 's'),
             ('random_', 'Random', style.VIOLET, 'o'), ('legacy','Predefined', '#96472D', '*'),
             ('random_mean','Random mean', '#664784','D')]
    maximum = max(10, int(np.ceil(100 * points.fdp.max() / 10) * 10))
    for i, (ax, (ref,target)) in enumerate(zip([hero] + [fig.add_subplot(c) for c in bottom], route_order)):
        g = points[points.reference_source.eq(ref) & points.target_source.eq(target)]
        for prefix, label, color, marker in shape:
            mask = g.base_subset.str.startswith(prefix) if prefix.endswith('_') else g.base_subset.eq(prefix)
            if prefix == 'random_': mask &= ~g.base_subset.eq('random_mean')
            v = g[mask]
            ax.plot(v.proportion_rmse, 100*v.fdp, ls='none', marker=marker, mec=color,
                    mfc='white' if prefix=='random_mean' else color, mew=1.3 if prefix=='random_mean' else .8,
                    ms={'legacy':9, 'random_mean':7}.get(prefix,4.6), zorder=5, label=label)
        ax.axhline(5,color=style.MUTED,lw=.8,ls=':')
        ax.set_ylim(-1, maximum*1.04); ax.yaxis.set_major_locator(MaxNLocator(4))
        ax.xaxis.set_major_locator(MaxNLocator(3 if i==0 else 2, prune='both'))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.4f')); ax.margins(x=.09)
        if i in (0,1): ax.set_ylabel('Global-null FDR (%)')
        if i not in (0,1): ax.tick_params(axis='y',labelleft=False)
        ax.set_xlabel('Proportion RMSE')
        style.route_title(ax,ref,target); core.letter(ax,'abcd'[i])
        if i==0: focus=g.set_index('base_subset')
    ledger.axis('off'); ledger.set(xlim=(0,1),ylim=(0,1))
    for x, word, alignment in [(0,'Subset','left'),(.69,'RMSE','right'),(1,'FDR','right')]:
        ledger.text(x,.94,word,ha=alignment,va='top',weight='bold',fontsize=8)
    ledger.plot([0,1],[.80,.80],c=style.INK,lw=.7)
    for y,key,label,color in [(.68,'legacy','Predefined','#96472D'),(.47,'random_mean','Random mean','#664784')]:
        row=focus.loc[key]
        ledger.text(0,y,label,ha='left',va='center',fontsize=8,color=color,weight='bold')
        ledger.text(.69,y,f'{row.proportion_rmse:.4f}',ha='right',va='center',fontsize=9)
        ledger.text(1,y,f'{100*row.fdp:.2f}%',ha='right',va='center',fontsize=9)
    ledger.plot([0,1],[.32,.32],c=style.GRID,lw=.9)
    random=focus[focus.index.str.match(r'random_\d+$')].sort_index()
    ledger.text(0,.19,'Individual random-subset FDR',fontsize=8,va='center')
    ledger.text(0,.035,'  /  '.join(f'{v:.2f}%' for v in 100*random.fdp),fontsize=8,va='center',color='#664784')
    fig.legend(*hero.get_legend_handles_labels(),ncol=5,loc='lower center',bbox_to_anchor=(.5,.012),
               columnspacing=1.05,handlelength=1.1)
    finish(fig,'Fig6')


def figure8():
    f=read('clinical_coefficient_pairs.csv.gz'); s=read('clinical_stability.csv')
    fig=plt.figure(figsize=(6.85,5.85))
    grid=fig.add_gridspec(2,2,left=.115,right=.97,bottom=.12,top=.925,
                          height_ratios=[1.55,1],wspace=.43,hspace=.59)
    axs=[fig.add_subplot(cell) for cell in grid]
    refs=['LIBD','SEAAD']; counts={}
    keys={'Neither':('#D6DDE1','.'),
          'Both':('#68747C','D'),'TMM only':(style.BLUE,'^'),'MRN only':(style.TEAL,'s')}
    for ax,ref in zip(axs[:2],refs):
        g=f[f.reference.eq(ref)]
        cats,counts[ref]=fd.clinical_status(g,select(s,reference=ref,method='marker_matched'))
        for i,(key,(color,marker)) in enumerate(keys.items()):
            v=g[cats==key]
            ax.scatter(v.tmm_beta,v.mrn_beta,s=[3.2,2.9,4.7,6.8][i],marker=marker,
                       c=color,edgecolors='none',alpha=[.70,.66,.87,.96][i],rasterized=True,zorder=2+i)
        low=g[['tmm_beta','mrn_beta']].min().min(); high=g[['tmm_beta','mrn_beta']].max().max()
        pad=.04*(high-low); lim=[low-pad,high+pad]
        ax.plot(lim,lim,color=style.MUTED,ls=':',lw=.8,zorder=1)
        ax.set(xlim=lim,ylim=lim,xlabel='TMM disease coefficient',ylabel='MRN disease coefficient',
               title=ref.replace('SEAAD','SEA-AD')+' reference',aspect='equal')
    ax=axs[2];ax.axis('off');ax.set(xlim=(0,1),ylim=(0,1))
    for x,label,alignment in [(0,'BH status','left'),(.67,'LIBD','right'),(1,'SEA-AD','right')]:
        ax.text(x,.965,label,ha=alignment,va='top',weight='bold',fontsize=8)
    ax.plot([0,1],[.81,.81],c=style.INK,lw=.7)
    for y,status in zip([.68,.49,.30,.11],fd.CLASSES):
        color,marker=keys[status]
        ax.plot(.028,y,marker=marker,c=color,ls='none',ms=6)
        ax.text(.11,y,status,va='center',fontsize=8)
        for x,ref in zip([.67,1],refs):
            ax.text(x,y,f'{counts[ref][status]:,}',ha='right',va='center',fontsize=9,
                    weight='bold' if 'only' in status else 'normal')
    ax=axs[3]
    for ref,c,ls in [('LIBD',style.INK,'-'),('SEAAD',style.ROSE,'--')]:
        g=f[f.reference.eq(ref)];ratio=np.sort((g.mrn_se/g.tmm_se).to_numpy())
        ax.plot(ratio,np.arange(1,len(ratio)+1)/len(ratio),color=c,ls=ls,label=ref.replace('SEAAD','SEA-AD'))
    ax.axvline(1,color=style.MUTED,ls=':',lw=.8);ax.set_xlabel('MRN / TMM standard error');ax.set_ylabel('Fraction of genes')
    ax.set_xscale('log');ax.set_xticks([1,1.2,1.4,1.6,1.8]);ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.xaxis.set_minor_locator(NullLocator());ax.legend(loc='lower right')
    for ax,l in zip(axs,'abcd'):core.letter(ax,l)
    finish(fig,'Fig8')


def figure9():
    f=read('clinical_finite.csv')
    measured='clinical_probes.csv' in VIEWS
    fig,axs=plt.subplots(2,2,figsize=(6.85,6.35))
    fig.subplots_adjust(left=.135,right=.968,bottom=.185 if measured else .12,top=.925,hspace=.92,wspace=.54)
    dose=read('clinical_probes.csv') if measured else None
    if measured:
        dose=fd.fusiform_doses(f,dose)
    refs=['LIBD','SEAAD']
    arms=[('tmm','o','TMM','-'),('tmm_exclude_factor','^','Factor-only','--'),
          ('tmm_exclude_both','s','Three invariant arms',':')]
    for j,ref in enumerate(refs):
        for i,metric in enumerate(['threshold_crossing_fraction','coefficient_change_rms'][:1 if measured else 2]):
            ax=axs[i,j]; sub=f[f.reference.eq(ref)]
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
    fig.legend(*axs[0,0].get_legend_handles_labels(),ncol=3,loc='center',bbox_to_anchor=(.55,.512),
               handlelength=1.8,columnspacing=1.1)
    if measured:
        ax=axs[1,0]
        for k,level in enumerate(('low','middle','high')):
            for i,ref in enumerate(refs):
                for direction,mk,offset in [('up','^',-.065),('down','v',.065)]:
                    row=select(dose,probe=f'single_{level}_{direction}')
                    value=row[ref+'_singleton_beta_change']*1e5
                    if not np.isfinite(value): raise ValueError('Missing single-gene coefficient change')
                    ax.scatter(k+(i-.5)*.36+offset,value,marker=mk,facecolors=style.INK if i==0 else 'white',
                               s=35,edgecolors=style.INK if i==0 else '#994F63',linewidths=1,zorder=3)
        ax.axhline(0,color=style.MUTED,lw=.8);ax.set_xlim(-.55,2.55)
        ax.set_xticks(range(3),['Low','Middle','High']); ax.set_xlabel('Single gene by count share')
        ax.set_ylabel(r'Own coefficient change ($10^{-5}$ units)')
        own=dose[dose.probe.str.startswith('single_')][[ref+'_own_BH_crossings' for ref in refs]]
        if own.shape!=(6,2) or not own.isin([0,1]).all().all():
            raise ValueError('Expected six paired singleton interventions under both references')
        ax.set_title(f'Own BH status changed: {int(own.to_numpy().sum())} of {own.size} fits')
        h=[Line2D([],[],marker=m,color=style.MUTED,ls='none',ms=5,label=t)
           for m,t in [('^','Increase'),('v','Decrease')]]
        ax.legend(handles=h,loc='lower center',ncol=2,bbox_to_anchor=(.5,-.51),columnspacing=.8)
        ax=axs[1,1]
        for j,(group,color,mk) in enumerate([('AD',style.INK,'o'),('control',style.ROSE,'s')]):
            for k,direction in enumerate(('down','up')):
                row=select(dose,probe='test_family_'+direction)
                q1,med,q3=[row[f'{group}_effective_{v}_pct'] for v in ('q25','median','q75')]
                ax.errorbar(k+(j-.5)*.18,med,yerr=[[med-q1],[q3-med]],color=color,fmt=mk,ms=5,
                            capsize=3,label=('AD' if j==0 else 'Control') if not k else None)
                if direction=='up' and j==0:
                    ax.annotate(f'{med:+.2f}%',(k-.09,q3),xytext=(0,4),textcoords='offset points',
                                ha='center',fontsize=8,weight='bold')
        ticks=[select(dose,probe='test_family_'+t).AD_input_change_pct for t in ('down','up')]
        ax.set_xticks([0,1],[f'{x:+.2f}%' for x in ticks]);ax.set_xlim(-.55,1.55)
        yl=ax.get_ylim();ax.set_ylim(yl[0],yl[1]+.16*(yl[1]-yl[0]))
        ax.axhline(0,color=style.MUTED,lw=.8);ax.set_xlabel('Full-family count perturbation in AD')
        ax.set_ylabel('Effective TMM denominator Δ (%)')
        share=select(dose,probe='test_family_up').AD_count_share_pct
        ax.set_title(f'Full family: {share:.2f}% of AD counts')
        ax.legend(loc='lower center',ncol=2,bbox_to_anchor=(.5,-.51),columnspacing=.8)
        # References and diagnostic groups are different units, even when both use two symbols.
        ax_l=axs[1,0]
        ax_l.text(.5,-.57,'LIBD: filled   /   SEA-AD: open',ha='center',transform=ax_l.transAxes,fontsize=8)
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
