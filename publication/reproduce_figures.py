"""Reproduce all manuscript figures from the supplied numerical source tables."""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline.publication import load_frame
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

HERE = Path(__file__).resolve().parent
CAPTIONS = json.loads((HERE/'figure_captions.json').read_text(encoding='utf-8'))
COL = style.METHOD_COLORS
METHODS = ['tmm','tmm_exclude_factor','tmm_exclude_both','marker_denominator_only','marker_matched']
LABELS = ['TMM','Factor-only','Complete\nexclusion','Denom.','MRN']
ROUTES, SHORT = core.ROUTES, core.SHORT
MANIFEST = {}


def read(name):
    return load_frame(DATA, name)


def finish(fig, name, sources):
    style.clean_axes(fig)
    style.check_text_bounds(fig)
    for fmt in FORMATS:
        kwargs = {'dpi':300 if fmt=='png' else 600, 'facecolor':'white'}
        if fmt=='tiff':
            kwargs['pil_kwargs']={'compression':'tiff_lzw'}
        fig.savefig(OUT/(name+'.'+fmt), **kwargs)
        if fmt=='tiff':
            target=OUT/(name+'.'+fmt)
            with Image.open(target) as source:
                rgb=source.convert('RGB')
            rgb.save(target,compression='tiff_lzw',dpi=(600,600))
    catalog=json.loads((DATA/'catalog.json').read_text(encoding='utf-8'))['views']
    MANIFEST[name]={'caption':CAPTIONS[name], 'sources':sorted({catalog[s]['file'] for s in sources}),
                    'views':{s:catalog[s] for s in sources},
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
    finish(fig,'Fig5',['common_support_summary.csv'])


def figure6():
    f=read('common_support_donor_subsets.csv')
    f=f[f.method.eq('marker_matched')&f.outcome.eq('null')&
        (f.support_rule.eq('positive_only')|f.base_subset.eq('zero'))]
    f=f.groupby(['reference_source','target_source','base_subset'],as_index=False)[['fdp','proportion_rmse']].mean()
    fig,axs=plt.subplots(2,2,figsize=(6.85,5.05))
    fig.subplots_adjust(left=.11,right=.98,bottom=.18,top=.91,hspace=.70,wspace=.40)
    styles=[('zero','Unperturbed',style.GRAY,'x'),('legacy','Predefined',style.ROSE,'*'),
            ('matched','Matched',style.TEAL,'D'),('random','Random',style.VIOLET,'o')]
    for ax,route,tag,letter in zip(axs.flat,ROUTES,SHORT,'abcd'):
        g=f[f.reference_source.eq(route[0])&f.target_source.eq(route[1])]
        for prefix,label,color,marker in styles:
            d=g[g.base_subset.str.startswith(prefix)]
            ax.scatter(d.proportion_rmse,100*d.fdp,label=label,c=color,marker=marker,
                       s=80 if prefix=='legacy' else 30,zorder=4)
        ax.set_xlabel('Proportion RMSE'); ax.set_ylabel('Global-null FDR (%)')
        ax.set_ylim(0,min(100,max(45,110*f.fdp.max()))); ax.axhline(5,color='#777777',ls=':',lw=.8)
        ax.xaxis.set_major_locator(MaxNLocator(4));ax.xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
        style.route_title(ax,*route);core.letter(ax,letter)
    fig.legend(*axs[0,0].get_legend_handles_labels(),ncol=4,loc='lower center',frameon=False)
    finish(fig,'Fig6',['common_support_donor_subsets.csv'])


def figure8():
    f=read('clinical_coefficient_pairs.csv.gz'); s=read('clinical_stability.csv')
    fig,axs=plt.subplots(2,2,figsize=(6.85,5.8))
    fig.subplots_adjust(left=.13,right=.975,bottom=.13,top=.92,hspace=.72,wspace=.52)
    for ax,ref in zip(axs[0],['LIBD','SEAAD']):
        g=f[f.reference.eq(ref)]
        cats=np.where(g.tmm_q.le(.05)&g.mrn_q.le(.05),'Both',np.where(g.tmm_q.le(.05),'TMM only',np.where(g.mrn_q.le(.05),'MRN only','Neither')))
        for cat,c,marker in [('Neither','#DCE1E5','o'),('TMM only',COL['tmm'],'^'),('MRN only',COL['marker_matched'],'s'),('Both','#525B65','D')]:
            d=g[cats==cat]
            ax.scatter(d.tmm_beta,d.mrn_beta,s=2.1,c=c,marker=marker,rasterized=True,edgecolors='none',label=cat)
        limits=[g[['tmm_beta','mrn_beta']].min().min(),g[['tmm_beta','mrn_beta']].max().max()]
        span=limits[1]-limits[0];lim=[limits[0]-.04*span,limits[1]+.04*span]
        ax.plot(lim,lim,color='#777777',ls=':',lw=.8);ax.set_xlim(lim);ax.set_ylim(lim)
        ax.set_xlabel('TMM disease coefficient');ax.set_ylabel('MRN disease coefficient')
        ax.set_title(ref.replace('SEAAD','SEA-AD')+' reference')
        ax.set_aspect('equal',adjustable='box')
    ax=axs[1,0]
    for j,ref in enumerate(['LIBD','SEAAD']):
        values=[select(s,reference=ref,method=m).discoveries for m in METHODS]
        y=np.arange(5)+(j-.5)*.20
        ax.plot(values,y,marker='o' if j==0 else 's',color=style.REFERENCE_COLORS[ref],
                ls='none',label=ref.replace('SEAAD','SEA-AD'))
    ax.set_yticks(range(5),['TMM','Factor-only','Complete','Denom.','MRN'])
    ax.invert_yaxis();ax.set_xlim(0,3900);ax.set_xticks([0,1000,2000,3000])
    ax.set_xlabel('Discoveries (q ≤ 0.05)')
    ax.legend(frameon=False,fontsize=8,loc='lower left',bbox_to_anchor=(-.08,-.43),ncol=2)
    ax=axs[1,1]
    for ref,color,ls in [('LIBD',style.BLUE,'-'),('SEAAD',style.ROSE,'--')]:
        g=f[f.reference.eq(ref)];ratio=np.sort((g.mrn_se/g.tmm_se).to_numpy())
        ax.plot(ratio,np.arange(1,len(ratio)+1)/len(ratio),color=color,ls=ls,label=ref.replace('SEAAD','SEA-AD'))
    ax.axvline(1,color='#777777',ls=':',lw=.8);ax.set_xlabel('MRN / TMM standard error');ax.set_ylabel('Fraction of genes')
    ax.set_xscale('log');ax.set_xticks([1,1.2,1.4,1.6,1.8]);ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.xaxis.set_minor_locator(NullLocator());ax.legend(frameon=False,loc='lower right',fontsize=7)
    for ax,l in zip(axs.flat,'abcd'):core.letter(ax,l)
    fig.legend(*axs[0,0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.55,.535),
               ncol=4,frameon=False,markerscale=3)
    finish(fig,'Fig8',['clinical_coefficient_pairs.csv.gz','clinical_stability.csv'])


def figure9():
    f=read('clinical_finite.csv')
    fig,axs=plt.subplots(2,2,figsize=(6.85,5.15))
    fig.subplots_adjust(left=.12,right=.975,bottom=.23,top=.91,hspace=.65,wspace=.44)
    for i,ref in enumerate(['LIBD','SEAAD']):
        for j,metric in enumerate(['threshold_crossing_fraction','coefficient_change_rms']):
            ax=axs[j,i]
            for m,marker,label in [('tmm','o','TMM'),('tmm_exclude_factor','^','Factor-only exclusion'),
                                   ('tmm_exclude_both','s','Complete / denominator-only / MRN')]:
                g=f[f.reference.eq(ref)&f.method.eq(m)].sort_values('log_count_step')
                factor=100 if j==0 else 1
                ax.plot(100*np.expm1(g.log_count_step),factor*g[metric],marker=marker,ms=4,color=COL[m],label=label)
            ax.set_xlabel('Tested-count increase (%)');ax.set_xticks([0,5,10])
            ax.set_ylabel('Threshold crossings (%)' if j==0 else 'Coefficient-change RMS')
            ax.set_ylim(bottom=-.7 if j==0 else -.002)
            if j==0:ax.set_ylim(top=26)
            if j==0:ax.set_title(ref.replace('SEAAD','SEA-AD')+' reference')
    for ax,l in zip(axs.flat,'abcd'):core.letter(ax,l)
    fig.legend(*axs[0,0].get_legend_handles_labels(),ncol=1,loc='lower center',frameon=False,fontsize=8)
    finish(fig,'Fig9',['clinical_finite.csv'])


def figureS2():
    f=read('dose_retention.csv');c=read('dose_by_class.csv')
    fig,axs=plt.subplots(2,1,figsize=(6.85,4.5),gridspec_kw={'height_ratios':[1,1.5]})
    fig.subplots_adjust(left=.12,right=.9,bottom=.12,top=.93,hspace=.62)
    vals=[100*select(f,reference_source=r,target_source=t).retained_fraction for r,t in ROUTES]
    axs[0].bar(range(4),vals,color=[style.BLUE,style.ROSE,style.BLUE,style.ROSE],width=.52)
    axs[0].set_xticks(range(4),SHORT);axs[0].set_ylim(0,110);axs[0].set_ylabel('Retained counts (%)')
    for i,y in enumerate(vals):axs[0].text(i,y+2,f'{y:.2f}',ha='center',fontsize=8)
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
    finish(fig,'FigS2',['dose_retention.csv','dose_by_class.csv'])


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
    finish(fig,'FigS3',['common_support_contrasts.csv'])


def figureS4():
    f=read('clinical_design.csv');fig,axs=plt.subplots(1,2,figsize=(6.85,3.25))
    fig.subplots_adjust(left=.11,right=.97,bottom=.32,top=.9,wspace=.4)
    for ax,metric,yl,l in zip(axs,['exposure_r2','exposure_residual_ss'],['Disease–nuisance R²','Residual exposure sum of squares'],'ab'):
        for ref,c,marker in [('LIBD',style.BLUE,'o'),('SEAAD',style.ROSE,'s')]:
            vals=[select(f,reference=ref,method=m)[metric] for m in METHODS]
            ax.plot(range(5),vals,color=c,marker=marker,label=ref.replace('SEAAD','SEA-AD'))
        ax.set_xticks(range(5),LABELS,rotation=35,ha='right',fontsize=8);ax.set_ylabel(yl);core.letter(ax,l)
    fig.legend(*axs[0].get_legend_handles_labels(),loc='lower center',ncol=2,frameon=False)
    finish(fig,'FigS4',['clinical_design.csv'])


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
    finish(fig,'FigS7',['diagnostic_qualification.csv','clinical_derivatives.csv'])


def main():
    global DATA,OUT,FORMATS
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data',type=Path,default=HERE/'data')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--formats',nargs='+',choices=['pdf','png','eps','tiff'],default=['pdf'])
    args=p.parse_args();DATA=args.data;OUT=separate(args.output,DATA);FORMATS=tuple(dict.fromkeys(args.formats));OUT.mkdir(parents=True,exist_ok=True)
    core.DATA=DATA;core.OUT=OUT;core.finish=finish;core.read=read
    for fn in [core.figure1,core.figure2,core.figure3,core.figure4,figure5,figure6,core.figure7,figure8,figure9,
               core.figureS1,figureS2,figureS3,figureS4,core.figureS5,core.figureS6,figureS7]:
        fn();print('Rendered',len(MANIFEST),'of 16',flush=True)
    assert set(MANIFEST)==set(CAPTIONS)
    (OUT/'figure_sources.json').write_text(json.dumps(MANIFEST,ensure_ascii=False,indent=2),encoding='utf-8')


if __name__=='__main__':
    main()
