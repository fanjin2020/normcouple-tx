"""Shared figure panels; the renderer supplies data and export callbacks."""

from pathlib import Path

import argparse

import json

import string

import numpy as np

import pandas as pd

import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt

from matplotlib import colors

from matplotlib.patches import FancyBboxPatch

from matplotlib.ticker import MaxNLocator, FormatStrFormatter, NullLocator

import figure_style as style

plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.labelsize':9,'xtick.labelsize':8,
                     'ytick.labelsize':8,'legend.fontsize':8,'axes.spines.top':False,'axes.spines.right':False,
                     'axes.linewidth':.7,'lines.linewidth':1.2,'pdf.fonttype':42,'ps.fonttype':42,
                     'svg.fonttype':'none','savefig.facecolor':'white','axes.unicode_minus':True})

style.apply()
COL=style.METHOD_COLORS

LAB={'tmm':'TMM','marker_matched':'MRN','marker_denominator_only':'Denominator only','oracle':'Oracle'}

ROUTES=[('LIBD','LIBD'),('SEAAD','SEAAD'),('LIBD','GSE174367'),('SEAAD','GSE174367')]

SHORT=['L-in','S-in','L-ext','S-ext']

CAPTIONS={}

FIGURES={}

FORMATS=('pdf','png')

CAPTION_TEXT=json.loads((Path(__file__).parent/'figure_captions.json').read_text(encoding='utf-8'))

def read(name):
    return pd.read_csv(DATA/name,keep_default_na=False,na_values=['','NaN','nan'])

def letter(ax,label):
    style.panel(ax,label)

def choose(f,**kwargs):
    part=f
    for key,value in kwargs.items(): part=part[part[key].eq(value)]
    if len(part)!=1: raise ValueError(('Expected one data row',kwargs,len(part)))
    return part.iloc[0]

def routes_axes(ax):
    ax.set_xticks(range(4),SHORT)
    ax.set_xlim(-.5,3.5)

def dots(ax,frame,metric,methods,mean='mean',low='low',high='high',factor=1.):
    for j,method in enumerate(methods):
        ys=[];lo=[];hi=[]
        for a,b in ROUTES:
            r=choose(frame,reference_source=a,target_source=b,method=method,metric=metric)
            ys.append(r[mean]*factor);lo.append((r[mean]-r[low])*factor);hi.append((r[high]-r[mean])*factor)
        x=np.arange(4)+(j-(len(methods)-1)/2)*.15
        ax.errorbar(x,ys,yerr=[lo,hi],color=COL.get(method,'#555555'),marker=['o','s','^','D'][j],
                    ms=4,capsize=2,linestyle='none',label=LAB.get(method,method))
    routes_axes(ax)

def figure1():
    fig=plt.figure(figsize=(6.85,4.35))
    gs=fig.add_gridspec(2,2,height_ratios=[.9,1.4],hspace=.58,wspace=.48,left=.11,right=.97,top=.94,bottom=.14)
    a=fig.add_subplot(gs[0,:]);a.set_axis_off()
    a.text(-.03,1.06,'a',transform=a.transAxes,fontweight='bold',fontsize=11,va='bottom')
    xs=np.linspace(.09,.91,5)
    labels=['Non-marker\ncounts','Shared\ndenominator','Normalized\nmarker input','Estimated\ncomposition','Adjusted\ncoefficient']
    fills=['#F4F5F7','#E8F1F5','#ECF4F3','#F2EDF5','#F8ECEF']
    for x,label,fill in zip(xs,labels,fills):
        box=FancyBboxPatch((x-.077,.44),.154,.40,boxstyle='round,pad=0.009,rounding_size=0.018',
                            transform=a.transAxes,facecolor=fill,edgecolor='#75818B',lw=.6)
        a.add_patch(box);a.text(x,.64,label,transform=a.transAxes,ha='center',va='center',fontsize=8)
    for x0,x1 in zip(xs[:-1],xs[1:]):
        a.annotate('',xy=(x1-.083,.64),xytext=(x0+.083,.64),xycoords='axes fraction',
                   arrowprops=dict(arrowstyle='->',color='#444444',lw=1))
    a.annotate('Marker-restricted denominator',xy=(xs[1],.40),xytext=(xs[1],.08),xycoords='axes fraction',
               ha='center',fontsize=8,color=style.INK,
               arrowprops=dict(arrowstyle='-|>',color=style.TEAL,lw=1.2))
    b=fig.add_subplot(gs[1,0]);letter(b,'b')
    curves=read('analytic_curves.csv')
    signal=curves.signal_over_pseudocount.to_numpy()
    b.plot(signal,curves.derivative_log_denominator,color=COL['tmm'])
    b.set_xscale('log');b.set_xlabel('Marker signal / pseudocount');b.set_ylabel('dZ / d ln D')
    b.axhline(-1/np.log(2),ls=':',color='#777777',lw=.7)
    c=fig.add_subplot(gs[1,1]);letter(c,'c')
    for shift,color,ls in zip([.1,.25,.5],[style.BLUE,style.TEAL,style.ROSE],['-','--','-.']):
        change=curves['change_'+str(shift)]
        c.plot(signal,change,color=color,ls=ls,label=f'{shift:.0%}')
    c.set_xscale('log');c.set_xlabel('Marker signal / pseudocount');c.set_ylabel('Marker-input change (log2)')
    c.legend(title='Denominator increase',title_fontsize=8,loc='lower left',frameon=False)
    finish(fig, 'Fig1', ['analytic_curves.csv'])

def figure2():
    f=read('normalization.csv')
    fig,axs=plt.subplots(1,3,figsize=(6.85,3.25))
    fig.subplots_adjust(left=.085,right=.985,bottom=.29,top=.88,wspace=.58)
    for ax,l,metric,label,factor in zip(axs,'abc',['fdp','null_coverage','proportion_rmse'],
                                      ['FDR (%)','Unaffected-gene coverage (%)','Proportion RMSE'],[100,100,1]):
        letter(ax,l);dots(ax,f,metric,['Observed total','Fixed total','TMM','DESeq2 poscounts'],factor=factor)
        ax.set_ylabel(label)
        if metric=='fdp': ax.axhline(5,color='#555555',ls=':',lw=.7);ax.set_ylim(bottom=0)
        if metric=='null_coverage':ax.axhline(95,color='#555555',ls=':',lw=.7)
    handles,labels=axs[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',ncol=2,frameon=False,bbox_to_anchor=(.5,.015))
    finish(fig, 'Fig2', ['normalization.csv'])

def figure3():
    f=read('pseudocount.csv')
    fig,axs=plt.subplots(2,2,figsize=(6.85,4.95))
    fig.subplots_adjust(left=.11,right=.975,bottom=.18,top=.91,hspace=.73,wspace=.38)
    maximum=f.loc[f.metric.eq('fdp'),'mean'].max()*100
    for ax,route,label,tag in zip(axs.ravel(),ROUTES,'abcd',SHORT):
        letter(ax,label)
        for denom in ['Observed total','Fixed total']:
            g=f[f.reference_source.eq(route[0])&f.target_source.eq(route[1])&f.metric.eq('fdp')&f.denominator.eq(denom)].sort_values('ratio')
            ax.errorbar(g.ratio,g['mean']*100,yerr=np.vstack([(g['mean']-g.low)*100,(g.high-g['mean'])*100]),
                        marker='o' if denom=='Observed total' else 's',ms=3.5,capsize=2,color=COL[denom],label=denom)
        ax.set_xscale('log');ax.set_ylim(0,maximum*1.16);ax.set_ylabel('FDR (%)');ax.set_xlabel('Pseudocount / scale')
        style.route_title(ax,*route)
        ax.axhline(5,color='#777777',ls=':',lw=.7)
    fig.legend(*axs[0,0].get_legend_handles_labels(),loc='lower center',ncol=2,frameon=False,bbox_to_anchor=(.5,.005))
    finish(fig, 'Fig3', ['pseudocount.csv'])

GRID=['f05_a025','f05_a050','f05_a100','f10_a025','f10_a050','f10_a100','f20_a025','f20_a050','f20_a100']

GRIDLAB=[f'{f}%\n{a}' for f in [5,10,20] for a in ['0.25','0.5','1.0']]

def heat(ax,data,metric,outcome,method,paired=False,limits=None,cmap='viridis',ylabel=''):
    matrix=np.empty((4,9))
    for i,(ref,target) in enumerate(ROUTES):
        for j,grid in enumerate(GRID):
            matrix[i,j]=choose(data,reference_source=ref,target_source=target,grid_id=grid,
                               metric=metric,outcome=outcome,method=method)['mean']*100
    if limits is None: limits=(float(matrix.min()),float(matrix.max()))
    image=ax.imshow(matrix,aspect='auto',cmap=cmap,vmin=limits[0],vmax=limits[1],interpolation='nearest')
    ax.set_xticks(range(9),['0.25','0.5','1.0']*3)
    ax.set_yticks(range(4),SHORT);ax.set_ylabel(ylabel,labelpad=8)
    for center,label in zip([1,4,7],['5% markers','10% markers','20% markers']):
        ax.text(center,1.045,label,transform=ax.get_xaxis_transform(),ha='center',va='bottom',fontsize=8)
    ax.set_xticks(np.arange(-.5,9,1),minor=True)
    ax.set_yticks(np.arange(-.5,4,1),minor=True)
    ax.grid(which='minor',color='white',linewidth=.7)
    ax.tick_params(which='minor',length=0)
    for i in range(4):
        for j in range(9):
            pos=(matrix[i,j]-limits[0])/max(1e-12,limits[1]-limits[0])
            rgb=np.array(image.cmap(pos)[:3]);linear=np.where(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)**2.4);lum=float(linear@np.array([.2126,.7152,.0722]))
            shown=0. if abs(matrix[i,j])<.05 else matrix[i,j]
            ax.text(j,i,f'{shown:.1f}',ha='center',va='center',fontsize=8,color='black' if (lum+.05)/.05 >= 1.05/(lum+.05) else 'white')
    return image

def figure4():
    s,p=read('dose_grid.csv'),read('dose_contrasts.csv')
    fig,axs=plt.subplots(3,1,figsize=(6.85,6.15))
    fig.subplots_adjust(left=.13,right=.90,bottom=.10,top=.92,hspace=.76)
    specs=[(p,'fdp','null',True,(-90,90),style.DIVERGING,'Δ FDR (pp)'),
           (s,'fdp','null',False,(0,100),style.SEQUENTIAL,'MRN FDR (%)'),
           (p,'power','up05',True,(-20,20),style.DIVERGING,'Δ power (pp)')]
    for ax,label,(f,metric,outcome,paired,lim,cmap,ylabel) in zip(axs,'abc',specs):
        letter(ax,label);im=heat(ax,f,metric,outcome,'marker_matched',limits=lim,cmap=cmap,ylabel=ylabel)
        cb=fig.colorbar(im,ax=ax,pad=.025,fraction=.028);cb.ax.tick_params(labelsize=8)
        cb.set_ticks([lim[0],(lim[0]+lim[1])/2,lim[1]])
    axs[-1].set_xlabel('Relative-count amplitude',labelpad=7)
    finish(fig, 'Fig4', ['dose_grid.csv', 'dose_contrasts.csv'])

def figure7():
    f=read('power_trajectories.csv');c=read('power_components.csv');d=read('power_component_donors.csv')
    fig,axs=plt.subplots(2,2,figsize=(6.85,5.2))
    fig.subplots_adjust(left=.115,right=.975,bottom=.20,top=.91,hspace=.68,wspace=.45)
    for ax,ref in zip(axs[0],['LIBD','SEAAD']):
        for method in ['tmm','marker_matched']:
            g=f[f.reference_source.eq(ref)&f.method.eq(method)].copy()
            g['effect']=g.condition.map({'up01':.1,'up025':.25,'up05':.5});g=g.sort_values('effect')
            ax.errorbar(g.effect,100*g['mean'],yerr=np.vstack([100*(g['mean']-g.low),100*(g.high-g['mean'])]),
                        marker='o' if method=='tmm' else 's',ms=4,capsize=2,color=COL[method],label=LAB[method])
        ax.set_ylim(0,103);ax.set_xticks([.1,.25,.5]);ax.set_xlabel('Direct effect');ax.set_ylabel('Power (%)')
        style.route_title(ax,ref,ref)
    ax=axs[1,0];part=c[c.reference_source.eq('SEAAD')].set_index('method')
    methods=['tmm','marker_matched','marker_denominator_only'];x=np.arange(3)
    base=np.zeros(3)
    for metric,col,lab in [('leakage_variance','#82B4AB','Composition'),('noise_variance','#A4B9D2','Noise')]:
        val=part.loc[methods,metric].to_numpy()*1000
        ax.bar(x,val,bottom=base,color=col,width=.52,edgecolor='white',linewidth=.6,label=lab);base+=val
    cov=part.loc[methods,'twice_covariance'].to_numpy()*1000
    ax.bar(x,cov,color='#D7A7B3',width=.52,edgecolor='white',linewidth=.6,label='2 × covariance')
    ax.scatter(x,part.loc[methods,'error_variance']*1000,marker='D',s=22,color='#222222',label='Total',zorder=4)
    ax.axhline(0,color='#555555',lw=.7);ax.set_xticks(x,['TMM','MRN','Denom.']);ax.set_ylabel('Variance (×0.001)')
    ax=axs[1,1];part=d[d.reference_source.eq('SEAAD')]
    for donor,g in part.groupby('target_donor'):
        q=g.set_index('method').reindex(['tmm','marker_matched'])
        ax.plot([0,1],q.mean_beta,color='#BCC4CB',marker='o',ms=3,lw=.8)
    for j,method in enumerate(['tmm','marker_matched']):
        ax.scatter(j,float(c[(c.reference_source.eq('SEAAD'))&c.method.eq(method)].mean_beta.iloc[0]),
                   color=COL[method],s=45,zorder=4)
    ax.axhline(.25,ls=':',color='#555555',lw=.8);ax.set_xlim(-.3,1.3);ax.set_xticks([0,1],['TMM','MRN'])
    ax.set_ylabel('Mean coefficient');ax.text(.05,.08,'Truth = 0.25',transform=ax.transAxes,va='bottom',fontsize=8)
    for ax,label in zip(axs.ravel(),'abcd'):letter(ax,label)
    h,l=axs[0,0].get_legend_handles_labels();fig.legend(h,l,ncol=2,loc='lower center',bbox_to_anchor=(.26,.02),frameon=False)
    h,l=axs[1,0].get_legend_handles_labels();fig.legend(h,l,ncol=2,loc='lower center',bbox_to_anchor=(.75,.005),frameon=False)
    finish(fig, 'Fig7', ['power_trajectories.csv', 'power_components.csv', 'power_component_donors.csv'])

def figureS1():
    s,p=read('dose_grid.csv'),read('dose_contrasts.csv')
    fig,axs=plt.subplots(3,1,figsize=(6.85,6.15));fig.subplots_adjust(left=.13,right=.90,bottom=.10,top=.92,hspace=.76)
    specs=[(s,'fdp',(0,100),style.SEQUENTIAL,'MRN FDR (%)'),(s,'coverage',(0,100),style.SEQUENTIAL,'MRN coverage (%)'),
           (p,'coverage',(-45,45),style.DIVERGING,'Δ coverage (pp)')]
    for ax,label,(f,metric,lim,cm,yl) in zip(axs,'abc',specs):
        letter(ax,label);im=heat(ax,f,metric,'up05','marker_matched',limits=lim,cmap=cm,ylabel=yl)
        cb=fig.colorbar(im,ax=ax,pad=.025,fraction=.028)
        cb.set_ticks([lim[0],(lim[0]+lim[1])/2,lim[1]])
    axs[-1].set_xlabel('Relative-count amplitude')
    finish(fig, 'FigS1', ['dose_grid.csv', 'dose_contrasts.csv'])

def figureS5():
    f=read('clinical_markers.csv');fig,axs=plt.subplots(1,2,figsize=(6.85,3.4))
    fig.subplots_adjust(left=.11,right=.98,bottom=.26,top=.88,wspace=.39)
    classes=['Astro','Excit','Inhib','Micro','OPC','Oligo']
    rng=np.random.default_rng(27)
    for ax,ref,label in zip(axs,['LIBD','SEAAD'],'ab'):
        for j,ct in enumerate(classes):
            g=f[f.reference.eq(ref)&f.cell_type.eq(ct)&f.available]
            ax.scatter(j+rng.uniform(-.17,.17,len(g)),g.log2FC,s=14,color=style.CLASS_COLORS[j],alpha=.75,
                       edgecolors='white',linewidths=.25)
            if len(g):ax.plot([j-.19,j+.19],[g.log2FC.median()]*2,color='#222222',lw=1)
        ax.axhline(0,color='#777777',lw=.7)
        for y in [-np.log2(1.5),np.log2(1.5)]:ax.axhline(y,color='#777777',ls=':',lw=.7)
        ax.set_xticks(range(6),classes,rotation=35,ha='right');ax.set_ylabel('AD − control log2 fold change')
        ax.text(.5,1.035,ref.replace('SEAAD','SEA-AD')+' panel',transform=ax.transAxes,va='bottom',ha='center',fontsize=9);letter(ax,label)
    finish(fig, 'FigS5', ['clinical_markers.csv'])

def figureS6():
    f=read('actual_span.csv');fig,axs=plt.subplots(2,2,figsize=(6.85,4.95))
    fig.subplots_adjust(left=.11,right=.98,bottom=.19,top=.91,hspace=.60,wspace=.43)
    for ax,label,metric,yl in zip(axs.ravel(),'abcd',['availability','coverage','delivered_coverage','fdp'],
                                ['Test availability (%)','Conditional coverage (%)','Delivered coverage (%)','FDR (%)']):
        for j,rule in enumerate(['full_six','adaptive']):
            means=[];lo=[];hi=[]
            for ref,target in ROUTES:
                r=choose(f,reference_source=ref,target_source=target,rule=rule,metric=metric)
                means.append(r['mean']*100);lo.append((r['mean']-r.low)*100);hi.append((r.high-r['mean'])*100)
            ax.errorbar(np.arange(4)+(j-.5)*.16,means,yerr=[lo,hi],fmt='o' if j==0 else 's',capsize=2,
                        color=style.GRAY if j==0 else style.VIOLET,label='Six-contrast rule' if j==0 else 'Actual-span rule')
        routes_axes(ax);ax.set_ylabel(yl);letter(ax,label)
        if metric in ['coverage','delivered_coverage']:ax.axhline(95,color='#999999',ls=':',lw=.7)
        if metric=='fdp':ax.axhline(5,color='#999999',ls=':',lw=.7)
    fig.legend(*axs[0,0].get_legend_handles_labels(),ncol=2,loc='lower center',frameon=False,bbox_to_anchor=(.5,.005))
    finish(fig, 'FigS6', ['actual_span.csv'])
