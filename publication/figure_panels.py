"""Shared figure panels; the renderer supplies data and export callbacks."""

from decimal import Decimal, ROUND_HALF_UP
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import figure_style as style

style.apply()
COL=style.METHOD_COLORS

LAB={'tmm':'TMM','marker_matched':'MRN','marker_denominator_only':'Denominator only','oracle':'Oracle'}

ROUTES=[('LIBD','LIBD'),('SEAAD','SEAAD'),('LIBD','GSE174367'),('SEAAD','GSE174367')]

SHORT=['L-in','S-in','L-ext','S-ext']

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
    fig=plt.figure(figsize=(6.85,4.88))
    gs=fig.add_gridspec(2,2,height_ratios=[1.36,1],hspace=.62,wspace=.48,left=.108,right=.973,top=.928,bottom=.135)
    a=fig.add_subplot(gs[0,:]);a.set_axis_off(); a.set(xlim=(0,1),ylim=(0,1))
    letter(a,'a')
    centers=np.linspace(.085,.912,5)
    paths=[
        (.70,.94,'Whole-library TMM',style.BLUE,style.ICE,
         ['Test\ncounts','Effective\ndenominator','Normalized\nmarkers','Estimated\nfractions','Disease\ncoefficient']),
        (.25,.48,'Marker-restricted normalization',style.TEAL,style.MINT,
         ['Test\ncounts','Marker-only\ndenominator','Same marker\ninput','Same\nfractions','Same\ncoefficient'])
    ]
    for y,h,title,c,bg,labels in paths:
        a.text(0,h,title,ha='left',va='center',color=c,weight='bold',fontsize=9)
        for j,(x,label) in enumerate(zip(centers,labels)):
            blocked=y<.5 and j==0
            dark=(y>.5 and j in (0,4)) or (y<.5 and j==1)
            fc='#F0F2F4' if blocked else c if dark else bg
            a.add_patch(plt.Rectangle((x-.081,y-.145),.162,.29,fc=fc,ec='none',zorder=1))
            a.text(x,y,label,ha='center',va='center',fontsize=8.4,weight='bold' if dark else 'normal',
                   color='white' if dark else style.MUTED if blocked else style.INK)
            if j<4:
                mid=(x+centers[j+1])/2
                if blocked:
                    a.plot(mid,y,marker='x',mew=1.8,ms=6,color=style.ROSE)
                else:
                    a.annotate('',(centers[j+1]-.089,y),(x+.089,y),
                               arrowprops=dict(arrowstyle='-|>',color=c,lw=1.2,shrinkA=0,shrinkB=0))
    a.text(.0,-.015,'Held fixed: marker counts, expression responses, samples and tested genes.',
           va='bottom',fontsize=8.3,color=style.MUTED)
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
    finish(fig, 'Fig1')

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
    finish(fig, 'Fig2')

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
    finish(fig, 'Fig3')

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
    for center,label,tint in zip([1,4,7],['5% markers','10% markers','20% markers'],
                                ['#EDF3F9','#DEEAF5','#CDDFF0']):
        ax.text(center,1.049,label,transform=ax.get_xaxis_transform(),ha='center',va='bottom',fontsize=8,
                weight='bold',bbox=dict(facecolor=tint,edgecolor='none',pad=2.5))
    ax.set_xticks(np.arange(-.5,9,1),minor=True)
    ax.set_yticks(np.arange(-.5,4,1),minor=True)
    ax.grid(which='minor',color='white',linewidth=.7)
    ax.tick_params(which='minor',length=0)
    for boundary in (2.5,5.5):
        ax.axvline(boundary,c='white',lw=3.8)
    if metric=='fdp' and limits==(0,100):
        ii,jj=np.where(matrix>5+1e-10)
        ax.scatter(jj-.33,ii-.24,s=12,marker='o',facecolors='none',
                   edgecolors=style.INK,linewidths=.65,zorder=4)
    for i in range(4):
        for j in (2,5,8):
            pos=(matrix[i,j]-limits[0])/max(1e-12,limits[1]-limits[0])
            rgb=np.array(image.cmap(pos)[:3]);linear=np.where(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)**2.4);lum=float(linear@np.array([.2126,.7152,.0722]))
            # Discard floating-point noise at decimal ties before rounding the displayed label.
            shown = Decimal(f'{matrix[i,j]:.12f}').quantize(Decimal('0.1'), rounding=ROUND_HALF_UP)
            ax.text(j,i,format(shown if shown else abs(shown), '.1f'),ha='center',va='center',fontsize=8,color='black' if (lum+.05)/.05 >= 1.05/(lum+.05) else 'white')
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
        if label in ('a','b'):
            ax.add_patch(plt.Rectangle((7.53,2.53),.94,.94,fill=False,ec='white',lw=3))
            ax.add_patch(plt.Rectangle((7.53,2.53),.94,.94,fill=False,ec=style.INK,lw=1.1))
        cb=fig.colorbar(im,ax=ax,pad=.025,fraction=.028);cb.ax.tick_params(labelsize=8)
        cb.set_ticks([0,25,50,75,100] if label=='b' else [lim[0],0,lim[1]])
    axs[-1].set_xlabel('Relative-count amplitude',labelpad=7)
    fig.text(.49,.014,'Circle: absolute FDR above 5%    |    Differences: MRN minus TMM',ha='center',fontsize=8)
    finish(fig, 'Fig4')

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
    old,new=part.loc['tmm'],part.loc['marker_matched']
    variance_change=new.error_variance-old.error_variance
    components=['leakage_variance','noise_variance','twice_covariance']
    changes=np.array([new[z]-old[z] for z in components])*1000
    if not np.isclose(changes.sum()/1000,variance_change,rtol=1e-7,atol=1e-12):
        raise ValueError('Variance increment does not equal its signed components')
    start=1000*old.error_variance
    levels=np.r_[start,start+np.cumsum(changes)]
    ax.bar(0,start,color=style.BLUE,width=.55)
    for j,(lo,hi,delta) in enumerate(zip(levels[:-1],levels[1:],changes),1):
        ax.bar(j,abs(delta),bottom=min(lo,hi),color=style.TEAL if j==1 else style.ROSE,width=.55)
    ax.bar(4,1000*new.error_variance,color=style.TEAL,width=.55)
    for j,v in enumerate(levels):
        ax.plot([j+.275,j+.725],[v,v],color=style.GRAY,lw=.7,ls='--',zorder=3)
    for j,val in enumerate([start,*changes,levels[-1]]):
        height=max(levels[max(0,j-1)],levels[min(j,3)])
        text=f'{val:+.2f}' if 0<j<4 else f'{val:.2f}'
        ax.annotate(text,(j,height),xytext=(0,4),textcoords='offset points',ha='center',fontsize=7.7)
    ax.set_xticks(range(5),['TMM','Comp.','Noise','2 × Cov.','MRN'])
    ax.set_ylabel('Error variance (×0.001)')
    ax.set_ylim(0,levels.max()*1.34)
    if variance_change>0:
        ax.text(.025,.97,f'{100*changes[0]/(1000*variance_change):.2f}% of increase: composition',
                transform=ax.transAxes,va='top',fontsize=8,weight='bold')
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
    fig.text(.74,.038,'Signed TMM-to-MRN increments',ha='center',color=style.MUTED,fontsize=8)
    finish(fig, 'Fig7')

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
    fig.text(.49,.015,'Open circle: absolute FDR estimate above 5% (panel a).',ha='center',fontsize=8)
    finish(fig, 'FigS1')

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
    finish(fig, 'FigS5')

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
    finish(fig, 'FigS6')
