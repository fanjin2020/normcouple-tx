"""Figures from exported numbers; no analysis is silently rerun."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COLORS={"observed":"#D55E00","fixed":"#777777","tmm":"#0072B2",
        "marker_matched":"#009E73","marker_unmatched_scale":"#CC79A7","marker_denominator_only":"#AA7733"}


def save(fig, out, name):
    fig.savefig(out/(name+".pdf"),bbox_inches="tight")
    fig.savefig(out/(name+".png"),dpi=200,bbox_inches="tight")
    plt.close(fig)


def figures(root):
    root=Path(root); out=root/"figures";out.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42})
    f=pd.read_csv(root/"summary.csv",keep_default_na=False,na_values=[""])
    groups=f[["reference_source","target_source"]].drop_duplicates().itertuples(index=False,name=None)
    for source,target in groups:
        x=f[f.reference_source.eq(source)&f.target_source.eq(target)&f.condition.eq("up05")&f.rule.eq("full_six")]
        fig,axes=plt.subplots(1,3,figsize=(10,3),layout="constrained")
        for ax,metric,label in zip(axes,["fdp","null_coverage","availability"],["Mean FDP (%)","Unspiked coverage (%)","Available tests (%)"]):
            for i,method in enumerate(COLORS):
                row=x[x.method.eq(method)&x.metric.eq(metric)]
                if row.empty:continue
                r=row.iloc[0]
                if np.isfinite(r["mean"]):
                    ax.errorbar(i,100*r["mean"],yerr=196*r.mcse if np.isfinite(r.mcse) else None,fmt="o",capsize=3,color=COLORS[method])
            ax.set_xticks(range(len(COLORS)),["Observed","Fixed","TMM","Marker\nmatched","Marker\nscale control","Denominator\nonly"],rotation=35,ha="right")
            ax.set_ylabel(label);ax.grid(axis="y",alpha=.2)
            if metric in ("fdp","null_coverage"): ax.axhline(5 if metric=="fdp" else 95,color=".5",ls="--")
        fig.suptitle(source+" reference to "+target+" targets")
        save(fig,out,"normalization_"+source+"_"+target)
    p=pd.read_csv(root/"paired_summary.csv",keep_default_na=False,na_values=[""])
    x=p[p.condition.eq("up05")&p.rule.eq("adaptive")&p.reference_rule.eq("full_six")&p.method.eq("mead_reference_converted")]
    fig,axes=plt.subplots(1,3,figsize=(9,3),layout="constrained")
    for ax,metric in zip(axes,["availability","delivered_coverage","fdp"]):
        d=x[x.metric.eq(metric)]
        if len(d):
            ax.errorbar(d["mean"]*100,np.arange(len(d)),xerr=d.mcse*196,fmt="o",color="#0072B2",capsize=3)
            ax.set_yticks(range(len(d)),d.reference_source+" to "+d.target_source)
        ax.axvline(0,color=".5",ls="--");ax.set_xlabel("Adaptive minus full-six\n"+metric+" (percentage points)")
    save(fig,out,"MEAD_column_space_sensitivity")
    d=pd.read_csv(root/"donor_summary.csv",keep_default_na=False,na_values=[""])
    d=d[d.condition.eq("up05")&d.metric.eq("fdp")&d.method.isin(["observed","tmm","marker_matched"])]
    fig,ax=plt.subplots(figsize=(7,4),layout="constrained")
    methods=["observed","tmm","marker_matched"]
    for _,f in d.groupby(["reference_source","target_source","target_donor"]):
        vals=f.set_index("method")["mean"].reindex(methods).to_numpy(float)*100
        ax.plot(range(3),vals,color=".65",alpha=.35,lw=.7,marker=".")
    ax.set_xticks(range(3),methods);ax.set_ylabel("Donor mean FDP (%)");ax.axhline(5,ls="--",color=".5")
    save(fig,out,"paired_donor_heterogeneity")


def bulk_figures(root):
    root=Path(root);out=root/"figures";out.mkdir(exist_ok=True)
    a=pd.read_csv(root/"disease_associations.csv")
    for source in ("LIBD","SEAAD"):
        base=a[a.method.eq(source+"_observed")].set_index("gene_id")
        fig,axes=plt.subplots(1,2,figsize=(8,3.5),layout="constrained")
        for ax,arm in zip(axes,("tmm","marker_matched")):
            f=a[a.method.eq(source+"_"+arm)].set_index("gene_id").reindex(base.index)
            good=base.available & f.available
            ax.scatter(base.beta[good],f.beta[good],s=4,alpha=.2,color=COLORS[arm],rasterized=True)
            vals=np.r_[base.beta[good],f.beta[good]]
            if len(vals): ax.plot([vals.min(),vals.max()],[vals.min(),vals.max()],color=".4",ls="--")
            ax.set_xlabel("Observed-total input: AD coefficient")
            ax.set_ylabel(arm+" input: AD coefficient")
        fig.suptitle(source+" reference; same bulk cohort and expression outcome")
        save(fig,out,"bulk_coefficient_stability_"+source)


if __name__=="__main__":
    import argparse
    p=argparse.ArgumentParser();p.add_argument("results",type=Path)
    figures(p.parse_args().results)
