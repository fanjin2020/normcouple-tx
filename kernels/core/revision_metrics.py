"""All-planned denominators, FWL audit, and paired method contrasts."""
import numpy as np
import pandas as pd
from scipy import stats
from revision_core import M, require

BASELINE = "dt_observed_a1_c1000000"


def evaluate(data, props, reason=""):
    fit = M.model_fit(data["y"], data["group"], props, data["nuisance"], data["noise"], data["truth"])
    valid = fit["good"] & fit["basic"]
    truth = data["truth"]
    null = truth == 0
    p = np.ones(len(truth))
    p[valid] = 2 * stats.t.sf(abs(fit["beta"][valid] / fit["se"][valid]), fit["df"])
    q = M.bh(p)
    hit = q <= .05
    critical = stats.t.ppf(.975, fit["df"])
    covered = valid & (abs(fit["beta"] - truth) <= critical * fit["se"])
    row = dict(availability=float(valid.mean()), complete_fit=int(valid.all()),
               fdp=float((hit & null).sum()/max(1, hit.sum())), any_false=float((hit & null).any()),
               power=float((hit & ~null).sum()/(~null).sum()) if (~null).any() else np.nan,
               coverage=float(covered[valid].mean()) if valid.any() else np.nan,
               delivered_coverage=float(covered.mean()), null_coverage=float(covered[valid & null].mean()) if (valid & null).any() else np.nan,
               null_delivered_coverage=float(covered[null].mean()) if null.any() else np.nan,
               null_availability=float(valid[null].mean()) if null.any() else np.nan,
               gene_error_mean=float((fit["beta"]-truth)[valid].mean()) if valid.any() else np.nan,
               abs_leakage=float(abs(fit["leakage"][valid]).mean()) if valid.any() else np.nan,
               mean_se=float(fit["se"][valid].mean()) if valid.any() else np.nan,
               identity_error=fit["identity_error"], df=fit["df"], leverage=fit["leverage"],
               design_condition=fit["condition"], reason="; ".join(x for x in (reason,fit["reason"]) if x),
               allocation=data["allocation"], n_tests=len(truth), n_alternatives=int((~null).sum()))
    row.update(M.prop_diagnostics(props,data["p"],data["group"]) if props is not None else dict(
        proportion_rmse=np.nan,group_error_contrast=np.nan,boundary_fraction=np.nan))
    contrast = np.full(7,np.nan)
    if props is not None and np.isfinite(props).all():
        error = props-data["p"]
        contrast = error[data["group"]==1].mean(0)-error[data["group"]==0].mean(0)
    row.update({"error_contrast_"+ct:float(value) for ct,value in zip(M.CT,contrast)})
    arrays = dict(truth=truth, beta=np.where(valid,fit["beta"],np.nan), se=np.where(valid,fit["se"],np.nan),
                  leakage=np.where(valid,fit["leakage"],np.nan), noise_projection=np.where(valid,fit["noise_projection"],np.nan),
                  available=valid.astype(float),covered=covered.astype(float),p=p,q=q)
    return row, arrays


ALL_PLANNED = ["availability", "fdp", "any_false", "power", "delivered_coverage", "null_delivered_coverage", "null_availability"]
CONDITIONAL = ["coverage", "null_coverage", "gene_error_mean", "abs_leakage", "proportion_rmse", "mean_se"]


def pair_one(records, arrays):
    """One condition/replicate: method A minus original dtangle, same outcomes."""
    rows = []
    specs = [(method,BASELINE) for method in records if method!=BASELINE]
    specs += [(method,method.replace("dt_observed_","dt_fixed_")) for method in records
              if method.startswith("dt_observed_") and method.replace("dt_observed_","dt_fixed_") in records]
    specs += [(method,method.replace("mead_","mead_lambda0_",1)) for method in ("mead_native","mead_reference_converted")
              if method in records and method.replace("mead_","mead_lambda0_",1) in records]
    for method, reference in specs:
        record,baseline = records[method],records[reference]
        common_complete = baseline["complete_fit"] and record["complete_fit"]
        for metric in ALL_PLANNED + CONDITIONAL:
            allow = metric in ALL_PLANNED or common_complete
            rows.append(dict(method_a=method,method_b=reference,metric=metric,
                delta=record[metric]-baseline[metric] if allow else np.nan,
                common_complete=int(common_complete),denominator="all_planned" if metric in ALL_PLANNED else "common_complete_fits"))
        common = (arrays[method]["available"] > 0) & (arrays[reference]["available"] > 0)
        for metric in ("beta","leakage","se"):
            delta = arrays[method][metric]-arrays[reference][metric]
            rows.append(dict(method_a=method,method_b=reference,metric=metric+"_common_gene_mean",
                delta=float(delta[common].mean()) if common.any() else np.nan,
                common_complete=int(common_complete),denominator="common_available_genes",common_genes=int(common.sum())))
    return rows


def per_gene_summary(frame, planned):
    out = []
    for (condition, method, gene), f in frame.groupby(["condition","method","gene_id"],sort=False):
        require(len(f)==planned and f.rep.nunique()==planned, "Missing/duplicate gene repetitions")
        valid = f.available.eq(1)
        mean = f.loc[valid,"beta"].mean()
        sd = f.loc[valid,"beta"].std(ddof=1)
        out.append(dict(condition=condition,method=method,gene_id=gene,truth=f.truth.iloc[0],
            available=int(valid.sum()),planned=planned,bias=mean-f.truth.iloc[0],empirical_sd=sd,
            mean_se=f.loc[valid,"se"].mean(),rms_se=np.sqrt((f.loc[valid,"se"]**2).mean()),
            mean_leakage=f.loc[valid,"leakage"].mean(),mean_noise_projection=f.loc[valid,"noise_projection"].mean(),
            coverage=f.loc[valid,"covered"].mean(),delivered_coverage=f.covered.mean(),availability=valid.mean()))
    return pd.DataFrame(out)


def mc_summary(values):
    x = np.asarray(values,dtype=float)
    x = x[np.isfinite(x)]
    mean = float(x.mean()) if len(x) else np.nan
    se = float(x.std(ddof=1)/np.sqrt(len(x))) if len(x)>1 else np.nan
    return dict(mean=mean,mcse=se,mc_low=mean-1.96*se,mc_high=mean+1.96*se,contributing_replicates=len(x))
