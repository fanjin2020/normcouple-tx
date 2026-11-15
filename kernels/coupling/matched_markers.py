"""Reference-only matching and paired integer budgets with explicit support rules."""
import hashlib
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from common import require

MATCH_SEED = 2026103101
FEATURES = ['log_abundance', 'detection', 'specificity']


def match_indices(features, owners, genes, selected, seed):
    features = np.asarray(features, float)
    owners, genes, selected = np.asarray(owners), np.asarray(genes), np.asarray(selected, int)
    require(features.shape == (len(owners), 3) and np.isfinite(features).all(), 'Invalid reference features')
    rng = np.random.default_rng(seed)
    chosen, pairs = [], []
    for k in range(1, 8):
        ix = np.flatnonzero(owners == k)
        source = np.array(sorted(selected[owners[selected] == k], key=lambda i: genes[i]))
        candidates = np.array(sorted(np.setdiff1d(ix, selected), key=lambda i: genes[i]))
        require(len(source) > 0 and len(candidates) >= len(source), 'Insufficient distinct class-matched candidates')
        scale = features[ix].std(0)
        scale = np.where(scale > 1e-12, scale, 1.)
        gap = (features[source, None, :] - features[candidates][None, :, :]) / scale
        distance = np.sum(gap**2, axis=2)
        # Small, prespecified random jitter supplies replicate optimal assignments.
        rows, cols = linear_sum_assignment(distance + rng.uniform(0, .1, size=distance.shape))
        chosen.extend(candidates[cols].tolist())
        for a, b in zip(rows, cols):
            record = dict(cell_class=k, original_gene=str(genes[source[a]]),
                          matched_gene=str(genes[candidates[b]]), distance=float(np.sqrt(distance[a, b])),
                          weak_match=bool(np.max(abs(gap[a, b])) > 1.))
            record.update({f+'__standardized_gap':float(-gap[a,b,j]) for j,f in enumerate(FEATURES)})
            pairs.append(record)
    result = np.sort(chosen)
    require(len(result) == len(selected) and len(np.unique(result)) == len(result), 'Invalid matching assignment')
    require(not set(result) & set(selected), 'Matched set overlaps predefined set')
    return result, pairs


def reference_features(job):
    import revision_core as C
    marker = np.asarray(job['panels']['balanced'])
    source = job['source']
    owners = C.E.official_markers(job).cell_type.to_numpy(int)
    means = np.asarray(job['base']['mean_counts'])[marker]
    own = means[np.arange(len(marker)), owners-1]
    specificity = own / np.maximum(means.sum(1), 1e-12)
    detection = np.zeros(len(marker))
    for k in range(1, 8):
        ix = np.flatnonzero(owners == k)
        fractions = []
        for donor in job['refs']:
            cells = source['masks'][donor][k-1]
            require(len(cells) > 0, 'Reference donor missing a required class')
            fractions.append(np.asarray((source['x'][marker[ix]][:, cells] > 0).mean(axis=1)).ravel())
        detection[ix] = np.mean(fractions, axis=0)
    return np.column_stack([np.log1p(own), detection, specificity]), owners


def freeze_matched(jobs, number, output):
    import revision_core as C
    features_all, pairs_all, balances, manifests = [], [], [], []
    for index, job in enumerate(jobs):
        features, owners = reference_features(job)
        marker = job['panels']['balanced']
        genes = job['source']['genes'][marker]
        original = np.asarray(job['stress_subsets']['legacy'], int)
        context = dict(job_index=index, reference_source=job['source']['label'],
                       target_source=job.get('target_source',job['source'])['label'])
        for i,gene in enumerate(genes):
            features_all.append(dict(context,gene_id=str(gene),cell_class=int(owners[i]),
                                     **dict(zip(FEATURES,features[i].tolist()))))
        for j in range(number):
            key = C.M.digest([MATCH_SEED,context,genes.tolist(),j])
            seed = int(key[:16],16)
            selected, pairs = match_indices(features,owners,genes,original,seed)
            name = 'matched_%03d' % (j+1)
            job['stress_subsets'][name] = selected
            pairs_all.extend(dict(context,subset_id=name,**row) for row in pairs)
            manifests.append(dict(context,subset_id=name,selection_seed=str(seed),
                selected_gene_ids=genes[selected].tolist(),selected_class_ids=owners[selected].tolist(),
                selection_rule='within-class reference feature matching; predefined genes excluded; no outcome access'))
        for name,selected in job['stress_subsets'].items():
            if name == 'zero':
                continue
            for k in range(1,8):
                ix = np.flatnonzero(owners == k)
                a = original[owners[original] == k]
                b = selected[owners[selected] == k]
                sd = features[ix].std(0)
                sd = np.where(sd>1e-12,sd,1.)
                gap = (features[b].mean(0)-features[a].mean(0))/sd
                balances.append(dict(context,subset_id=name,cell_class=k,selected_count=len(b),
                    overlap_with_predefined=len(set(b)&set(a)),
                    **{f+'__standardized_mean_difference':float(gap[i]) for i,f in enumerate(FEATURES)}))
    output.mkdir(parents=True,exist_ok=True)
    for name,rows in [('reference_marker_features',features_all),('matching_pairs',pairs_all),
                       ('marker_matching_balance',balances)]:
        pd.DataFrame(rows).to_csv(output/(name+'.csv'),index=False)
    return manifests


def variants(subsets):
    return ['zero'] + [name+'__'+rule for name in subsets if name!='zero'
                      for rule in ('allow_zero','positive_only')]


def split_variant(variant):
    return ('zero','none') if variant=='zero' else tuple(variant.rsplit('__',1))


def allocate(weights, selected, owners, budget, observed, rule, seed):
    weights, observed = np.asarray(weights,float), np.asarray(observed)
    selected, owners, budget = np.asarray(selected,int),np.asarray(owners,int),np.asarray(budget)
    require(rule in ('allow_zero','positive_only'),'Unknown support rule')
    require(observed.ndim==2 and observed.shape[1]==len(weights),'Invalid observed shape')
    require(np.isfinite(observed).all() and (observed>=0).all(),'Invalid observed counts')
    require(budget.shape==(len(observed),7) and np.isfinite(budget).all() and
            (budget>=0).all() and np.array_equal(budget,np.rint(budget)) and
            (budget<np.iinfo(np.int64).max).all(),'Invalid integer budget')
    require(np.isfinite(weights).all() and (weights>=0).all() and owners.shape==weights.shape and
            np.isin(owners,np.arange(1,8)).all(),'Invalid fixed reference weights/classes')
    require(len(selected)==len(set(selected)) and (selected>=0).all() and
            (selected<len(weights)).all(),'Invalid selected genes')
    unsupported = []
    value = np.zeros_like(observed,dtype=np.int64)
    rng = np.random.default_rng(seed)
    for k in range(1,8):
        ix = np.sort(selected[owners[selected]==k])
        require(len(ix)>0,'Selected set lacks a cell class')
        for sample in np.flatnonzero(budget[:,k-1]>0):
            eligible = ix[weights[ix]>0]
            if rule=='positive_only':
                eligible = eligible[observed[sample,eligible]>0]
            if not len(eligible):
                unsupported.append((int(sample),k,int(budget[sample,k-1])))
                continue
            w = weights[eligible]
            value[sample,eligible] = rng.multinomial(int(budget[sample,k-1]),w/w.sum())
    if unsupported:
        return None,unsupported
    actual = np.column_stack([value[:,owners==k].sum(1) for k in range(1,8)])
    require(np.array_equal(actual,budget),'Realized class budget changed')
    require(rule!='positive_only' or not value[observed==0].any(),'Zero-support restriction violated')
    return value,[]


def stress(job,base,variant):
    import revision_core as C
    import selection as S
    from reference_allocation import frozen_weights
    subset,rule = split_variant(variant)
    marker = job['panels']['balanced']
    observed = base['draw']['counts'][:,marker]
    owners = C.E.official_markers(job).cell_type.to_numpy(int)
    selected = np.asarray(job['stress_subsets'][subset],int)
    legacy,_ = S.stress(job,base,'legacy')
    budget = np.column_stack([legacy[:,owners==k].sum(1) for k in range(1,8)])
    if subset=='zero':
        budget[:] = 0
        value,unsupported = np.zeros_like(observed,dtype=np.int64),[]
    else:
        seed = int(C.M.digest([MATCH_SEED,'allocation',base['key'],subset,rule])[:16],16)
        value,unsupported = allocate(frozen_weights(job,owners,'reference'),selected,owners,budget,
                                     observed,rule,seed)
    case = base['group']==1
    selected_mass = float(observed[case][:,selected].sum())
    selected_case = observed[case][:,selected]
    selected_control = observed[~case][:,selected]
    record = dict(base_subset=subset,subset_kind=subset.split('_')[0],support_rule=rule,
        intervention_feasible=value is not None,selected_markers=len(selected),
        budget_sha256=hashlib.sha256(budget.astype('<i8').tobytes()).hexdigest(),
        requested_budget_total=int(budget.sum()),unsupported_sample_class_count=len(unsupported),
        unsupported_sample_count=len({v[0] for v in unsupported}),
        unsupported_budget_total=sum(v[2] for v in unsupported),
        selected_baseline_case_mass=selected_mass,
        selected_baseline_control_mass=float(selected_control.sum()),
        target_mixture_mean_selected_count=float(selected_case.mean()) if selected_case.size else np.nan,
        target_mixture_selected_zero_fraction=float((selected_case==0).mean()) if selected_case.size else np.nan,
        target_control_selected_zero_fraction=float((selected_control==0).mean()) if selected_control.size else np.nan,
        effective_relative_amplitude=float(budget.sum()/selected_mass) if selected_mass else np.nan,
        total_added_counts=int(value.sum()) if value is not None else np.nan,
        added_counts_on_observed_zero=int(value[observed==0].sum()) if value is not None else np.nan,
        infeasibility_reason='NO_POSITIVE_SELECTED_SUPPORT' if unsupported else '')
    if value is not None:
        require(not value[~case].any(),'Controls received perturbation')
    return value,record
