"""Coefficient-change identities and local sensitivities, not causal bias estimates."""
import numpy as np


def geometry(x, p, covariates=None):
    x, p = np.asarray(x, float), np.asarray(p, float)
    if p.shape != (len(x), 7) or not np.isfinite(p).all():
        raise ValueError("Nonfinite or invalid composition matrix")
    pieces = [np.ones((len(x), 1))]
    if covariates is not None:
        pieces.append(np.asarray(covariates, float))
    # On the simplex the last column is redundant with the intercept.
    pieces.append(p[:, :6])
    w = np.column_stack(pieces)
    u, s, _ = np.linalg.svd(w, full_matrices=False)
    rank = int((s > max(1e-12, s[0] * 1e-10)).sum())
    m = np.eye(len(x)) - u[:, :rank] @ u[:, :rank].T
    z = m @ x
    zz = float(z @ z)
    if zz <= max(1e-12, np.sum((x-x.mean())**2)*1e-10):
        raise ValueError("Exposure is not estimable")
    return w, m, z, zz, rank


def coefficient(x, y, p, covariates=None):
    _, _, z, zz, _ = geometry(x, p, covariates)
    return z @ np.asarray(y) / zz


def geometry_change(x, p0, p1, covariates=None):
    _, m0, z0, ss0, r0 = geometry(x, p0, covariates)
    _, m1, z1, ss1, r1 = geometry(x, p1, covariates)
    return dict(nuisance_rank_before=r0, nuisance_rank_after=r1,
        nuisance_space_distance=float(np.linalg.norm(m1-m0, 'fro')),
        residualized_exposure_change=float(np.linalg.norm(z1-z0)),
        coefficient_weight_change=float(np.linalg.norm(z1/ss1-z0/ss0)),
        exposure_information_before=ss0, exposure_information_after=ss1)


def exact_change(x, y, p0, p1, covariates=None):
    _, m0, z0, zz0, _ = geometry(x, p0, covariates)
    _, m1, z1, zz1, _ = geometry(x, p1, covariates)
    b0, b1 = z0 @ y / zz0, z1 @ y / zz1
    predicted = x @ (m1-m0) @ (y-x[:, None]*b0) / zz1
    return b1-b0, predicted


def coefficient_derivative(x, y, p, dp, covariates=None):
    w, m, z, zz, rank = geometry(x, p, covariates)
    if rank != w.shape[1]:
        raise ValueError("Local derivative requires full, stable nuisance rank")
    dw = np.zeros_like(w)
    dw[:, -6:] = dp[:, :6]
    plus = np.linalg.pinv(w, rcond=1e-10)
    dm = -(m @ dw @ plus + plus.T @ dw.T @ m)
    beta = z @ y / zz
    derivative = x @ dm @ (y-x[:, None]*beta) / zz
    return derivative


def local_report(x, y, baseline, probes, h, se, covariates=None):
    """No true composition, effect, noise or simulated truth enters these scores."""
    p0 = baseline
    w, _, z, zz, rank = geometry(x, p0, covariates)
    ranks = [geometry(x, p, covariates)[-1] for t,p in probes.items() if abs(t)<=h]
    stable = rank == w.shape[1] and all(r == rank for r in ranks)
    dp = (probes[h]-probes[-h])/(2*h)
    dph = (probes[h/2]-probes[-h/2])/h
    scale = np.maximum(np.asarray(se, float), 1e-12)
    report = dict(stable_nuisance_rank=stable, nuisance_rank=rank,
                  exposure_information=zz,
                  composition_direction_norm=float(np.linalg.norm(dph)/np.sqrt(len(x))),
                  exposure_aligned_composition=float(np.linalg.norm(z @ dph)/zz))
    if not stable:
        report.update(score=np.nan, derivative_agreement=np.nan,
                      status="UNSTABLE_OR_DEFICIENT_NUISANCE_RANK")
        return report, []
    db = coefficient_derivative(x, y, p0, dp, covariates)
    dbh = coefficient_derivative(x, y, p0, dph, covariates)
    observed_derivative = (coefficient(x,y,probes[h/2],covariates)-
                           coefficient(x,y,probes[-h/2],covariates))/h
    agreement = float(np.linalg.norm(dbh-db)/max(np.linalg.norm(dbh), 1e-10))
    report.update(score=float(np.quantile(abs(dbh)/scale, .90)),
                  derivative_agreement=agreement,
                  derivative_chain_absolute_error=float(np.max(abs(dbh-observed_derivative))),
                  local_smoothness_warning=agreement > .10,
                  status="LOCAL_SENSITIVITY")
    rows = []
    for t, props in probes.items():
        if t <= h:
            continue
        try:
            probe_rank=geometry(x,props,covariates)[-1]
            actual, exact = exact_change(x, y, p0, props, covariates)
        except (ValueError,np.linalg.LinAlgError) as exc:
            rows.append(dict(log_count_step=t,status='UNAVAILABLE:'+str(exc)))
            continue
        linear = t*dbh
        rows.append(dict(log_count_step=t,status='ESTIMABLE',probe_nuisance_rank=probe_rank,
            **geometry_change(x,p0,props,covariates),
            probe_rank_changed=probe_rank!=rank,
            exact_identity_max_error=float(np.max(abs(actual-exact))),
            observed_change_rms=float(np.sqrt(np.mean(actual**2))),
            predicted_change_rms=float(np.sqrt(np.mean(linear**2))),
            linear_error_rms=float(np.sqrt(np.mean((linear-actual)**2))),
            relative_linear_error=float(np.linalg.norm(linear-actual)/max(np.linalg.norm(actual),1e-10)),
            standardized_change_p90=float(np.quantile(abs(actual)/scale,.90))))
    return report, rows
