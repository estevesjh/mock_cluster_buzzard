#!/usr/bin/env python
"""Cluster miscentering model.

Forward-models the effect of optical-finder (redMaPPer) center misidentification
on the stacked weak-lensing profile, following the standard literature
construction (Yang+2006; Johnston+2007 Eq. 32; Ford+2014; McClintock+2019
Sec. 5.3-5.4):

- A fraction ``f_cen`` of clusters are correctly centered; the remainder carry
  a 2D projected offset R_s drawn from a Rayleigh distribution whose scale
  sigma_s = tau * R_lambda(lambda, z) tracks the richness-radius relation
  (Rykoff+2012 / McClintock+2019 Eq. 25).
- The miscentered Sigma(R) is the azimuthal average of the correctly-centered
  profile over the offset-circle, then averaged over the R_s population.
- Delta Sigma(R) is *re-derived* from the miscentered Sigma(R) via the usual
  mean-enclosed-density definition -- it is NOT the same convolution kernel
  applied directly to Delta Sigma_cen (convolution does not commute with the
  enclosed-mass integral).

Default (f_cen, tau) are illustrative DES Y1-redMaPPer-like values (loosely
following McClintock+2019), not calibrated against this Buzzard mock -- same
caveat as the other "known simplification" systematics in this repo.
"""
import numpy as np

# ---------------------------------------------------------------------------
# Illustrative DES Y1 redMaPPer-like defaults (McClintock+2019 Sec. 5.3-5.4),
# NOT fit to this Buzzard mock. One (f_cen, tau) pair per richness bin,
# [20,30), [30,45), [45,60), [60,inf).
# ---------------------------------------------------------------------------
F_CEN_DEFAULT = np.array([0.75, 0.80, 0.85, 0.88])
TAU_DEFAULT = 0.21  # sigma_s = tau * R_lambda


def R_lambda(lam, z, h=0.70):
    """Richness-radius relation (Rykoff+2012 / McClintock+2019 Eq. 25).

    Returns the characteristic cluster radius in physical Mpc.
    """
    lam = np.asarray(lam, dtype=float)
    z = np.asarray(z, dtype=float)
    R_cMpc_h = 1.0 * (lam / 100.0) ** 0.2
    return R_cMpc_h / h / (1.0 + z)


def _loglog_interp_extrap(R_query, R_grid, Y_grid):
    """Log-log linear interpolation with log-log linear extrapolation.

    Unlike ``np.interp`` (which clamps to the edge value outside the grid),
    this extends the nearest-edge power law -- needed because the
    miscentering convolution probes radii below the tabulated grid's R_min.
    """
    R_query = np.asarray(R_query, dtype=float)
    lnR_grid = np.log(R_grid)
    lnY_grid = np.log(Y_grid)
    lnR_q = np.log(R_query)

    lnY_q = np.interp(lnR_q, lnR_grid, lnY_grid)

    below = lnR_q < lnR_grid[0]
    if np.any(below):
        slope_lo = (lnY_grid[1] - lnY_grid[0]) / (lnR_grid[1] - lnR_grid[0])
        lnY_q = np.where(
            below, lnY_grid[0] + slope_lo * (lnR_q - lnR_grid[0]), lnY_q
        )

    above = lnR_q > lnR_grid[-1]
    if np.any(above):
        slope_hi = (lnY_grid[-1] - lnY_grid[-2]) / (lnR_grid[-1] - lnR_grid[-2])
        lnY_q = np.where(
            above, lnY_grid[-1] + slope_hi * (lnR_q - lnR_grid[-1]), lnY_q
        )

    return np.exp(lnY_q)


def _sigma_mis_fixed_offset(R_eval, Rs, sigma_cen_func, n_theta=60):
    """Azimuthally-averaged Sigma(R) for a single fixed offset Rs.

    Johnston+2007 Eq. 32:
        Sigma_mis(R|Rs) = (1/2pi) * int_0^2pi Sigma_cen(sqrt(R^2+Rs^2-2*R*Rs*cos(theta))) dtheta
    Exploits the theta -> 2pi-theta symmetry to integrate over [0, pi] only.
    """
    R_eval = np.atleast_1d(np.asarray(R_eval, dtype=float))
    theta = np.linspace(0.0, np.pi, n_theta)
    Rm = R_eval[:, None]
    cos_t = np.cos(theta)[None, :]
    Rprime = np.sqrt(np.clip(Rm ** 2 + Rs ** 2 - 2.0 * Rm * Rs * cos_t, 0.0, None))
    Sigma_vals = sigma_cen_func(Rprime.ravel()).reshape(Rprime.shape)
    return np.trapz(Sigma_vals, theta, axis=1) / np.pi


def sigma_miscentered_avg(R_eval, sigma_cen_func, sigma_s, n_theta=60, n_Rs=40,
                           Rs_max_sigma=5.0):
    """<Sigma_mis(R)> averaged over the Rayleigh offset population P(Rs)."""
    R_eval = np.atleast_1d(np.asarray(R_eval, dtype=float))
    if sigma_s <= 0:
        return sigma_cen_func(R_eval)

    Rs_grid = np.linspace(1e-6, Rs_max_sigma * sigma_s, n_Rs)
    P_Rs = (Rs_grid / sigma_s ** 2) * np.exp(-0.5 * (Rs_grid / sigma_s) ** 2)

    vals = np.array([
        _sigma_mis_fixed_offset(R_eval, Rs, sigma_cen_func, n_theta=n_theta)
        for Rs in Rs_grid
    ])  # (n_Rs, n_R)

    norm = np.trapz(P_Rs, Rs_grid)
    return np.trapz(vals * P_Rs[:, None], Rs_grid, axis=0) / norm


def _cumtrapz(y, x):
    dx = np.diff(x)
    seg = 0.5 * (y[1:] + y[:-1]) * dx
    return np.concatenate([[0.0], np.cumsum(seg)])


def delta_sigma_from_sigma(R_grid, Sigma_grid):
    """Delta Sigma(R) = Sigma_bar(<R) - Sigma(R) from a tabulated Sigma(R).

    R_grid must start close to 0 (the enclosed-mass integral below R_grid[0]
    is treated as negligible, valid as long as R_grid[0] << R_grid[-1]).
    """
    R_grid = np.asarray(R_grid, dtype=float)
    Sigma_grid = np.asarray(Sigma_grid, dtype=float)
    cum = _cumtrapz(R_grid * Sigma_grid, R_grid)
    Sigma_bar = 2.0 * cum / R_grid ** 2
    return Sigma_bar - Sigma_grid


def apply_miscentering(R_native, Sigma_cen_native, DSigma_cen_native,
                        f_cen, tau, lam_rep, z_rep, R_eval, h=0.70,
                        n_theta=60, n_Rs=40, n_fine=300,
                        R_fine_min=1e-3, R_fine_max=30.0):
    """Mix correctly-centered and miscentered profiles; return obs + biases.

    Returns
    -------
    Sigma_obs, DSigma_obs : (n_R,) arrays on ``R_eval``
    B_mis_Sigma, B_mis_DSigma : (n_R,) multiplicative bias ratios obs/cen
    """
    sigma_cen_func = lambda R: _loglog_interp_extrap(R, R_native, Sigma_cen_native)

    R_eval = np.atleast_1d(np.asarray(R_eval, dtype=float))
    Sigma_cen_eval = sigma_cen_func(R_eval)
    DSigma_cen_eval = _loglog_interp_extrap(R_eval, R_native, DSigma_cen_native)

    if f_cen >= 1.0:
        Sigma_obs = Sigma_cen_eval
        DSigma_obs = DSigma_cen_eval
    else:
        Rlam = R_lambda(lam_rep, z_rep, h)
        sigma_s = tau * Rlam

        R_fine = np.geomspace(R_fine_min, R_fine_max, n_fine)
        Sigma_mis_fine = sigma_miscentered_avg(
            R_fine, sigma_cen_func, sigma_s, n_theta=n_theta, n_Rs=n_Rs,
        )
        DSigma_mis_fine = np.clip(
            delta_sigma_from_sigma(R_fine, Sigma_mis_fine), 1e-30, None,
        )

        Sigma_mis_eval = _loglog_interp_extrap(R_eval, R_fine, Sigma_mis_fine)
        DSigma_mis_eval = _loglog_interp_extrap(R_eval, R_fine, DSigma_mis_fine)

        Sigma_obs = f_cen * Sigma_cen_eval + (1.0 - f_cen) * Sigma_mis_eval
        DSigma_obs = f_cen * DSigma_cen_eval + (1.0 - f_cen) * DSigma_mis_eval

    with np.errstate(divide='ignore', invalid='ignore'):
        B_mis_Sigma = np.where(Sigma_cen_eval != 0,
                                Sigma_obs / Sigma_cen_eval, np.nan)
        B_mis_DSigma = np.where(DSigma_cen_eval != 0,
                                 DSigma_obs / DSigma_cen_eval, np.nan)

    return Sigma_obs, DSigma_obs, B_mis_Sigma, B_mis_DSigma
