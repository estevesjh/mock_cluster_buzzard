# mock_cluster_buzzard — Buzzard cluster mock with optical-selection effects

Forward model for a DES Y3-like cluster data vector built on Buzzard halos,
with the Costanzi+2026 (C26) optical-selection systematics layered on. The
notebook constructs the lensing data vector along **two parallel routes**
that live on the same DES Y1 binning, so a downstream user can compare
them directly:

- **Analytical (Matteo / C26 Eq. C1).** Bin halos on the C19 forward-model
  $(\lambda^\mathrm{ob}, z^\mathrm{ob})$ — kernel-driven, no per-halo LSS
  information — then apply the C26 Appendix C fitting function
  $\mathcal{B}^{\Delta\Sigma}_\mathrm{C1}(R)$ to the stacked $\Delta\Sigma$.
  Five constants per richness/redshift bin, no integration. Cheap enough
  for an MCMC chain.
- **Empirical (Heidi / mass-matched ratio).** Bin halos on the catalog
  redMaPPer outputs $(\mathtt{LAMBDA\_CHISQ}, \mathtt{Z\_LAMBDA})$ — these
  carry per-halo LSS information — then divide by a
  $(\log M, z_\mathrm{true})$-matched random reference (Wu+2022 method
  iii). No fitting function, both numerator and denominator come straight
  from `mock['DeltaSigma']`.

Both routes share the same C19 forward model upstream (mass–richness +
projection kernel) and the same lensing geometry $\Sigma_\mathrm{crit}^{-1}(z_l)$.

Member-galaxy contamination is injected from the published DES Y1 boost-factor
profiles in `data/boost_factor/profiles/` (McClintock+2019); the saved final
shear `gamma_t_mock_obs` mimics the raw uncorrected DES Y1 data vector.

Cluster miscentering (optical-finder center offset from the true halo
center) is injected as a third multiplicative systematic, `src/miscentering.py`
(§8.4): a fraction of clusters have their lensing center offset by a 2D
distance drawn from a Rayleigh distribution tied to the richness-radius
relation, following Johnston+2007 / McClintock+2019 §5.3-5.4.

## Outputs

Saved by the final notebook cell to `output/MockDataVector.npz`. With
$N_\lambda = 4$, $N_z = 3$, $N_R = 11$ (DES Y1 boost-factor R grid,
$0.166$ to $15.80$ Mpc):

| key                    | shape          | meaning                                                          |
|------------------------|----------------|------------------------------------------------------------------|
| `scinv_bin_C19`        | $(4, 3)$       | $\Sigma_\mathrm{crit}^{-1}(\langle z^\mathrm{ob}\rangle_\mathrm{bin})$ [$\mathrm{pc}^2/M_\odot$], physical; lets `gamma_t_mock_obs_C1`/`gamma_t_stack_C19`/`gamma_t_obs_C1` be inverted to $\Delta\Sigma$ without replaying the forward model. **C19/analytical route only** — not self-consistent with the empirical/RM route (see units note below) |
| `z_rep_C19`            | $(4, 3)$       | bin-mean $\langle z^\mathrm{ob}\rangle$ used to evaluate `scinv_bin_C19`. C19/analytical route only |
| `h_buzzard`             | scalar         | Buzzard's fixed $h = 0.70$ (DeRose+2019); tag so downstream code never hardcodes it |
| `scinv_bin_C19_littleh`| $(4, 3)$       | $\Sigma_\mathrm{crit}^{-1}$ in the little-h convention, $\mathrm{pc}^2/(h\,M_\odot)$ = `scinv_bin_C19 / h_buzzard`. Dividing the (unchanged, dimensionless) `gamma_t_mock_obs_C1` by this recovers $\Delta\Sigma$ directly in $h\,M_\odot/\mathrm{pc}^2$ |
| `radii_cMpc_h_C19`     | $(4, 3, 11)$   | comoving Mpc/h coordinate of each of the 11 `radii_phys_mpc` points, per bin: `radii_phys_mpc * (1+z_rep_C19) * h_buzzard`. A coordinate relabeling of the same points `gamma_t_mock_obs_C1` is already indexed by — not a re-grid |
| `NC`                   | $(4, 3)$       | $N(\lambda^\mathrm{ob}, z^\mathrm{ob})$, Buzzard counts rescaled by $\Omega_\mathrm{Y1}(z)$ |
| `NC_cov`               | $(12, 12)$     | diagonal Poisson covariance of `NC`                              |
| `gamma_t_stack_C19`    | $(4, 3, 11)$   | analytical: stacked $\gamma_t$, no correction                    |
| `gamma_t_obs_C1`       | $(4, 3, 11)$   | analytical: stack $\times \mathcal{B}^{\Delta\Sigma}_\mathrm{C1}$ |
| `B_sel_C1`             | $(4, 3, 11)$   | analytical: $\mathcal{B}^{\Delta\Sigma}_\mathrm{C1}(R)$          |
| `shear_C19_stack_cov`  | $(4, 3, 11, 11)$| jackknife covariance of `gamma_t_stack_C19` (§8.3)              |
| `gamma_t_stack_RM`     | $(4, 3, 11)$   | empirical: $(\log M, z)$-matched random reference                |
| `gamma_t_obs_RM`       | $(4, 3, 11)$   | empirical: redMaPPer-selected stack                              |
| `B_sel_emp`            | $(4, 3, 11)$   | empirical: bootstrap-mean $\mathcal{B}_\mathrm{sel}^\mathrm{emp}(R)$ |
| `shear_data_vector_cov`| $(4, 3, 11, 11)$| jackknife covariance of `gamma_t_obs_RM` (§8.3)                 |
| `B_data`               | $(4, 3, 11)$   | DES Y1 measured boost factor $B(R)$                              |
| `B_data_err`           | $(4, 3, 11)$   | $1\sigma$ error per radial bin                                  |
| `B_data_cov`           | $(4, 3, 11, 11)$| full radial covariance per bin                                  |
| `B_mis_DSigma`         | $(4, 3, 11)$   | miscentering bias $\mathcal{B}_\mathrm{mis}^{\Delta\Sigma}(R) = \Delta\Sigma_\mathrm{obs}/\Delta\Sigma_\mathrm{cen}$ (§8.4) |
| `B_mis_Sigma`          | $(4, 3, 11)$   | miscentering bias on $\Sigma(R)$, same construction             |
| `f_cen_mis`            | $(4, 3)$       | correctly-centered fraction used per bin (illustrative DES Y1 redMaPPer-like default, not Buzzard-calibrated) |
| `tau_mis`              | scalar         | offset-scale parameter, $\sigma_s = \mathtt{tau\_mis}\times R_\lambda(\lambda,z)$ |
| `gamma_t_mock_obs_C1`  | $(4, 3, 11)$   | analytical route $\times \mathcal{B}^{\Delta\Sigma}_\mathrm{C1} \times \mathcal{B}_\mathrm{mis}^{\Delta\Sigma} / B_\mathrm{data}$ |
| `gamma_t_mock_obs_RM`  | $(4, 3, 11)$   | empirical route $\times \mathcal{B}_\mathrm{mis}^{\Delta\Sigma} / B_\mathrm{data}$ |

Plus bin metadata: `radii_phys_mpc`, `lambda_bins`, `z_bin_min`, `z_bin_max`.

### Units & h-convention

All `gamma_t_*` arrays are dimensionless ($\gamma_t = \Delta\Sigma \times
\Sigma_\mathrm{crit}^{-1}$) computed self-consistently in physical units
throughout — they never need an h-conversion and are identical regardless
of which convention below you read `scinv_bin_C19`/`radii_phys_mpc` in.

- **Physical** (`radii_phys_mpc`, `scinv_bin_C19`): proper Mpc and
  $\mathrm{pc}^2/M_\odot$, at Buzzard's fixed cosmology $h = 0.70$. This is
  an absolute-unit quantity, not a little-h-parametrized one — Buzzard's
  cosmology isn't a free parameter of the mock.
- **Little-h** (`radii_cMpc_h_C19`, `scinv_bin_C19_littleh`): the standard
  weak-lensing convention, comoving Mpc/h and $\mathrm{pc}^2/(h\,M_\odot)$.
  `gamma_t_mock_obs_C1 / scinv_bin_C19_littleh` recovers $\Delta\Sigma$
  directly in $h\,M_\odot/\mathrm{pc}^2$.
- Both `_C19` fields are valid **only** for the analytical/Matteo route
  (`gamma_t_*_C1`/`gamma_t_stack_C19`) — the empirical/RM route
  (`gamma_t_*_RM`) uses its own internal bin-mean redshift (from the
  redMaPPer catalog, §8.2) that isn't currently exposed in the npz.
- `output/MockDataVector_targetR.npz` (a separate, secondary artifact) is
  **not** covered by any of the above and remains in its previous,
  unconverted physical-Mpc convention.

The headline plots are also written to `output/figs/` as PNGs and embedded
in `docs/MockDataVector.pdf`.

## Cosmology & binning (DES Y1)

Cosmology is the Buzzard $\Lambda$CDM (DeRose+2019, arXiv:1901.02401):
`FlatLambdaCDM(H0=70, Om0=0.286, Ob0=0.046, Tcmb0=2.725)`, with $\sigma_8 = 0.82$,
$n_s = 0.96$ passed to `hmf` in §6.

Counts and lensing share the DES Y1 cluster binning (Costanzi+2019,
McClintock+2019):

- `LBDBINS = [20, 30, 45, 60, 500]` (4 richness bins; `500` caps $[60, \infty)$).
- `ZMIN_LIST = [0.20, 0.37, 0.50]`, `ZMAX_LIST = [0.33, 0.50, 0.65]` — the
  inner edges are pulled to 0.33/0.37 (from a nominal 0.35) to skip the Buzzard
  $0.33 \le z \le 0.37$ box-junction seam, dropped at halo-load time.

Halo selection mirrors `0-MakeMock.ipynb`: `pid == -1`, $0 \le \cos i \le 1$,
drop the Buzzard $0.33 \le z < 0.37$ box seam, then require
$0.20 \le z \le 0.65$ and $\log_{10} M_\mathrm{vir} \ge 13$. Yields
$\sim 5.83 \times 10^5$ halos.

## DES-Y1 NC+3×2pt mass–richness parameters

Hard-coded in `src/costanzi_selection.py`.

| Symbol     | Value             |
|------------|-------------------|
| `M_min`    | $10^{11.385}\,M_\odot$ |
| `α`        | 0.859             |
| `M_1`      | $10^{12.696}\,M_\odot$ |
| `M_pivot`  | $M_1 - M_\mathrm{min}$ |
| `σ_intr`   | 0.181             |
| `ε`        | 0.284             |
| `z_piv`    | 0.4544            |

Projection coefficients ($a_\tau, b_\tau, a_{f_\mathrm{prj}}, b_{f_\mathrm{prj}}, \dots$)
are loaded as the posterior mean of the 15 rows in
`data/prj_params_DESY3_lss_lin_dep_getdist_v1.txt` (downloaded from
`MCostanzi/SelectionBias`).

## Notebook structure

| Section | Content |
|---|---|
| §1 | Inputs (cosmology, DES Y1 binning, RNG seed) |
| §2 | Load halo catalogue + per-halo $\Sigma$, $\Delta\Sigma$ profiles; photo-$z$ proxy |
| §3 | Sample $\lambda_\mathrm{true}$ from C26 Eq. 15 |
| §4 | Project $\lambda_\mathrm{true} \to \lambda_\mathrm{obs}$ via C19 Eq. 6 mixture |
| §4.5 | Stacked $\Sigma$ in $\lambda^\mathrm{tr}$ vs $\lambda^\mathrm{ob}$ bins (mass-dilution diagnostic) |
| §5 | Output 1 — $N(\lambda_\mathrm{obs}, z)$ |
| §5.5 | Sanity: C19 sampler $\lambda^\mathrm{ob}$ vs catalog `LAMBDA_CHISQ` |
| §6 | Halo mass function: Buzzard vs analytical Tinker08 (Bryan–Norman virial) |
| **§7** | **Lensing data vector — ANALYTICAL route (C26 Eq. C1)** |
| §7.1 | Selection-bias ratio $\mathcal{B}^{\Delta\Sigma}_\mathrm{C1}(R)$ |
| §7.2 | Tangential shear, with vs without C1 correction |
| §7.3 | Boost factor $B(R)$ — McClintock+2019 parametric fit (sanity) |
| §7.3.bis | Boost factor $B(R)$ — measured DES Y1 profiles; **saved** mock-observed shear |
| **§8** | **Lensing data vector — EMPIRICAL route (mass-matched ratio)** |
| §8.1 | Selection-bias ratio $\mathcal{B}_\mathrm{sel}^\mathrm{emp}(R)$ with bootstrap bands; null test |
| §8.2 | Tangential shear: redMaPPer-selected vs mass-matched random |
| §8.3 | Shear-stack covariance via spatial jackknife ($K=50$ KMeans patches) |
| §8.4 | Miscentering injection — $\mathcal{B}_\mathrm{mis}^{\Delta\Sigma}(R)$ (`src/miscentering.py`) |
| §9 | Save the data vector to `output/MockDataVector.npz` |

## How to run

1. Open `MockDataVector.ipynb` (at the repo root) and run top to
   bottom. The first cell adds `src/` to `sys.path` so the helper
   modules import directly. Expected wall-clock time: a few minutes
   (the §8.1 bootstrap of the mass-matched reference is the slow part).
2. Outputs land in `output/MockDataVector.npz` plus headline plots
   in `output/figs/`.

## Repository layout

```
mock_cluster_buzzard/
├── README.md
├── MockDataVector.ipynb        ← main entry, run top to bottom
├── src/                            ← Python helpers (added to sys.path by the notebook)
│   ├── costanzi_selection.py       ← C26 Eq. 15 sampler + C19 Eq. 6 mixture
│   ├── stacked_profile_weighted_by_mass_redshift.py  ← Wu+22 (M, z) match
│   ├── radial_bins_phys_mpc.py     ← rp_phys_mpc 15-bin radial grid
│   ├── miscentering.py             ← offset-convolution miscentering model (§8.4)
│   └── fileLoc.py                  ← machine-resolved paths to Buzzard inputs
├── data/
│   └── prj_params_DESY3_lss_lin_dep_getdist_v1.txt   ← C19 projection coefs
├── docs/
│   ├── MockDataVector.tex      ← technical reference
│   └── MockDataVector.pdf      ← compiled PDF (embeds vital plots)
└── output/                         ← produced by running the notebook
    ├── MockDataVector.npz
    └── figs/                       ← headline PNGs (7 figures)
```

## Files the notebook reads

External Buzzard inputs, all resolved via `FileLocs(machine='nersc')` in
`src/fileLoc.py`:

- `halo_run_fname` — raw halo catalogue (Mvir, z, RA, DEC, pid, cosi, …).
- `profile_output_fname` — per-halo $\Sigma(R)$ and $\Delta\Sigma(R)$ profiles.
- `mock_boost_factor_1d` → `beta_table_zl_y1_like.npz` ($\beta_\mathrm{eff}$,
  $z_\mathrm{lens}$ grid used to build $\Sigma_\mathrm{crit}^{-1}$).

Repo-local:

- `data/prj_params_DESY3_lss_lin_dep_getdist_v1.txt` — Costanzi posterior-sample
  table for the projection model (15 rows, 10 coefficients).
- `data/boost_factor/profiles/full-unblind-v2-mcal-zmix_y1clust_l{l}_z{z}_zpdf_boost{,_cov}.dat`
  — DES Y1 measured boost factor $B(R)$ profiles + radial covariance per
  $(\lambda, z)$ bin. Used by §7.3 of the notebook.

## Known simplifications / deferred work

- **Boost factor uses DES Y1, not Y3, measurements.** §7.3 injects
  $B(R)$ from McClintock+2019 Y1 profiles, the closest released
  measurement matching our cluster binning. A Y3 boost-factor
  measurement is not yet released; when it is, swap files in
  `data/boost_factor/profiles/`.
- Projection parameters use the **posterior mean** of the C19 chain;
  row-by-row covariance propagation is not done.
- The photo-$z$ kernel is a Gaussian $\sigma_z = 0.01\,(1+z_\mathrm{tr})$
  proxy rather than the full DES Y3 redshift-kernel of Myles+2021.
- Masking fraction $f_\mathrm{msk}$ (an area-masking systematic, loaded
  from the C19 projection file as `a_fmsk`/`b_fmsk`) is still unused.
  Centering mis-identification is now injected (§8.4), but its
  $(f_\mathrm{cen}, \tau)$ parameters are illustrative DES Y1
  redMaPPer-like defaults (`F_CEN_DEFAULT`, `TAU_DEFAULT` in
  `src/miscentering.py`), not calibrated against this Buzzard mock.
- The §8 empirical ratio is computed on **one** Y3 realisation; Wu+2022
  Fig. 2 averages 12 Y1 + 1 Y3 realisations and is correspondingly
  tighter.

## References

- Costanzi et al. 2026, *Forward analytical model for the optical selection
  bias on galaxy cluster lensing profiles* (C26). Eq. 15 (mass–richness),
  Appendix C Eq. C1 (selection-bias fit), Sec. III (Poisson–Gaussian
  convolution).
- Costanzi et al. 2019, arXiv:1807.07072. Eq. 6 (Gaussian + EMG mixture),
  Appendix A Eq. A8 (per-parameter forms).
- Wu et al. 2022, arXiv:2203.05416. Fig. 2 (selection-bias grid), Appendix B
  (mass–redshift weighted reference).
- DeRose et al. 2019, arXiv:1901.02401. Buzzard cosmology and catalog.
- McClintock et al. 2019, arXiv:1805.00039. DES Y1 boost-factor profiles.
- `MCostanzi/SelectionBias` — upstream reference notebook + projection-param file.
- `0-MakeMock.ipynb` — Buzzard halo selection recipe (mirrored).
