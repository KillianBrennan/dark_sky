"""Step 4: light-dome model.

KERNEL. For a unit point source on the ground at distance D and an observer looking at
azimuth offset da, elevation h, the artificial radiance is the single-scattering
line-of-sight integral (Garstang 1986; Cinzano et al. 2000):

  L = int ds [b_m(z) p_m(th) + b_a(z) p_a(th)] * I(psi)/rho^2 * exp(-tau_sp - tau_po)

  b_m: Rayleigh, sea-level 0.0116/km at 550 nm, scale height 8 km, ~lambda^-4
  b_a: aerosol, AOD/H_a, scale height 1.5 km, Angstrom 1.3, Henyey-Greenstein g = 0.65
  I(psi): Garstang emission function, normalised to the nadir value VIIRS measures:
          I(psi)/I(0) = cos(psi) + 0.554 F psi^4 / (2 G (1 - F))
          (the psi^4 term is the near-horizontal emission that VIIRS cannot see)
  Earth curvature for both the line of sight and the source position; source vertical
  tilted by D/R. Extinction along both legs uses the closed form for an exponential
  atmosphere.

TERRAIN BLOCKING (two-component, partial). Each kernel entry is binned by the
elevation angle eps of the scattering point as seen from the source. Terrain between
source and observer defines a shadow angle beta (terrain.py): scattering points with
eps < beta are not directly lit by that source (the *blockable, near-surface*
component), points above are (the *unblockable, high-altitude* component). The split
therefore follows from geometry instead of being imposed: a ridge near the observer
shades the low air near the observer, which dominates low-elevation views towards the
city, while the zenith column above ~1-2 km stays lit. A SHADOW_LEAK fraction of the
shadowed part is kept to stand in for multiple scattering.

The table K[band, D, da, h, dz, eps] stores cumulative contributions up to eps.
"""
import math
import time

import numpy as np
from numba import njit, prange

from . import config as C
from .viirs import FAR_RES, NEAR_RES

D_KM = np.geomspace(0.2, 300.0, 64)
DA_DEG = np.arange(0, 181, C.AZ_BIN)
H_DEG = np.array(C.ELEVATIONS, float)
DZ_M = np.array([-600, -300, 0, 300, 600, 900, 1200, 1600, 2000, 2600], float)  # z_obs - z_src
EPS_DEG = np.array([-2, -1, -0.5, 0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 11, 15, 20, 30, 90], float)
R = C.R_EARTH


def _tau(beta0, H, z1, z2, length):
    dz = z2 - z1
    small = np.abs(dz) < 1.0
    with np.errstate(divide="ignore", invalid="ignore"):
        slant = beta0 * length * H * (np.exp(-z1 / H) - np.exp(-z2 / H)) / np.where(small, 1.0, dz)
    return np.where(small, beta0 * length * np.exp(-0.5 * (z1 + z2) / H), slant)


def build_kernel(aod=C.AOD_550, F=C.GARSTANG_F, bands=tuple(C.BANDS), zenith_only=False, n_s=240):
    """Returns K[band, D, da, h, dz, eps] (cumulative over eps), arbitrary units."""
    hs = np.array([90.0]) if zenith_only else H_DEG
    das = np.array([0.0]) if zenith_only else DA_DEG
    K = np.zeros((len(bands), len(D_KM), len(das), len(hs), len(DZ_M), len(EPS_DEG)), np.float64)
    emis_c = 0.554 * F / (2 * C.GARSTANG_G * (1 - F))
    eps_edges = np.radians(EPS_DEG)
    for ib, band in enumerate(bands):
        lam = C.BANDS[band]
        bm0 = C.RAYLEIGH_BETA0_550 * (550.0 / lam) ** 4
        ba0 = (aod / C.H_AER) * (550.0 / lam) ** C.ANGSTROM
        for ih, h in enumerate(hs):
            hr = math.radians(h)
            smax = min(250e3, 45e3 / max(math.sin(hr), 0.15))
            s = np.geomspace(5.0, smax, n_s)
            ds = np.gradient(s)
            for iz, dz in enumerate(DZ_M):
                zs = C.SOURCE_ALT
                zo = zs + dz
                for ia, da in enumerate(das):
                    ar = math.radians(da)
                    los = np.array([math.cos(hr) * math.cos(ar), math.cos(hr) * math.sin(ar), math.sin(hr)])
                    P = s[:, None] * los[None, :]
                    P[:, 2] += zo
                    zP = P[:, 2] + (s * math.cos(hr)) ** 2 / (2 * R)          # true altitude
                    dens_m = bm0 * np.exp(-zP / C.H_MOL)
                    dens_a = ba0 * np.exp(-zP / C.H_AER)
                    tau_po = _tau(bm0, C.H_MOL, zo, zP, s) + _tau(ba0, C.H_AER, zo, zP, s)
                    for idd, Dkm in enumerate(D_KM):
                        D = Dkm * 1e3
                        S = np.array([D, 0.0, zs - D * D / (2 * R)])
                        v = P - S
                        rho = np.linalg.norm(v, axis=1)
                        tilt = D / R
                        ns = np.array([math.sin(tilt), 0.0, math.cos(tilt)])
                        cpsi = (v @ ns) / rho
                        psi = np.arccos(np.clip(cpsi, -1, 1))
                        I = np.where(cpsi > 0, cpsi + emis_c * psi ** 4, 0.0)
                        cth = -(v @ los) / rho          # d_in . d_out, d_out = -los
                        pm = 3 / (16 * math.pi) * (1 + cth ** 2)
                        g = C.HG_G
                        pa = (1 - g * g) / (4 * math.pi * (1 + g * g - 2 * g * cth) ** 1.5)
                        tau_sp = _tau(bm0, C.H_MOL, zs, zP, rho) + _tau(ba0, C.H_AER, zs, zP, rho)
                        dL = (dens_m * pm + dens_a * pa) * I / rho ** 2 * np.exp(-tau_sp - tau_po) * ds
                        eps = np.pi / 2 - psi
                        # cumulative contribution with eps <= edge
                        idx = np.searchsorted(eps_edges, eps)             # bin index
                        hist = np.bincount(np.minimum(idx, len(EPS_DEG)), weights=dL, minlength=len(EPS_DEG) + 1)
                        cum = np.cumsum(hist)[:len(EPS_DEG)]
                        cum[-1] = hist.sum()                             # last edge = everything
                        K[ib, idd, ia, ih, iz, :] = cum
    return K


def kernel(aod=C.AOD_550, F=C.GARSTANG_F):
    f = C.INTERIM / f"kernel_aod{aod:.3f}_F{F:.2f}.npz"
    if f.exists():
        return np.load(f)["K"]
    t = time.time()
    K = build_kernel(aod, F)
    np.savez_compressed(f, K=K)
    print(f"kernel aod={aod} F={F} built in {time.time() - t:.0f}s")
    return K


# ---------------------------------------------------------------- site summation
@njit(cache=True, fastmath=True)
def _interp_log(xs_log, x):
    lx = math.log(x)
    n = xs_log.shape[0]
    if lx <= xs_log[0]:
        return 0, 0.0
    if lx >= xs_log[n - 1]:
        return n - 2, 1.0
    i = 0
    while xs_log[i + 1] < lx:
        i += 1
    return i, (lx - xs_log[i]) / (xs_log[i + 1] - xs_log[i])


@njit(cache=True, fastmath=True)
def _interp_lin(xs, x):
    n = xs.shape[0]
    if x <= xs[0]:
        return 0, 0.0
    if x >= xs[n - 1]:
        return n - 2, 1.0
    i = 0
    while xs[i + 1] < x:
        i += 1
    return i, (x - xs[i]) / (xs[i + 1] - xs[i])


@njit(cache=True, fastmath=True)
def _accumulate(K, d_log, dz_grid, eps_grid, beta_dlog, beta_row, use_terrain, leak, want, az_bin,
                sx, sy, sz, qx, qy, qz, qw, soft, B, Bb):
    nh = K.shape[2]
    n_az = want.shape[1]
    n_da = K.shape[1]
    ne = eps_grid.shape[0]
    nbaz = beta_row.shape[0]
    dx = qx - sx
    dy = qy - sy
    # softened distance: a source cell is an area, not a point (soft = cell / sqrt(6),
    # the RMS radius of a uniform square); removes spikes at cell centroids
    D = math.sqrt(dx * dx + dy * dy + soft * soft)
    az = math.degrees(math.atan2(dx, dy)) % 360.0
    i_d, f_d = _interp_log(d_log, D / 1000.0)
    i_z, f_z = _interp_lin(dz_grid, sz - qz)
    beta = -90.0
    if use_terrain:
        ib = int(round(az / (360.0 / nbaz))) % nbaz
        j, fj = _interp_log(beta_dlog, D / 1000.0)
        beta = beta_row[ib, j] * (1 - fj) + beta_row[ib, j + 1] * fj
    shadowed = beta > eps_grid[0]
    i_e, f_e = _interp_lin(eps_grid, beta)
    for ih in range(nh):
        for a in range(n_az):
            if not want[ih, a]:
                continue
            da = abs(((a * az_bin) - az + 180.0) % 360.0 - 180.0)   # bin centred on a*az_bin
            ia = int(round(da / az_bin))
            if ia > n_da - 1:
                ia = n_da - 1
            tot = 0.0
            sh = 0.0
            for u in range(2):
                wd = f_d if u else 1 - f_d
                for v in range(2):
                    w = wd * (f_z if v else 1 - f_z)
                    tot += w * K[i_d + u, ia, ih, i_z + v, ne - 1]
                    if shadowed:
                        sh += w * (K[i_d + u, ia, ih, i_z + v, i_e] * (1 - f_e)
                                   + K[i_d + u, ia, ih, i_z + v, i_e + 1] * f_e)
            blocked = (1.0 - leak) * sh
            B[ih, a] += qw * (tot - blocked)
            Bb[ih, a] += qw * blocked


@njit(parallel=True, cache=True, fastmath=True)
def _sum_sites(K, d_log, dz_grid, eps_grid, beta_dlog, site_x, site_y, site_z, beta_tab,
               nx, ny, nz, nw, fx, fy, fz, fw, fstart, fcount, want, az_bin, leak, use_terrain,
               max_d, near_d, soft_near, soft_far):
    """B[site, h, az] artificial radiance (kernel units) and the part removed by terrain."""
    n_site = site_x.shape[0]
    nh = K.shape[2]
    n_az = want.shape[1]
    B = np.zeros((n_site, nh, n_az))
    Bb = np.zeros((n_site, nh, n_az))
    for p in prange(n_site):
        sx = site_x[p]
        sy = site_y[p]
        sz = site_z[p]
        brow = beta_tab[p] if use_terrain else beta_tab[0]
        Bp = np.zeros((nh, n_az))
        Bbp = np.zeros((nh, n_az))
        for q in range(fx.shape[0]):
            dx = fx[q] - sx
            dy = fy[q] - sy
            d = math.sqrt(dx * dx + dy * dy)
            if d > max_d:
                continue
            if d < near_d:
                for k in range(fstart[q], fstart[q] + fcount[q]):
                    _accumulate(K, d_log, dz_grid, eps_grid, beta_dlog, brow, use_terrain, leak, want,
                                az_bin, sx, sy, sz, nx[k], ny[k], nz[k], nw[k], soft_near, Bp, Bbp)
            else:
                _accumulate(K, d_log, dz_grid, eps_grid, beta_dlog, brow, use_terrain, leak, want,
                            az_bin, sx, sy, sz, fx[q], fy[q], fz[q], fw[q], soft_far, Bp, Bbp)
        B[p] = Bp
        Bb[p] = Bbp
    return B, Bb


def site_brightness(K_band, sites_xyz, beta_tab, beta_d_km, near, far, want, use_terrain=True,
                    leak=C.SHADOW_LEAK, max_d_km=C.SOURCE_RADIUS_KM, near_d_km=15.0):
    """K_band: K[band] (D, da, h, dz, eps). `want` boolean [n_h, 360/AZ_BIN] selects the
    (elevation, azimuth-bin) outputs to compute. Returns (B, Bblocked), (n_site, n_h, n_az)."""
    sx, sy, sz = (np.ascontiguousarray(a, dtype=np.float64) for a in sites_xyz)
    bt = np.ascontiguousarray(beta_tab, dtype=np.float32)
    return _sum_sites(np.ascontiguousarray(K_band), np.log(D_KM), DZ_M, EPS_DEG, np.log(beta_d_km),
                      sx, sy, sz, bt, near["x"], near["y"], near["z"], near["w"],
                      far["x"], far["y"], far["z"], far["w"], far["start"], far["count"],
                      np.ascontiguousarray(want), float(C.AZ_BIN), leak, use_terrain,
                      max_d_km * 1e3, near_d_km * 1e3, NEAR_RES / math.sqrt(6), FAR_RES / math.sqrt(6))


def natural_rel(h_deg, k=C.V_EXTINCTION_K):
    """Natural sky brightness at elevation h relative to zenith (van Rhijn-type airmass
    brightening with extinction; Krisciunas & Schaefer 1991)."""
    z = np.radians(90 - np.asarray(h_deg, float))
    X = 1 / np.sqrt(1 - 0.96 * np.sin(z) ** 2)
    return X * 10 ** (-0.4 * k * (X - 1))


def to_mpsas(art_ratio, h_deg=90.0):
    return C.NATURAL_ZENITH_MPSAS - 2.5 * np.log10(natural_rel(h_deg) + art_ratio)
