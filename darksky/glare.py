"""Direct glare from nearby lamps as an equivalent veiling luminance.

A lamp in view scatters light inside the eye and lays a veil over the retinal image that
acts like extra background luminance (CIE 146:2002 general disability-glare equation):

  L_veil / E_eye = 10/th^3 + (5/th^2 + 0.1 p/th) [1 + (A/62.5)^4] + 0.0025 p   [1/sr, th in deg]

E_eye from each lamp is taken in the sky model's own calibration, so L_veil comes out in
units of the natural zenith sky (like the artificial sky-glow ratio) with no absolute
photometry:
  E_eye / L_nat = scale * w * g(psi) * T / d^2
  scale  = fitted model scale (calibration_allsky.json) x LED factor (V band)
  w      = source weight (local_sources: hectare or VIIRS point)
  g(psi) = Garstang emission relative to the nadir value VIIRS sees,
           cos(psi) + 0.554 F psi^4 / (2 G (1 - F)), evaluated at the emission angle psi
           towards the eye (psi = 90 deg + elevation of the lamp seen from the eye; eyes
           below a lamp use the horizontal value). The psi^4 term is the near-horizontal
           light VIIRS cannot see.
  T      = extinction along the path (Rayleigh + fitted aerosol at the path's mean height)
Line of sight: per 100 m pixel, the running-maximum terrain elevation angle (Copernicus
DSM, so trees and buildings block) along 72 rays at GLARE_LOS_D_KM distances, stored as
int8 in 0.25 deg steps (-32 .. +31.75 deg, clamped; ~1.5 GB for the 100 m grid); a lamp
5 m above the ground counts with a soft visibility ramp of +/-0.25 deg around that angle.
Gaze: the veil is evaluated for the same sky directions as the all-sky layer (rings at
30/45/60 deg, 12 azimuths, plus zenith) so it adds directly to the sky brightness there.
"""
import math

import numpy as np
from numba import njit, prange

from . import allsky as A
from . import config as C
from .terrain import _bilinear, load_dem

LOS_D = np.array(C.GLARE_LOS_D_KM) * 1000.0
N_AZ = 72
LOS_STEP = 0.25   # deg per int8 count: covers -32 .. +31.75 deg (clamped), finer than the DEM


@njit(cache=True, fastmath=True)
def _q(deg, step):
    v = int(round(deg / step))
    return max(-128, min(127, v))


@njit(parallel=True, cache=True, fastmath=True)
def _los_table(z, x0, y0, res, xs, ys, zeye, nodes, n_az, re_eff, step):
    n = xs.shape[0]
    nd = nodes.shape[0]
    out = np.empty((n, n_az, nd), np.int8)          # units of `step` degrees
    dmax = nodes[nd - 1]
    for p in prange(n):
        for a in range(n_az):
            az = 2.0 * math.pi * a / n_az
            sa = math.sin(az)
            ca = math.cos(az)
            best = -1e9
            s = res
            k = 0
            while k < nd:
                # record the running max for every node already passed
                while k < nd and s > nodes[k]:
                    out[p, a, k] = _q(math.degrees(math.atan(best)) if best > -1e8 else -90.0, step)
                    k += 1
                if s > dmax:
                    break
                fx = (xs[p] + s * sa - x0) / res - 0.5
                fy = (y0 - (ys[p] + s * ca)) / res - 0.5
                zz = _bilinear(z, fx, fy)
                if not math.isnan(zz):
                    t = (zz - s * s / (2.0 * re_eff) - zeye[p]) / s
                    if t > best:
                        best = t
                s += max(res, s / 150.0)
            while k < nd:
                out[p, a, k] = _q(math.degrees(math.atan(best)) if best > -1e8 else -90.0, step)
                k += 1
    return out


def los_table(xs, ys, z0):
    z, tr = load_dem()
    return _los_table(z, tr.c, tr.f, tr.a, np.asarray(xs, float), np.asarray(ys, float),
                      np.asarray(z0, float) + 1.6, LOS_D, N_AZ, C.K_REFRACTION * C.R_EARTH, LOS_STEP)


@njit(cache=True, fastmath=True)
def _cie(theta_deg, age, p):
    th = max(theta_deg, 1.0)
    return 10.0 / th ** 3 + (5.0 / th ** 2 + 0.1 * p / th) * (1.0 + (age / 62.5) ** 4) + 0.0025 * p


@njit(parallel=True, cache=True, fastmath=True)
def _glare(sx, sy, sz, los, nodes_log, px, py, pz, pw, bkt_start, bkt_count, bkt_x0, bkt_y0, bkt_res, nbx, nby,
           radius, scale, kappa, beta_m0, h_m, beta_a0, h_a, re_eff, gaze_h, gaze_az, age, pig, los_step):
    n = sx.shape[0]
    ng = gaze_h.shape[0]
    n_az = los.shape[1]
    nd = nodes_log.shape[0]
    V = np.zeros((n, ng), np.float32)
    E = np.zeros(n, np.float32)
    nb = int(math.ceil(radius / bkt_res))
    sgh = np.sin(gaze_h)
    cgh = np.cos(gaze_h)
    for p in prange(n):
        acc = np.zeros(ng)
        etot = 0.0
        zeye = sz[p] + 1.6
        bi = int((sx[p] - bkt_x0) // bkt_res)
        bj = int((sy[p] - bkt_y0) // bkt_res)
        for jj in range(max(bj - nb, 0), min(bj + nb + 1, nby)):
            for ii in range(max(bi - nb, 0), min(bi + nb + 1, nbx)):
                k0 = bkt_start[jj, ii]
                for q in range(k0, k0 + bkt_count[jj, ii]):
                    dx = px[q] - sx[p]
                    dy = py[q] - sy[p]
                    d = math.sqrt(dx * dx + dy * dy)
                    if d >= radius or d < 100.0:
                        continue
                    zl = pz[q] + 5.0
                    alpha = math.atan((zl - zeye - d * d / (2.0 * re_eff)) / d)      # lamp elevation from eye
                    az = math.atan2(dx, dy) % (2.0 * math.pi)
                    # obstruction angle just short of the lamp: interpolate in azimuth and log-distance
                    fa = az / (2.0 * math.pi) * n_az
                    a0 = int(fa) % n_az
                    a1 = (a0 + 1) % n_az
                    wa = fa - int(fa)
                    ld = math.log(0.97 * d)
                    if ld <= nodes_log[0]:
                        k = 0
                        wk = 0.0
                    elif ld >= nodes_log[nd - 1]:
                        k = nd - 2
                        wk = 1.0
                    else:
                        k = 0
                        while nodes_log[k + 1] < ld:
                            k += 1
                        wk = (ld - nodes_log[k]) / (nodes_log[k + 1] - nodes_log[k])
                    obs0 = los[p, a0, k] * (1 - wk) + los[p, a0, k + 1] * wk
                    obs1 = los[p, a1, k] * (1 - wk) + los[p, a1, k + 1] * wk
                    obs = math.radians(los_step * (obs0 * (1 - wa) + obs1 * wa))
                    vis = (alpha - obs) / math.radians(0.5) + 0.5
                    if vis <= 0.0:
                        continue
                    if vis > 1.0:
                        vis = 1.0
                    # emission angle from the lamp towards the eye
                    psi = 0.5 * math.pi + alpha
                    if psi > 0.5 * math.pi:
                        psi = 0.5 * math.pi
                    g = math.cos(psi) + kappa * psi ** 4
                    zm = 0.5 * (zl + zeye)
                    T = math.exp(-(beta_m0 * math.exp(-zm / h_m) + beta_a0 * math.exp(-zm / h_a)) * d)
                    e = scale * pw[q] * g * T * vis / (d * d)
                    etot += e
                    sa = math.sin(alpha)
                    ca = math.cos(alpha)
                    for k2 in range(ng):
                        cth = sgh[k2] * sa + cgh[k2] * ca * math.cos(gaze_az[k2] - az)
                        if cth > 1.0:
                            cth = 1.0
                        th = math.degrees(math.acos(cth))
                        acc[k2] += e * _cie(th, age, pig)
        for k2 in range(ng):
            V[p, k2] = acc[k2]
        E[p] = etot
    return V, E


def veil(H, sites, los, scale, aod, F=C.GARSTANG_F, age=C.CIE_AGE, pig=C.CIE_PIGMENT):
    """Veiling luminance (natural-zenith units) for the all-sky gaze directions, and the
    total lamp illuminance at the eye (same units x sr)."""
    order, start, count, bx0, by0, nbx, nby = A.buckets(H["px"], H["py"])
    px, py, pz, pw = (np.ascontiguousarray(H[k][order]) for k in ("px", "py", "pz", "pw"))
    wh, wa = A.want_list()
    gaze_h = np.radians(A.HS[wh])
    gaze_az = np.radians(wa * C.AZ_BIN)
    kappa = 0.554 * F / (2 * C.GARSTANG_G * (1 - F))
    sx, sy, sz = (np.ascontiguousarray(a, np.float64) for a in sites)
    return _glare(sx, sy, sz, los, np.log(LOS_D), px, py, pz, pw, start, count, bx0, by0, 500.0, nbx, nby,
                  C.GLARE_RADIUS_M, scale, kappa, C.RAYLEIGH_BETA0_550, C.H_MOL, aod / C.H_AER, C.H_AER,
                  C.K_REFRACTION * C.R_EARTH, gaze_h, gaze_az, age, pig, LOS_STEP)
