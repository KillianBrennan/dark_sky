"""All-sky (> 30 deg) brightness with local sources: two-scale light-dome model.

Sources: local_sources.build() (VIIRS downscaled to BFS hectares + fill for VIIRS-dark
inhabited hectares; VIIRS unchanged abroad), organised in three levels:
hectare points (soft 41 m) -> 1 km cells (soft 408 m) -> 5 km blocks (soft 2041 m).

FAR term (250 m grid, terrain shadow from the ray-cast beta table): every source, with
the hectare points nearer than LOCAL_RADIUS_M down-weighted by a smoothstep taper
w_far(D) (0 below R - T, 1 above R + T).
LOCAL term (100 m grid): hectare points with weight 1 - w_far(D), each with its own
terrain shadow angle traced along the source -> observer path in the DEM (lamp height
5 m, first 60 m skipped so a lamp is not shadowed by its own roof). The two terms add
up to exactly one copy of every source.

Sky sampling: zenith plus rings at 30, 45 and 60 deg elevation (12 azimuths each).
Mean over the sky above 30 deg uses solid-angle weights of the bands
[30, 37.5], [37.5, 52.5], [52.5, 75], [75, 90] deg -> 0.218, 0.369, 0.345, 0.068.
"""
import math

import numpy as np
from numba import njit, prange

from . import config as C
from . import lightdome as L
from .terrain import load_dem

HS = np.array([30.0, 45.0, 60.0, 90.0])
D_KM = np.geomspace(0.03, 300.0, 80)
RING_AZ_BINS = np.arange(0, 72, 6)              # every 30 deg (5-deg bins)
BAND_W = np.array([0.2176, 0.3692, 0.3450, 0.0682])
PT_RES, CELL_RES, BLOCK_RES = 100.0, 1000.0, 5000.0
VIIRS_PT_RES = 385.0     # sqrt(460 m x 320 m): a VIIRS pixel kept as a point source
LAMP_H = 5.0


def want_list():
    """(h index, az bin) pairs: zenith once, then each ring."""
    wh, wa = [3], [0]
    for ih in range(3):
        for a in RING_AZ_BINS:
            wh.append(ih)
            wa.append(int(a))
    return np.array(wh, np.int64), np.array(wa, np.int64)


def kernel(aod, F=C.GARSTANG_F):
    f = C.INTERIM / f"kernel_allsky_aod{aod:.3f}_F{F:.2f}.npz"
    if f.exists():
        return np.load(f)["K"]
    old = L.D_KM, L.H_DEG
    L.D_KM, L.H_DEG = D_KM, HS
    try:
        K = L.build_kernel(aod, F, bands=("V",))[0]
    finally:
        L.D_KM, L.H_DEG = old
    np.savez_compressed(f, K=K)
    return K


# ------------------------------------------------------------------ source hierarchy
def hierarchy(src, z_dem):
    """Sort points into cells and cells into blocks; returns flat arrays for numba."""
    z, tr = z_dem
    x, y, w = src["x"], src["y"], src["w"]
    ps = np.where(src["kind"] == 0, VIIRS_PT_RES, PT_RES) / math.sqrt(6.0)
    c = ((x - tr.c) / tr.a).astype(np.int64)
    r = ((y - tr.f) / tr.e).astype(np.int64)
    ok = (r >= 0) & (r < z.shape[0]) & (c >= 0) & (c < z.shape[1])
    pz = np.full(x.shape, C.SOURCE_ALT)
    pz[ok] = np.nan_to_num(z[r[ok], c[ok]], nan=C.SOURCE_ALT)
    # points are hectares (Swiss) or VIIRS pixel centres (abroad / no proxy)
    ck = np.floor(x / CELL_RES).astype(np.int64) * 10_000_000 + np.floor(y / CELL_RES).astype(np.int64)
    bk = (np.floor(x / BLOCK_RES).astype(np.int64) * 10_000_000 + np.floor(y / BLOCK_RES).astype(np.int64))
    order = np.lexsort((ck, bk))
    x, y, pz, w, ps, ck, bk = x[order], y[order], pz[order], w[order], ps[order], ck[order], bk[order]
    # cells
    uc, cstart, ccount = np.unique(ck, return_index=True, return_counts=True)
    order_c = np.argsort(cstart)
    uc, cstart, ccount = uc[order_c], cstart[order_c], ccount[order_c]
    agg = lambda v: np.add.reduceat(v, cstart)
    cw = agg(w)
    cx, cy, cz = agg(w * x) / cw, agg(w * y) / cw, agg(w * pz) / cw
    cb = bk[cstart]
    ub, bstart, bcount = np.unique(cb, return_index=True, return_counts=True)
    order_b = np.argsort(bstart)
    ub, bstart, bcount = ub[order_b], bstart[order_b], bcount[order_b]
    aggb = lambda v: np.add.reduceat(v, bstart)
    bw = aggb(cw)
    bx, by, bz = aggb(cw * cx) / bw, aggb(cw * cy) / bw, aggb(cw * cz) / bw
    return dict(px=x, py=y, pz=pz, pw=w, ps=ps, cx=cx, cy=cy, cz=cz, cw=cw, cstart=cstart, ccount=ccount,
                bx=bx, by=by, bz=bz, bw=bw, bstart=bstart, bcount=bcount)


# ------------------------------------------------------------------ numba core
@njit(cache=True, fastmath=True)
def _w_far(d, r_loc, taper):
    if r_loc <= 0.0:
        return 1.0
    t = (d - (r_loc - taper)) / (2.0 * taper)
    if t <= 0.0:
        return 0.0
    if t >= 1.0:
        return 1.0
    return t * t * (3.0 - 2.0 * t)


@njit(cache=True, fastmath=True)
def _add(K, d_log, dz_grid, eps_grid, wh, wa, az_bin, sx, sy, sz, qx, qy, qz, qw, soft, beta, leak, out):
    dx = qx - sx
    dy = qy - sy
    D = math.sqrt(dx * dx + dy * dy + soft * soft)
    az = math.degrees(math.atan2(dx, dy)) % 360.0
    i_d, f_d = L._interp_log(d_log, D / 1000.0)
    i_z, f_z = L._interp_lin(dz_grid, sz - qz)
    ne = eps_grid.shape[0]
    shadowed = beta > eps_grid[0]
    i_e, f_e = L._interp_lin(eps_grid, beta)
    n_da = K.shape[1]
    for k in range(wh.shape[0]):
        ih = wh[k]
        if ih == K.shape[2] - 1:       # zenith: azimuth irrelevant
            ia = 0
        else:
            da = abs(((wa[k] * az_bin) - az + 180.0) % 360.0 - 180.0)
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
        out[k] += qw * (tot - (1.0 - leak) * sh)


@njit(cache=True, fastmath=True)
def _beta_from_table(brow, beta_dlog, az, D):
    nbaz = brow.shape[0]
    ib = int(round(az / (360.0 / nbaz))) % nbaz
    j, fj = L._interp_log(beta_dlog, D / 1000.0)
    return brow[ib, j] * (1 - fj) + brow[ib, j + 1] * fj


@njit(parallel=True, cache=True, fastmath=True)
def _far(K, d_log, dz_grid, eps_grid, wh, wa, az_bin, leak, use_terrain, sx, sy, sz, beta_tab, beta_dlog,
         px, py, pz, pw, ps, cx, cy, cz, cw, cstart, ccount, bx, by, bz, bw, bstart, bcount,
         max_d, near_d, r_loc, taper):
    n = sx.shape[0]
    nw = wh.shape[0]
    B = np.zeros((n, nw), np.float32)
    s_c, s_b = CELL_RES / math.sqrt(6.0), BLOCK_RES / math.sqrt(6.0)
    cell_reach = r_loc + taper + 1.42 * CELL_RES     # weighted centroid can sit near a cell corner
    for p in prange(n):
        acc = np.zeros(nw)
        brow = beta_tab[p] if use_terrain else beta_tab[0]
        for b in range(bx.shape[0]):
            d = math.hypot(bx[b] - sx[p], by[b] - sy[p])
            if d > max_d:
                continue
            if d >= near_d:
                beta = _beta_from_table(brow, beta_dlog, math.degrees(math.atan2(bx[b] - sx[p], by[b] - sy[p])) % 360.0, d) if use_terrain else -90.0
                _add(K, d_log, dz_grid, eps_grid, wh, wa, az_bin, sx[p], sy[p], sz[p], bx[b], by[b], bz[b], bw[b], s_b, beta, leak, acc)
                continue
            for c in range(bstart[b], bstart[b] + bcount[b]):
                dc = math.hypot(cx[c] - sx[p], cy[c] - sy[p])
                if r_loc > 0.0 and dc < cell_reach:
                    for q in range(cstart[c], cstart[c] + ccount[c]):
                        dq = math.hypot(px[q] - sx[p], py[q] - sy[p])
                        wf = _w_far(dq, r_loc, taper)
                        if wf <= 0.0:
                            continue
                        beta = _beta_from_table(brow, beta_dlog, math.degrees(math.atan2(px[q] - sx[p], py[q] - sy[p])) % 360.0, max(dq, 200.0)) if use_terrain else -90.0
                        _add(K, d_log, dz_grid, eps_grid, wh, wa, az_bin, sx[p], sy[p], sz[p], px[q], py[q], pz[q], pw[q] * wf, ps[q], beta, leak, acc)
                else:
                    beta = _beta_from_table(brow, beta_dlog, math.degrees(math.atan2(cx[c] - sx[p], cy[c] - sy[p])) % 360.0, max(dc, 200.0)) if use_terrain else -90.0
                    _add(K, d_log, dz_grid, eps_grid, wh, wa, az_bin, sx[p], sy[p], sz[p], cx[c], cy[c], cz[c], cw[c], s_c, beta, leak, acc)
        for k in range(nw):
            B[p, k] = acc[k]
    return B


@njit(cache=True, fastmath=True)
def _bilin(z, fx, fy):
    ny, nx = z.shape
    if fx < 0 or fy < 0 or fx >= nx - 1 or fy >= ny - 1:
        return np.nan
    i = int(fy)
    j = int(fx)
    dy = fy - i
    dx = fx - j
    return ((z[i, j] * (1 - dx) + z[i, j + 1] * dx) * (1 - dy) + (z[i + 1, j] * (1 - dx) + z[i + 1, j + 1] * dx) * dy)


@njit(parallel=True, cache=True, fastmath=True)
def _local(K, d_log, dz_grid, eps_grid, wh, wa, az_bin, leak, sx, sy, sz, zdem, x0, y0, res,
           px, py, pz, pw, ps, bkt_start, bkt_count, bkt_x0, bkt_y0, bkt_res, nbx, nby, r_loc, taper, re_eff):
    n = sx.shape[0]
    nw = wh.shape[0]
    B = np.zeros((n, nw), np.float32)
    reach = r_loc + taper
    nb = int(math.ceil(reach / bkt_res))
    for p in prange(n):
        acc = np.zeros(nw)
        bi = int((sx[p] - bkt_x0) // bkt_res)
        bj = int((sy[p] - bkt_y0) // bkt_res)
        for jj in range(max(bj - nb, 0), min(bj + nb + 1, nby)):
            for ii in range(max(bi - nb, 0), min(bi + nb + 1, nbx)):
                k0 = bkt_start[jj, ii]
                for q in range(k0, k0 + bkt_count[jj, ii]):
                    dq = math.hypot(px[q] - sx[p], py[q] - sy[p])
                    if dq >= reach:
                        continue
                    wl = 1.0 - _w_far(dq, r_loc, taper)
                    if wl <= 0.0:
                        continue
                    # terrain shadow along the lamp -> observer path
                    beta = -90.0
                    if dq >= 120.0:
                        ux = (sx[p] - px[q]) / dq
                        uy = (sy[p] - py[q]) / dq
                        zl = pz[q] + LAMP_H
                        best = -1e9
                        s = 60.0
                        while s < dq - 30.0:
                            fx = (px[q] + ux * s - x0) / res - 0.5
                            fy = (y0 - (py[q] + uy * s)) / res - 0.5
                            zz = _bilin(zdem, fx, fy)
                            if not math.isnan(zz):
                                t = (zz - s * s / (2.0 * re_eff) - zl) / s
                                if t > best:
                                    best = t
                            s += res
                        if best > -1e8:
                            beta = math.degrees(math.atan(best))
                    _add(K, d_log, dz_grid, eps_grid, wh, wa, az_bin, sx[p], sy[p], sz[p], px[q], py[q], pz[q] + LAMP_H,
                         pw[q] * wl, ps[q], beta, leak, acc)
        for k in range(nw):
            B[p, k] = acc[k]
    return B


def buckets(px, py, res=500.0):
    x0, y0 = px.min() - 1, py.min() - 1
    nbx = int((px.max() - x0) // res) + 1
    nby = int((py.max() - y0) // res) + 1
    bi = ((px - x0) // res).astype(np.int64)
    bj = ((py - y0) // res).astype(np.int64)
    key = bj * nbx + bi
    order = np.argsort(key, kind="stable")
    counts = np.bincount(key, minlength=nbx * nby).reshape(nby, nbx)
    start = np.concatenate([[0], np.cumsum(counts.ravel())[:-1]]).reshape(nby, nbx)
    return order, start.astype(np.int64), counts.astype(np.int64), x0, y0, nbx, nby


# ------------------------------------------------------------------ drivers
def far_term(K, H, sites, beta_tab, beta_d_km, use_terrain=True, r_loc=C.LOCAL_RADIUS_M, leak=C.SHADOW_LEAK):
    wh, wa = want_list()
    sx, sy, sz = (np.ascontiguousarray(a, np.float64) for a in sites)
    return _far(np.ascontiguousarray(K), np.log(D_KM), L.DZ_M, L.EPS_DEG, wh, wa, float(C.AZ_BIN), leak, use_terrain,
                sx, sy, sz, np.ascontiguousarray(beta_tab, np.float32), np.log(beta_d_km),
                H["px"], H["py"], H["pz"], H["pw"], H["ps"], H["cx"], H["cy"], H["cz"], H["cw"], H["cstart"], H["ccount"],
                H["bx"], H["by"], H["bz"], H["bw"], H["bstart"], H["bcount"],
                C.SOURCE_RADIUS_KM * 1e3, 15e3, r_loc, C.LOCAL_TAPER_M)


def local_term(K, H, sites, leak=C.SHADOW_LEAK):
    wh, wa = want_list()
    z, tr = load_dem()
    # only hectare-scale points matter locally; VIIRS points abroad are included too (harmless)
    order, start, count, bx0, by0, nbx, nby = buckets(H["px"], H["py"])
    px, py, pz, pw, ps = (np.ascontiguousarray(H[k][order]) for k in ("px", "py", "pz", "pw", "ps"))
    sx, sy, sz = (np.ascontiguousarray(a, np.float64) for a in sites)
    return _local(np.ascontiguousarray(K), np.log(D_KM), L.DZ_M, L.EPS_DEG, wh, wa, float(C.AZ_BIN), leak,
                  sx, sy, sz, z, tr.c, tr.f, tr.a, px, py, pz, pw, ps, start, count, bx0, by0, 500.0, nbx, nby,
                  C.LOCAL_RADIUS_M, C.LOCAL_TAPER_M, C.K_REFRACTION * C.R_EARTH)
