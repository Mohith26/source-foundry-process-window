// Adaptive process-window qualification: core module (no dependencies).
// Bayesian Bossung models per feature, a posterior over depth of focus at a
// required exposure latitude, and metrology-site selection that targets the
// uncertainty in that number. Runs in the browser (ES module) and in Node.

// ---------------------------------------------------------------- constants
export const SPEC = { target: 40, tol: 0.10 };            // CD 40 nm +/- 10%
export const EREF = 31.5;                                   // mJ/cm2, dose normalisation
export const USCALE = 0.2;                                  // u = ln(E/EREF)/USCALE
export const FSCALE = 100;                                  // v = F/FSCALE (nm)
export const FOCUS_COLS = [-120, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 120];
export const DOSE_ROWS = [26.0, 27.1, 28.2, 29.3, 30.4, 31.5, 32.6, 33.7, 34.8, 35.9, 37.0];
export const DIE = { w: 21, h: 24, waferR: 150, edge: 3 }; // mm
export const FEATURES = [
  { id: 'dense', name: 'Dense L/S', detail: '40 nm lines, 80 nm pitch' },
  { id: 'semi', name: 'Semi-dense', detail: '40 nm lines, 120 nm pitch' },
  { id: 'iso', name: 'Isolated line', detail: '40 nm line, >400 nm space' },
];
export const BASIS_NAMES = ['1', 'u', 'u2', 'v', 'v2', 'v3', 'v4', 'uv', 'uv2', 'uv3', 'u2v2', 'v5', 'v6', 'uv4'];
export const P = BASIS_NAMES.length;

// ---------------------------------------------------------------- random
export function rng(seed) {
  let a = (seed >>> 0) || 1;
  return function () {
    a = (a + 0x6D2B79F5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
export function gauss(r) {
  let u = 0, v = 0;
  while (u === 0) u = r();
  while (v === 0) v = r();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}
// Marsaglia-Tsang, shape >= 1
export function gammaSample(r, shape) {
  const d = shape - 1 / 3, c = 1 / Math.sqrt(9 * d);
  for (;;) {
    let x, v;
    do { x = gauss(r); v = 1 + c * x; } while (v <= 0);
    v = v * v * v;
    const u = r();
    if (u < 1 - 0.0331 * x * x * x * x) return d * v;
    if (Math.log(u) < 0.5 * x * x + d * (1 - v + Math.log(v))) return d * v;
  }
}

// ---------------------------------------------------------------- coordinates
export const toU = (E) => Math.log(E / EREF) / USCALE;
export const fromU = (u) => EREF * Math.exp(USCALE * u);
export const toV = (F) => F / FSCALE;
export function basis(u, v) {
  const v2 = v * v, v3 = v2 * v, u2 = u * u;
  const v4 = v2 * v2;
  return [1, u, u2, v, v2, v3, v4, u * v, u * v2, u * v3, u2 * v2, v4 * v, v4 * v2, u * v4];
}
// CD(u) at a fixed focus is A + B u + C u^2. Return [A,B,C] for coefficients th at focus v.
export function abcAt(th, v) {
  const v2 = v * v, v3 = v2 * v, v4 = v2 * v2;
  return [
    th[0] + th[3] * v + th[4] * v2 + th[5] * v3 + th[6] * v4 + th[11] * v4 * v + th[12] * v4 * v2,
    th[1] + th[7] * v + th[8] * v2 + th[9] * v3 + th[13] * v4,
    th[2] + th[10] * v2,
  ];
}
export const dot = (a, b) => { let s = 0; for (let i = 0; i < a.length; i++) s += a[i] * b[i]; return s; };

// ---------------------------------------------------------------- FEM layout
// Focus varies by column, dose by row (the usual FEM layout). Dies whose four
// corners do not fit inside the wafer minus edge exclusion are not exposed.
export function femSites() {
  const sites = [];
  const nc = FOCUS_COLS.length, nr = DOSE_ROWS.length;
  const rMax = DIE.waferR - DIE.edge;
  for (let r = 0; r < nr; r++) {
    for (let c = 0; c < nc; c++) {
      const x0 = (c - nc / 2) * DIE.w, y0 = (r - nr / 2) * DIE.h;
      const corners = [[x0, y0], [x0 + DIE.w, y0], [x0, y0 + DIE.h], [x0 + DIE.w, y0 + DIE.h]];
      const inside = corners.every(([x, y]) => Math.hypot(x, y) <= rMax);
      if (!inside) continue;
      sites.push({ id: sites.length, row: r, col: c, x: x0, y: y0, focus: FOCUS_COLS[c], dose: DOSE_ROWS[r], u: toU(DOSE_ROWS[r]), v: toV(FOCUS_COLS[c]) });
    }
  }
  return sites;
}

// ---------------------------------------------------------------- simulated tool (ground truth)
// A synthetic exposure tool + resist. CD of feature k:
//   CD = T + S l + q l^2 + (b + g l) w^2 + c w^3 + d w^4 - h max(0,|w|-0.9)^2
// with l = ln(E/Ek) and w = (F - F0k - Ftool)/100. The last (pattern top-loss)
// term is deliberately outside the fitted polynomial basis.
const NOMINAL = {
  dense: { Ek: 30.0, S: -38, q: 30, b: -9.0, g: 45, F0: 0, c: 2.0, d: -2.0, h: 40 },
  semi: { Ek: 31.4, S: -42, q: 30, b: -11.0, g: 60, F0: 6, c: 3.0, d: -2.5, h: 40 },
  iso: { Ek: 32.9, S: -47, q: 35, b: -16.0, g: 90, F0: -14, c: 5.0, d: -3.0, h: 50 },
};
export function makeTool(seed) {
  const r = rng(seed * 7919 + 13);
  const j = (s) => 1 + s * gauss(r);
  const Ftool = 12 * gauss(r);
  const sigma = 0.55 + 0.35 * r();
  const feats = FEATURES.map((f) => {
    const n = NOMINAL[f.id];
    return {
      Ek: n.Ek * j(0.015), S: n.S * j(0.08), q: n.q * j(0.2), b: n.b * j(0.12), g: n.g * j(0.15),
      F0: n.F0 + 4 * gauss(r), c: n.c * j(0.4) * (r() < 0.25 ? -1 : 1), d: n.d * j(0.3), h: n.h * j(0.2),
    };
  });
  return { seed, Ftool, sigma, feats };
}
export function trueCD(tool, k, E, F) {
  const p = tool.feats[k];
  const l = Math.log(E / p.Ek);
  const w = (F - p.F0 - tool.Ftool) / 100;
  const over = Math.max(0, Math.abs(w) - 0.9);
  return SPEC.target + p.S * l + p.q * l * l + (p.b + p.g * l) * w * w + p.c * w * w * w + p.d * w * w * w * w - p.h * over * over;
}
// Same model rewritten as A + B u + C u^2 at focus F (exact, used for the true window).
export function trueABC(tool, k, F) {
  const p = tool.feats[k];
  const l0 = Math.log(EREF / p.Ek);              // l = l0 + USCALE*u
  const w = (F - p.F0 - tool.Ftool) / 100;
  const over = Math.max(0, Math.abs(w) - 0.9);
  const A = SPEC.target + p.S * l0 + p.q * l0 * l0 + (p.b + p.g * l0) * w * w + p.c * w ** 3 + p.d * w ** 4 - p.h * over * over;
  const B = USCALE * (p.S + 2 * p.q * l0 + p.g * w * w);
  const C = USCALE * USCALE * p.q;
  return [A, B, C];
}
export function measure(tool, site, r) {
  return FEATURES.map((_, k) => trueCD(tool, k, site.dose, site.focus) + tool.sigma * gauss(r));
}

// ---------------------------------------------------------------- linear algebra (small, dense)
export function cholesky(M) {
  const n = M.length, L = Array.from({ length: n }, () => new Float64Array(n));
  for (let i = 0; i < n; i++) {
    for (let j = 0; j <= i; j++) {
      let s = M[i][j];
      for (let k = 0; k < j; k++) s -= L[i][k] * L[j][k];
      if (i === j) {
        if (s <= 0) s = 1e-12;
        L[i][j] = Math.sqrt(s);
      } else L[i][j] = s / L[j][j];
    }
  }
  return L;
}
export function invSPD(M) {
  const n = M.length, L = cholesky(M);
  const Li = Array.from({ length: n }, () => new Float64Array(n));
  for (let i = 0; i < n; i++) {
    Li[i][i] = 1 / L[i][i];
    for (let j = 0; j < i; j++) {
      let s = 0;
      for (let k = j; k < i; k++) s -= L[i][k] * Li[k][j];
      Li[i][j] = s / L[i][i];
    }
  }
  const out = Array.from({ length: n }, () => new Float64Array(n));
  for (let i = 0; i < n; i++) for (let j = 0; j <= i; j++) {
    let s = 0;
    for (let k = i; k < n; k++) s += Li[k][i] * Li[k][j];
    out[i][j] = s; out[j][i] = s;
  }
  return out;
}
const matVec = (M, x) => M.map((row) => dot(row, x));

// ---------------------------------------------------------------- Bayesian Bossung fit
// Normal-inverse-gamma conjugate regression per feature.
export const PRIOR = {
  mean: [SPEC.target, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
  sd: [10, 10, 5, 10, 10, 10, 10, 10, 10, 10, 5, 5, 5, 5],   // nm, in units of sigma=1 nm
  a0: 3, b0: 1.0,                                     // prior mean sigma^2 = 0.5 nm^2
};
export function fitFeature(rows, ys, prior = PRIOR) {
  const p = P;
  const V0inv = Array.from({ length: p }, (_, i) => { const r = new Float64Array(p); r[i] = 1 / (prior.sd[i] ** 2); return r; });
  const prec = V0inv.map((r) => Float64Array.from(r));
  const rhs = prior.mean.map((m, i) => m / (prior.sd[i] ** 2));
  let yy = 0;
  for (let n = 0; n < rows.length; n++) {
    const phi = rows[n], y = ys[n];
    for (let i = 0; i < p; i++) {
      rhs[i] += phi[i] * y;
      for (let j = 0; j < p; j++) prec[i][j] += phi[i] * phi[j];
    }
    yy += y * y;
  }
  const Vn = invSPD(prec);
  const mn = matVec(Vn, rhs);
  const m0Vm0 = prior.mean.reduce((s, m, i) => s + m * m / (prior.sd[i] ** 2), 0);
  const an = prior.a0 + rows.length / 2;
  const bn = prior.b0 + 0.5 * (yy + m0Vm0 - dot(mn, matVec(prec, mn)));
  return { mn, Vn, an, bn: Math.max(bn, 1e-6), n: rows.length };
}
// Sherman-Morrison: covariance after one more observation at phi (y-independent).
export function addPoint(Vn, phi) {
  const Vp = matVec(Vn, phi);
  const den = 1 + dot(phi, Vp);
  return Vn.map((row, i) => row.map((x, j) => x - Vp[i] * Vp[j] / den));
}
// Common random numbers so posteriors before/after a hypothetical site are comparable.
export function makeCRN(S, seed) {
  const r = rng(seed);
  return FEATURES.map(() => ({ z: Array.from({ length: S }, () => Array.from({ length: P }, () => gauss(r))), u: Array.from({ length: S }, () => r()) }));
}
export function drawThetas(fit, crn, Vn = fit.Vn) {
  const L = cholesky(Vn);
  const S = crn.z.length;
  // sigma^2 ~ InvGamma(an, bn) via a stratified gamma quantile approximation (Wilson-Hilferty)
  const out = new Array(S);
  for (let s = 0; s < S; s++) {
    const g = gammaQuantileWH(fit.an, crn.u[s]);
    const sig = Math.sqrt(fit.bn / g);
    const z = crn.z[s];
    const th = new Float64Array(P);
    for (let i = 0; i < P; i++) {
      let acc = 0;
      for (let k = 0; k <= i; k++) acc += L[i][k] * z[k];
      th[i] = fit.mn[i] + sig * acc;
    }
    out[s] = { th, sig };
  }
  return out;
}
// Gamma(shape,1) quantile via Wilson-Hilferty (accurate for shape >= 3)
export function gammaQuantileWH(shape, p) {
  const z = normInv(Math.min(Math.max(p, 1e-9), 1 - 1e-9));
  const t = 1 - 1 / (9 * shape) + z / (3 * Math.sqrt(shape));
  return shape * Math.max(t, 1e-6) ** 3;
}
export function normInv(p) {
  // Acklam's rational approximation
  const a = [-39.69683028665376, 220.9460984245205, -275.9285104469687, 138.357751867269, -30.66479806614716, 2.506628277459239];
  const b = [-54.47609879822406, 161.5858368580409, -155.6989798598866, 66.80131188771972, -13.28068155211781];
  const c = [-0.007784894002430293, -0.3223964580411365, -2.400758277161838, -2.549732539343734, 4.374664141464968, 2.938163982698783];
  const d = [0.007784695709041462, 0.3224671290700398, 2.445134137142996, 3.754408661907416];
  const pl = 0.02425;
  if (p < pl) { const q = Math.sqrt(-2 * Math.log(p)); return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1); }
  if (p > 1 - pl) { const q = Math.sqrt(-2 * Math.log(1 - p)); return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1); }
  const q = p - 0.5, r = q * q;
  return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
}

// ---------------------------------------------------------------- process window geometry
export const FOCUS_GRID = Array.from({ length: 131 }, (_, i) => -130 + 2 * i);
export const U_MIN = toU(DOSE_ROWS[0]), U_MAX = toU(DOSE_ROWS[DOSE_ROWS.length - 1]);

// Longest sub-interval of [lo,hi] where L <= A + B u + C u^2 <= H.
export function specInterval(A, B, C, L, H, lo = U_MIN, hi = U_MAX) {
  const pts = [lo, hi];
  const roots = (k) => {
    const a = C, b = B, c = A - k;
    if (Math.abs(a) < 1e-12) { if (Math.abs(b) > 1e-12) pts.push(-c / b); return; }
    const disc = b * b - 4 * a * c;
    if (disc < 0) return;
    const s = Math.sqrt(disc);
    pts.push((-b - s) / (2 * a), (-b + s) / (2 * a));
  };
  roots(L); roots(H);
  const xs = pts.filter((x) => x >= lo && x <= hi).sort((x, y) => x - y);
  let best = null;
  let cur = null;
  for (let i = 0; i < xs.length - 1; i++) {
    const m = 0.5 * (xs[i] + xs[i + 1]);
    const y = A + B * m + C * m * m;
    const ok = y >= L && y <= H;
    if (ok) {
      if (cur && Math.abs(cur[1] - xs[i]) < 1e-12) cur[1] = xs[i + 1];
      else cur = [xs[i], xs[i + 1]];
      if (!best || cur[1] - cur[0] > best[1] - best[0]) best = [cur[0], cur[1]];
    } else cur = null;
  }
  return best; // null if nowhere in spec
}
export const elOf = (ulo, uhi) => { const e1 = fromU(ulo), e2 = fromU(uhi); return (e2 - e1) / ((e1 + e2) / 2); };

// abcFns: array over features of (focusIndex) -> [A,B,C]
export function windowFromABC(abcGrid, elReq) {
  const nf = FOCUS_GRID.length;
  const LO = new Float64Array(nf), HI = new Float64Array(nf);
  const L = SPEC.target * (1 - SPEC.tol), H = SPEC.target * (1 + SPEC.tol);
  for (let i = 0; i < nf; i++) {
    let lo = U_MIN, hi = U_MAX;
    for (let k = 0; k < abcGrid.length; k++) {
      const [A, B, C] = abcGrid[k][i];
      const iv = specInterval(A, B, C, L, H);
      if (!iv) { lo = 1; hi = 0; break; }
      lo = Math.max(lo, iv[0]); hi = Math.min(hi, iv[1]);
    }
    LO[i] = lo; HI[i] = hi;
  }
  let best = { dof: 0, f1: null, f2: null, ulo: null, uhi: null, el: 0 };
  for (let i = 0; i < nf; i++) {
    let lo = -Infinity, hi = Infinity;
    for (let j = i; j < nf; j++) {
      lo = Math.max(lo, LO[j]); hi = Math.min(hi, HI[j]);
      if (hi <= lo) break;
      const el = elOf(lo, hi);
      if (el < elReq) break;
      const dof = FOCUS_GRID[j] - FOCUS_GRID[i];
      if (dof > best.dof || (best.f1 === null)) best = { dof, f1: FOCUS_GRID[i], f2: FOCUS_GRID[j], ulo: lo, uhi: hi, el };
    }
  }
  return { ...best, LO, HI };
}
// ED curve: DOF as a function of required EL.
export const EL_LEVELS = Array.from({ length: 31 }, (_, i) => i * 0.005); // 0..15%
export function edCurve(abcGrid) {
  return EL_LEVELS.map((el) => windowFromABC(abcGrid, el).dof);
}
export function abcGridFromThetas(thetas) {
  return thetas.map((th) => FOCUS_GRID.map((F) => abcAt(th, toV(F))));
}
export function trueWindow(tool, elReq) {
  const grid = FEATURES.map((_, k) => FOCUS_GRID.map((F) => trueABC(tool, k, F)));
  return windowFromABC(grid, elReq);
}
export function trueED(tool) {
  const grid = FEATURES.map((_, k) => FOCUS_GRID.map((F) => trueABC(tool, k, F)));
  return edCurve(grid);
}

// ---------------------------------------------------------------- posterior over the window
export function quantile(arr, q) {
  const a = Float64Array.from(arr).sort();
  const pos = (a.length - 1) * q, lo = Math.floor(pos), hi = Math.ceil(pos);
  return a[lo] + (a[hi] - a[lo]) * (pos - lo);
}
export function mean(a) { let s = 0; for (const x of a) s += x; return s / a.length; }

export function fitAll(sites, obs) {
  const rows = sites.map((s) => basis(s.u, s.v));
  return FEATURES.map((_, k) => fitFeature(rows, obs.map((o) => o[k])));
}
// samples[s] = per-feature theta draws, plus window metrics
export function posteriorWindow(fits, crns, elReq, Vs = null) {
  const draws = fits.map((f, k) => drawThetas(f, crns[k], Vs ? Vs[k] : f.Vn));
  const S = draws[0].length;
  const dof = new Float64Array(S), bestF = new Float64Array(S), bestE = new Float64Array(S);
  for (let s = 0; s < S; s++) {
    const grid = abcGridFromThetas(draws.map((d) => d[s].th));
    const w = windowFromABC(grid, elReq);
    dof[s] = w.dof;
    bestF[s] = w.f1 === null ? NaN : (w.f1 + w.f2) / 2;
    bestE[s] = w.ulo === null ? NaN : fromU((w.ulo + w.uhi) / 2);
  }
  return { draws, dof, bestF, bestE };
}
export function summarize(post, level = 0.9) {
  const a = (1 - level) / 2;
  const ok = (arr) => Array.from(arr).filter((x) => Number.isFinite(x));
  const bf = ok(post.bestF), be = ok(post.bestE);
  return {
    dofMed: quantile(post.dof, 0.5), dofLo: quantile(post.dof, a), dofHi: quantile(post.dof, 1 - a),
    bestFocus: bf.length ? quantile(bf, 0.5) : NaN, bestDose: be.length ? quantile(be, 0.5) : NaN,
    bestFocusLo: bf.length ? quantile(bf, a) : NaN, bestFocusHi: bf.length ? quantile(bf, 1 - a) : NaN,
  };
}

// ---------------------------------------------------------------- site selection
// Score = information a site's 3 CD readings carry about DOF, in a linear-Gaussian
// (ensemble-sensitivity) approximation: sum_k Cov(DOF, f_k(x))^2 / (Var f_k(x) + sigma_k^2).
export function scoreSites(post, fits, candidates) {
  const S = post.dof.length;
  const dm = mean(post.dof);
  let dv = 0; for (let s = 0; s < S; s++) dv += (post.dof[s] - dm) ** 2; dv /= S;
  return candidates.map((site) => {
    const phi = basis(site.u, site.v);
    let score = 0, predVar = 0;
    for (let k = 0; k < fits.length; k++) {
      const f = new Float64Array(S);
      let fm = 0, sig2 = 0;
      for (let s = 0; s < S; s++) { f[s] = dot(post.draws[k][s].th, phi); fm += f[s]; sig2 += post.draws[k][s].sig ** 2; }
      fm /= S; sig2 /= S;
      let cov = 0, vf = 0;
      for (let s = 0; s < S; s++) { cov += (post.dof[s] - dm) * (f[s] - fm); vf += (f[s] - fm) ** 2; }
      cov /= S; vf /= S;
      score += cov * cov / (vf + sig2);
      predVar += vf;
    }
    return { site, score: dv > 1e-9 ? score / dv : 0, predVar };
  });
}
// Greedy batch: after each pick, shrink each feature's covariance with the
// y-independent rank-one update and re-draw with the same random numbers.
export function selectBatch(fits, crns, elReq, candidates, k) {
  let Vs = fits.map((f) => f.Vn);
  const chosen = [];
  let pool = candidates.slice();
  const trace = [];
  for (let t = 0; t < k && pool.length; t++) {
    const post = posteriorWindow(fits, crns, elReq, Vs);
    const scored = scoreSites(post, fits, pool);
    // fallback when DOF posterior is degenerate: largest predictive variance
    const useVar = scored.every((s) => s.score === 0);
    scored.sort((a, b) => (useVar ? b.predVar - a.predVar : b.score - a.score));
    const pick = scored[0];
    chosen.push(pick.site);
    trace.push({ id: pick.site.id, score: pick.score, sd: Math.sqrt(varOf(post.dof)) });
    const phi = basis(pick.site.u, pick.site.v);
    Vs = Vs.map((V) => addPoint(V, phi));
    pool = pool.filter((s) => s.id !== pick.site.id);
  }
  return { chosen, trace };
}
export function varOf(a) { const m = mean(a); let v = 0; for (const x of a) v += (x - m) ** 2; return v / a.length; }

// Maximin space-filling design over normalised (focus, dose) coordinates.
export function spaceFilling(sites, k, start = []) {
  const chosen = start.slice();
  const norm = (s) => [s.col / (FOCUS_COLS.length - 1), s.row / (DOSE_ROWS.length - 1)];
  if (!chosen.length) {
    // start from the site nearest the matrix centre
    let best = sites[0], bd = Infinity;
    for (const s of sites) { const [a, b] = norm(s); const d = (a - 0.5) ** 2 + (b - 0.5) ** 2; if (d < bd) { bd = d; best = s; } }
    chosen.push(best);
  }
  while (chosen.length < k) {
    let best = null, bd = -1;
    for (const s of sites) {
      if (chosen.some((c) => c.id === s.id)) continue;
      const [a, b] = norm(s);
      let dmin = Infinity;
      for (const c of chosen) { const [x, y] = norm(c); dmin = Math.min(dmin, (a - x) ** 2 + (b - y) ** 2); }
      if (dmin > bd + 1e-12) { bd = dmin; best = s; }
    }
    if (!best) break;
    chosen.push(best);
  }
  return chosen;
}

// ---------------------------------------------------------------- full campaign
export const DEFAULTS = { seedSites: 16, rounds: 3, batch: 8, samples: 300, elReq: 0.05 };
export function runCampaign(tool, strategy = 'adaptive', opts = {}) {
  const o = { ...DEFAULTS, ...opts };
  const sites = femSites();
  const r = rng(tool.seed * 104729 + 7);
  // each site has one fixed measurement (the same wafer, measured once)
  const obsAll = sites.map((s) => measure(tool, s, r));
  const crns = makeCRN(o.samples, tool.seed + 99);
  const budget = o.seedSites + o.rounds * o.batch;
  let order;
  const history = [];
  const record = (ids, label) => {
    const ss = ids.map((i) => sites[i]);
    const fits = fitAll(ss, ids.map((i) => obsAll[i]));
    const post = posteriorWindow(fits, crns, o.elReq);
    const sum = summarize(post);
    history.push({ label, ids: ids.slice(), fits, post, sum });
    return { fits, post };
  };
  if (strategy === 'full') {
    record(sites.map((s) => s.id), 'Full FEM');
  } else if (strategy === 'grid') {
    order = spaceFilling(sites, budget).map((s) => s.id);
    record(order, 'Space-filling');
  } else if (strategy === 'random') {
    const seedIds = spaceFilling(sites, o.seedSites).map((s) => s.id);
    const rest = sites.filter((s) => !seedIds.includes(s.id)).map((s) => s.id);
    const rr = rng(tool.seed * 31 + 5);
    for (let i = rest.length - 1; i > 0; i--) { const j = Math.floor(rr() * (i + 1)); [rest[i], rest[j]] = [rest[j], rest[i]]; }
    record(seedIds.concat(rest.slice(0, budget - o.seedSites)), 'Random');
  } else {
    let ids = spaceFilling(sites, o.seedSites).map((s) => s.id);
    let cur = record(ids, 'Seed');
    for (let rd = 1; rd <= o.rounds; rd++) {
      const pool = sites.filter((s) => !ids.includes(s.id));
      const { chosen, trace } = selectBatch(cur.fits, crns, o.elReq, pool, o.batch);
      ids = ids.concat(chosen.map((s) => s.id));
      cur = record(ids, 'Round ' + rd);
      history[history.length - 1].picked = chosen.map((s) => s.id);
      history[history.length - 1].trace = trace;
    }
  }
  return { sites, obsAll, history, truth: trueWindow(tool, o.elReq), opts: o };
}

// ---------------------------------------------------------------- page helpers
// ED curve (DOF vs required EL) for one set of per-feature coefficient vectors;
// computes the per-focus dose intervals once and reuses them for every EL level.
export function edFromThetas(thetas, levels = EL_LEVELS) {
  const grid = abcGridFromThetas(thetas);
  const base = windowFromABC(grid, 1e9); // only to get LO/HI
  const nf = FOCUS_GRID.length, { LO, HI } = base;
  return levels.map((elReq) => {
    let dof = 0;
    for (let i = 0; i < nf; i++) {
      let lo = -Infinity, hi = Infinity;
      for (let j = i; j < nf; j++) {
        lo = Math.max(lo, LO[j]); hi = Math.min(hi, HI[j]);
        if (hi <= lo || elOf(lo, hi) < elReq) break;
        dof = Math.max(dof, FOCUS_GRID[j] - FOCUS_GRID[i]);
      }
    }
    return dof;
  });
}
// Posterior probability that every feature is in spec at each (focus, dose) cell.
export function inSpecMap(draws, focusAxis, doseAxis) {
  const L = SPEC.target * (1 - SPEC.tol), H = SPEC.target * (1 + SPEC.tol);
  const S = draws[0].length;
  const out = Array.from({ length: doseAxis.length }, () => new Float64Array(focusAxis.length));
  const us = doseAxis.map(toU);
  for (let s = 0; s < S; s++) {
    const abc = draws.map((d) => focusAxis.map((F) => abcAt(d[s].th, toV(F))));
    for (let i = 0; i < focusAxis.length; i++) {
      for (let j = 0; j < us.length; j++) {
        const u = us[j];
        let ok = true;
        for (let k = 0; k < abc.length; k++) {
          const [A, B, C] = abc[k][i];
          const y = A + B * u + C * u * u;
          if (y < L || y > H) { ok = false; break; }
        }
        if (ok) out[j][i] += 1 / S;
      }
    }
  }
  return out;
}
// Posterior mean and 90% band of CD for feature k along focus at dose E.
export function bossungBand(draws, k, E, focusAxis) {
  const u = toU(E);
  return focusAxis.map((F) => {
    const v = toV(F), phi = basis(u, v);
    const ys = draws[k].map((d) => dot(d.th, phi));
    return { F, mid: quantile(ys, 0.5), lo: quantile(ys, 0.05), hi: quantile(ys, 0.95) };
  });
}
