import test from 'node:test';
import assert from 'node:assert/strict';
import * as pw from '../core/pw.js';

test('FEM layout keeps only dies fully inside the wafer edge exclusion', () => {
  const sites = pw.femSites();
  assert.equal(sites.length, 111);
  const rMax = pw.DIE.waferR - pw.DIE.edge;
  for (const s of sites) {
    for (const [dx, dy] of [[0, 0], [pw.DIE.w, 0], [0, pw.DIE.h], [pw.DIE.w, pw.DIE.h]]) {
      assert.ok(Math.hypot(s.x + dx, s.y + dy) <= rMax + 1e-9);
    }
  }
  // the centre die (best guess dose/focus) is always exposed
  assert.ok(sites.some((s) => s.focus === 0 && s.dose === 31.5));
});

test('trueABC is the same function as trueCD rewritten in (u)', () => {
  const tool = pw.makeTool(3);
  for (let k = 0; k < 3; k++) for (const F of [-110, -37, 0, 52, 118]) for (const E of [26.3, 30, 33.3, 36.8]) {
    const [A, B, C] = pw.trueABC(tool, k, F);
    const u = pw.toU(E);
    assert.ok(Math.abs(A + B * u + C * u * u - pw.trueCD(tool, k, E, F)) < 1e-9);
  }
});

test('abcAt agrees with the basis expansion', () => {
  const r = pw.rng(5);
  const th = Array.from({ length: pw.P }, () => pw.gauss(r));
  for (const v of [-1.2, -0.3, 0, 0.7]) for (const u of [-0.8, 0.1, 0.9]) {
    const [A, B, C] = pw.abcAt(th, v);
    assert.ok(Math.abs(A + B * u + C * u * u - pw.dot(th, pw.basis(u, v))) < 1e-9);
  }
});

test('specInterval matches a brute-force scan', () => {
  const r = pw.rng(11);
  for (let t = 0; t < 300; t++) {
    const A = 40 + 6 * pw.gauss(r), B = -8 + 4 * pw.gauss(r), C = 3 * pw.gauss(r);
    const iv = pw.specInterval(A, B, C, 36, 44, -1, 1);
    // brute force longest run
    let best = 0, run = 0, step = 1e-4;
    for (let u = -1; u <= 1; u += step) {
      const y = A + B * u + C * u * u;
      if (y >= 36 && y <= 44) { run += step; best = Math.max(best, run); } else run = 0;
    }
    const len = iv ? iv[1] - iv[0] : 0;
    assert.ok(Math.abs(len - best) < 3e-3, `len ${len} vs ${best}`);
  }
});

test('window: a focus-independent process gives the full focus range; ED curve is non-increasing', () => {
  // CD = 40 - 8u at every focus: in spec for |u| <= 0.5
  const grid = [pw.FOCUS_GRID.map(() => [40, -8, 0])];
  const w = pw.windowFromABC(grid, 0.05);
  assert.equal(w.dof, pw.FOCUS_GRID.at(-1) - pw.FOCUS_GRID[0]);
  assert.ok(Math.abs(w.ulo + 0.5) < 1e-9 && Math.abs(w.uhi - 0.5) < 1e-9);
  const ed = pw.trueED(pw.makeTool(4));
  for (let i = 1; i < ed.length; i++) assert.ok(ed[i] <= ed[i - 1]);
});

test('window: parabolic Bossung matches the closed-form window', () => {
  // CD = 40 - 8u - 10 v^2 (v = F/100), spec 36..44  =>  in spec for
  // -0.5 - 1.25 v^2 <= u <= 0.5 - 1.25 v^2 at each focus.
  const grid = [pw.FOCUS_GRID.map((F) => [40 - 10 * (F / 100) ** 2, -8, 0])];
  for (const elReq of [0.03, 0.05, 0.08]) {
    const w = pw.windowFromABC(grid, elReq);
    let bestDof = 0;
    for (let i = 0; i < pw.FOCUS_GRID.length; i++) for (let j = i; j < pw.FOCUS_GRID.length; j++) {
      const a = pw.FOCUS_GRID[i] / 100, b = pw.FOCUS_GRID[j] / 100;
      const vmin2 = a <= 0 && b >= 0 ? 0 : Math.min(a * a, b * b);
      const vmax2 = Math.max(a * a, b * b);
      const lo = Math.max(pw.U_MIN, -0.5 - 1.25 * vmin2);
      const hi = Math.min(pw.U_MAX, 0.5 - 1.25 * vmax2);
      if (hi > lo && pw.elOf(lo, hi) >= elReq) bestDof = Math.max(bestDof, pw.FOCUS_GRID[j] - pw.FOCUS_GRID[i]);
    }
    assert.equal(w.dof, bestDof, `EL ${elReq}`);
    assert.ok(w.dof > 0);
  }
});

test('Bayesian fit recovers coefficients from noiseless data', () => {
  const r = pw.rng(2);
  const th = Array.from({ length: pw.P }, (_, i) => (i === 0 ? 40 : pw.gauss(r)));
  const sites = pw.femSites();
  const rows = sites.map((s) => pw.basis(s.u, s.v));
  const ys = rows.map((phi) => pw.dot(th, phi));
  const fit = pw.fitFeature(rows, ys, { ...pw.PRIOR, sd: pw.PRIOR.sd.map(() => 1e4) });
  for (let i = 0; i < pw.P; i++) assert.ok(Math.abs(fit.mn[i] - th[i]) < 1e-5, `coef ${i}`);
});

test('Sherman-Morrison update equals refitting with the extra row', () => {
  const sites = pw.femSites().slice(0, 30);
  const rows = sites.map((s) => pw.basis(s.u, s.v));
  const ys = rows.map(() => 40);
  const fit = pw.fitFeature(rows, ys);
  const extra = pw.femSites()[60];
  const phi = pw.basis(extra.u, extra.v);
  const V1 = pw.addPoint(fit.Vn, phi);
  const fit2 = pw.fitFeature(rows.concat([phi]), ys.concat([40]));
  for (let i = 0; i < pw.P; i++) for (let j = 0; j < pw.P; j++) assert.ok(Math.abs(V1[i][j] - fit2.Vn[i][j]) < 1e-9);
});

test('Wilson-Hilferty gamma quantiles have the right mean', () => {
  const shape = 30;
  let s = 0; const n = 4000;
  for (let i = 0; i < n; i++) s += pw.gammaQuantileWH(shape, (i + 0.5) / n);
  assert.ok(Math.abs(s / n - shape) / shape < 0.01);
});

test('adaptive selection beats random at the same budget and stays calibrated (12 tools)', () => {
  let errA = 0, errR = 0, covA = 0;
  const N = 12;
  for (let seed = 1; seed <= N; seed++) {
    const tool = pw.makeTool(seed);
    const truth = pw.trueWindow(tool, 0.05).dof;
    const a = pw.runCampaign(tool, 'adaptive').history.at(-1).sum;
    const b = pw.runCampaign(tool, 'random').history.at(-1).sum;
    errA += Math.abs(a.dofMed - truth); errR += Math.abs(b.dofMed - truth);
    if (truth >= a.dofLo && truth <= a.dofHi) covA++;
  }
  assert.ok(errA < errR, `adaptive ${errA / N} vs random ${errR / N}`);
  assert.ok(covA / N >= 0.75, `coverage ${covA / N}`);
});

test('campaigns are deterministic for a given tool seed', () => {
  const a = pw.runCampaign(pw.makeTool(9), 'adaptive');
  const b = pw.runCampaign(pw.makeTool(9), 'adaptive');
  assert.deepEqual(a.history.map((h) => h.ids), b.history.map((h) => h.ids));
  assert.equal(a.history.at(-1).sum.dofMed, b.history.at(-1).sum.dofMed);
});
