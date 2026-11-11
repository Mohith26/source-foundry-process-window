// Backtest: run every site-selection strategy on N simulated tools and record
// how close each one gets to the true depth of focus, and whether the 90%
// interval covers it. Writes data/backtest.json (read by the page).
// Usage: node scripts/backtest.mjs [N]
import * as pw from '../core/pw.js';
import { writeFileSync } from 'node:fs';

const N = +process.argv[2] || 100;
const EL = 0.05;
const budgets = [16, 24, 32, 40];
const strategies = ['adaptive', 'random', 'grid'];
const rows = [];
const t0 = Date.now();
for (let seed = 1; seed <= N; seed++) {
  const tool = pw.makeTool(seed);
  const truth = pw.trueWindow(tool, EL).dof;
  const row = { seed, truth, est: {} };
  // adaptive: one campaign, read every round
  const ad = pw.runCampaign(tool, 'adaptive', { elReq: EL });
  ad.history.forEach((h, i) => { row.est['adaptive@' + budgets[i]] = [h.sum.dofMed, h.sum.dofLo, h.sum.dofHi]; });
  for (const st of ['random', 'grid']) {
    for (let r = 1; r <= 3; r++) {
      const c = pw.runCampaign(tool, st, { elReq: EL, rounds: r });
      const h = c.history.at(-1).sum;
      row.est[st + '@' + (16 + 8 * r)] = [h.dofMed, h.dofLo, h.dofHi];
    }
  }
  const full = pw.runCampaign(tool, 'full', { elReq: EL });
  const hf = full.history.at(-1).sum;
  row.est['full@' + full.sites.length] = [hf.dofMed, hf.dofLo, hf.dofHi];
  row.nSites = full.sites.length;
  rows.push(row);
  if (seed % 10 === 0) console.error(`${seed}/${N} tools, ${((Date.now() - t0) / 1000).toFixed(0)} s`);
}
const keys = Object.keys(rows[0].est);
const summary = {};
for (const k of keys) {
  const err = rows.map((r) => Math.abs(r.est[k][0] - r.truth));
  const cov = rows.filter((r) => r.truth >= r.est[k][1] - 1e-9 && r.truth <= r.est[k][2] + 1e-9).length / rows.length;
  const width = rows.map((r) => r.est[k][2] - r.est[k][1]);
  summary[k] = {
    strategy: k.split('@')[0], sites: +k.split('@')[1],
    meanErr: pw.mean(err), medErr: pw.quantile(err, 0.5), p90Err: pw.quantile(err, 0.9),
    coverage90: cov, medWidth: pw.quantile(width, 0.5),
  };
}
const out = { generated: 'scripts/backtest.mjs', N, elReq: EL, nSites: rows[0].nSites, budgets, summary, rows };
writeFileSync(new URL('../data/backtest.json', import.meta.url), JSON.stringify(out));
console.log(JSON.stringify(summary, null, 1));
console.error(`done in ${((Date.now() - t0) / 1000).toFixed(0)} s`);
