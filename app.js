import * as pw from './core/pw.js';

const TOOL_SEED = 18;          // a typical tool from the 100-tool backtest (adaptive error = median)
const EL_REQ = 0.05;
const STOP_HALF_WIDTH = 5;     // nm, stop rule on the 90% interval
const NS = 'http://www.w3.org/2000/svg';
const tip = document.getElementById('tip');

// ------------------------------------------------------------- compute
const tool = pw.makeTool(TOOL_SEED);
const camp = pw.runCampaign(tool, 'adaptive', { elReq: EL_REQ });
const full = pw.runCampaign(tool, 'full', { elReq: EL_REQ });
const sites = camp.sites;
const truth = camp.truth;
const truthED = pw.trueED(tool);
const stages = camp.history.map((h, i) => ({ ...h, idx: i }));
stages.push({ ...full.history[0], idx: stages.length, label: 'Full wafer', isFull: true });
// next batch after the last adaptive round (only needed if the stop rule is not met)
const lastAd = stages[stages.length - 2];
const meets = (s) => (s.sum.dofHi - s.sum.dofLo) / 2 <= STOP_HALF_WIDTH;

const FOCUS_AX = Array.from({ length: 62 }, (_, i) => -122 + 4 * i);
const DOSE_AX = Array.from({ length: 56 }, (_, i) => 26 + i * 0.2);
const cache = new Map();
function derived(st) {
  if (cache.has(st.idx)) return cache.get(st.idx);
  const post = st.post;
  // representative draw: DOF equal to the median, best focus nearest the median
  const S = post.dof.length;
  const medDof = st.sum.dofMed;
  let rep = 0, bd = Infinity;
  for (let s = 0; s < S; s++) {
    const d = Math.abs(post.dof[s] - medDof) * 10 + Math.abs((post.bestF[s] || 0) - st.sum.bestFocus) + Math.abs((post.bestE[s] || 0) - st.sum.bestDose) * 4;
    if (d < bd) { bd = d; rep = s; }
  }
  const repWin = pw.windowFromABC(pw.abcGridFromThetas(post.draws.map((d) => d[rep].th)), EL_REQ);
  // window-edge spread
  const f1 = [], f2 = [];
  for (let s = 0; s < S; s++) {
    const w = pw.windowFromABC(pw.abcGridFromThetas(post.draws.map((d) => d[s].th)), EL_REQ);
    if (w.f1 !== null && w.dof > 0) { f1.push(w.f1); f2.push(w.f2); }
  }
  const map = pw.inSpecMap(post.draws, FOCUS_AX, DOSE_AX);
  const nED = 120;
  const eds = [];
  for (let s = 0; s < nED; s++) eds.push(pw.edFromThetas(post.draws.map((d) => d[s].th)));
  const edBand = pw.EL_LEVELS.map((_, i) => {
    const col = eds.map((e) => e[i]);
    return { lo: pw.quantile(col, 0.05), mid: pw.quantile(col, 0.5), hi: pw.quantile(col, 0.95) };
  });
  // the batch that this stage queues next
  let next = null;
  const nextStage = stages[st.idx + 1];
  if (!st.isFull && nextStage && !nextStage.isFull) next = { ids: nextStage.picked, trace: nextStage.trace };
  else if (!st.isFull && !meets(st)) {
    const pool = sites.filter((x) => !st.ids.includes(x.id));
    const b = pw.selectBatch(st.fits, pw.makeCRN(post.dof.length, TOOL_SEED + 99), EL_REQ, pool, 8);
    next = { ids: b.chosen.map((x) => x.id), trace: b.trace };
  }
  const out = { repWin, edges: f1.length ? { f1lo: pw.quantile(f1, 0.05), f1hi: pw.quantile(f1, 0.95), f2lo: pw.quantile(f2, 0.05), f2hi: pw.quantile(f2, 0.95) } : null, map, edBand, next };
  cache.set(st.idx, out);
  return out;
}

// ------------------------------------------------------------- state
const qs = new URLSearchParams(location.search);
let cur = qs.has('stage') ? Math.max(0, Math.min(stages.length - 1, +qs.get('stage') || 0)) : stages.length - 2; // default: last adaptive round
let showTruth = qs.get('truth') === '1';
if (showTruth) document.getElementById('btnTruth').classList.add('on');
let feat = 2;                   // isolated line limits the window most often

// ------------------------------------------------------------- svg helpers
function el(tag, attrs = {}, parent) {
  const e = document.createElementNS(NS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]);
  if (parent) parent.appendChild(e);
  return e;
}
function txt(parent, x, y, s, cls = 'tick', anchor = 'middle', extra = {}) {
  const t = el('text', { x, y, class: cls, 'text-anchor': anchor, ...extra }, parent);
  t.textContent = s;
  return t;
}
function frame(container, W, H, m) {
  container.innerHTML = '';
  const svg = el('svg', { viewBox: `0 0 ${W} ${H}` }, container);
  return { svg, W, H, m, iw: W - m.l - m.r, ih: H - m.t - m.b };
}
function scale(d0, d1, r0, r1) { const f = (x) => r0 + (x - d0) / (d1 - d0) * (r1 - r0); f.inv = (y) => d0 + (y - r0) / (r1 - r0) * (d1 - d0); return f; }
function axes(g, f, x, y, xt, yt, xl, yl, fmtX = (v) => v, fmtY = (v) => v) {
  const { m, iw, ih } = f;
  for (const v of yt) { el('line', { x1: m.l, x2: m.l + iw, y1: y(v), y2: y(v), class: 'grid-l' }, g); txt(g, m.l - 8, y(v) + 4, fmtY(v), 'tick', 'end'); }
  for (const v of xt) { txt(g, x(v), m.t + ih + 17, fmtX(v)); }
  el('line', { x1: m.l, x2: m.l + iw, y1: m.t + ih, y2: m.t + ih, class: 'ax' }, g);
  if (xl) txt(g, m.l + iw / 2, m.t + ih + 36, xl, 'axl');
  if (yl) txt(g, 14, m.t + ih / 2, yl, 'axl', 'middle', { transform: `rotate(-90 14 ${m.t + ih / 2})` });
}
function hover(node, html) {
  node.addEventListener('mousemove', (e) => { tip.innerHTML = html(); tip.style.left = e.clientX + 'px'; tip.style.top = e.clientY + 'px'; tip.classList.add('on'); });
  node.addEventListener('mouseleave', () => tip.classList.remove('on'));
}
const fmt = (v, d = 0) => (Number.isFinite(v) ? v.toFixed(d) : '–');
const sgn = (v, d = 0) => (v > 0 ? '+' : v < 0 ? '−' : '') + Math.abs(v).toFixed(d);

// ------------------------------------------------------------- rounds + KPIs
function renderRounds() {
  const box = document.getElementById('rounds');
  box.innerHTML = '';
  stages.forEach((s, i) => {
    if (s.isFull) { const sp = document.createElement('span'); sp.className = 'sep'; box.appendChild(sp); }
    const b = document.createElement('button');
    b.type = 'button';
    b.className = i === cur ? 'on' : '';
    b.innerHTML = `${s.label}<span class="n">${s.ids.length} sites</span>`;
    b.onclick = () => { cur = i; render(); };
    box.appendChild(b);
  });
}
function renderKPIs() {
  const s = stages[cur], sum = s.sum, half = (sum.dofHi - sum.dofLo) / 2;
  const ok = meets(s);
  const share = s.ids.length / sites.length;
  const k = document.getElementById('kpis');
  const truthLine = showTruth ? ` · truth ${truth.dof} nm` : '';
  k.innerHTML = `
    <div class="kpi"><div class="k-l">Depth of focus at 5% EL ${ok ? '<span class="pill good">Within ±5 nm</span>' : '<span class="pill">Measuring</span>'}</div>
      <div class="k-v">${fmt(sum.dofMed)}<small>nm</small></div><div class="k-s">90% interval ${fmt(sum.dofLo)}–${fmt(sum.dofHi)} nm (±${fmt(half, 0)})${truthLine}</div></div>
    <div class="kpi"><div class="k-l">Best focus</div><div class="k-v">${sgn(sum.bestFocus)}<small>nm</small></div><div class="k-s">90% interval ${sgn(sum.bestFocusLo)} to ${sgn(sum.bestFocusHi)} nm${showTruth ? ' · truth ' + sgn((truth.f1 + truth.f2) / 2) : ''}</div></div>
    <div class="kpi"><div class="k-l">Best dose</div><div class="k-v">${fmt(sum.bestDose, 2)}<small>mJ/cm²</small></div><div class="k-s">Centre of the median window${showTruth ? ' · truth ' + pw.fromU((truth.ulo + truth.uhi) / 2).toFixed(2) : ''}</div></div>
    <div class="kpi"><div class="k-l">CD-SEM sites measured</div><div class="k-v">${s.ids.length}<small>/ ${sites.length}</small></div><div class="k-s">${Math.round(share * 100)}% of the FEM · ${s.ids.length * 3} CD readings</div></div>`;
}

// ------------------------------------------------------------- process window chart
function renderPW() {
  const s = stages[cur], d = derived(s);
  const f = frame(document.getElementById('pwChart'), 760, 500, { l: 56, r: 14, t: 28, b: 46 });
  const x = scale(-130, 130, f.m.l, f.m.l + f.iw), y = scale(25.6, 37.4, f.m.t + f.ih, f.m.t);
  const g = f.svg;
  // heatmap: one pixel per (focus, dose) cell, upscaled with smoothing
  const cv = document.createElement('canvas');
  cv.width = FOCUS_AX.length; cv.height = DOSE_AX.length;
  const cx2 = cv.getContext('2d'), img = cx2.createImageData(cv.width, cv.height);
  for (let j = 0; j < DOSE_AX.length; j++) for (let i = 0; i < FOCUS_AX.length; i++) {
    const p = d.map[j][i], o = ((DOSE_AX.length - 1 - j) * cv.width + i) * 4;
    img.data[o] = img.data[o + 1] = img.data[o + 2] = 255;
    img.data[o + 3] = Math.round(255 * (p < 0.01 ? 0 : 0.05 + 0.6 * p));
  }
  cx2.putImageData(img, 0, 0);
  axes(g, f, x, y, [-120, -80, -40, 0, 40, 80, 120], [26, 28, 30, 32, 34, 36], 'Focus (nm)', 'Dose (mJ/cm²)', (v) => (v > 0 ? '+' + v : v));
  const x0i = x(FOCUS_AX[0] - 2), x1i = x(FOCUS_AX[FOCUS_AX.length - 1] + 2), y0i = y(DOSE_AX[DOSE_AX.length - 1] + 0.1), y1i = y(DOSE_AX[0] - 0.1);
  el('image', { href: cv.toDataURL(), x: x0i, y: y0i, width: x1i - x0i, height: y1i - y0i, preserveAspectRatio: 'none' }, g);
  // window edge uncertainty (top brackets)
  if (d.edges) {
    const yb = f.m.t - 10;
    for (const [a, b] of [[d.edges.f1lo, d.edges.f1hi], [d.edges.f2lo, d.edges.f2hi]]) {
      el('line', { x1: x(a), x2: x(b), y1: yb, y2: yb, stroke: '#fff', 'stroke-width': 2, 'stroke-linecap': 'round' }, g);
      el('line', { x1: x(a), x2: x(a), y1: yb - 4, y2: yb + 4, stroke: '#fff' }, g);
      el('line', { x1: x(b), x2: x(b), y1: yb - 4, y2: yb + 4, stroke: '#fff' }, g);
    }
    txt(g, f.m.l + f.iw, f.m.t - 14, 'window edges, 90%', 'tick', 'end');
  }
  // median window
  const w = d.repWin;
  if (w.dof > 0) {
    const e1 = pw.fromU(w.ulo), e2 = pw.fromU(w.uhi);
    el('rect', { x: x(w.f1), y: y(e2), width: x(w.f2) - x(w.f1), height: y(e1) - y(e2), fill: 'none', stroke: '#fff', 'stroke-width': 1.6 }, g);
    const lab = `Median window · DOF ${fmt(s.sum.dofMed)} nm · EL ${(w.el * 100).toFixed(1)}%`;
    const lw = lab.length * 5.6 + 16;
    el('rect', { x: x(w.f1), y: y(e2) - 26, width: lw, height: 20, rx: 10, fill: '#fff' }, g);
    txt(g, x(w.f1) + 8, y(e2) - 12.5, lab, 'tick', 'start', { fill: '#000' });
  }
  if (showTruth && truth.dof > 0) {
    const e1 = pw.fromU(truth.ulo), e2 = pw.fromU(truth.uhi);
    el('rect', { x: x(truth.f1), y: y(e2), width: x(truth.f2) - x(truth.f1), height: y(e1) - y(e2), fill: 'none', stroke: '#a3a3a3', 'stroke-width': 1.2, 'stroke-dasharray': '4 4' }, g);
  }
  // sites
  const measured = new Set(s.ids);
  const fresh = new Set(s.picked || []);
  const next = new Set(d.next ? d.next.ids : []);
  for (const site of sites) {
    const cx = x(site.focus), cy = y(site.dose);
    let node;
    if (fresh.has(site.id) && !s.isFull) { el('circle', { cx, cy, r: 7.5, fill: 'none', stroke: '#fff', 'stroke-width': 1, 'stroke-opacity': 0.6 }, g); node = el('circle', { cx, cy, r: 4.6, fill: '#fff', stroke: '#000', 'stroke-width': 1.4 }, g); }
    else if (measured.has(site.id)) node = el('circle', { cx, cy, r: 4, fill: '#8a8a8a', stroke: '#000', 'stroke-width': 1.2 }, g);
    else if (next.has(site.id)) node = el('circle', { cx, cy, r: 5, fill: '#000', stroke: '#fff', 'stroke-width': 1.4, 'stroke-dasharray': '2 2' }, g);
    else node = el('circle', { cx, cy, r: 2.2, fill: '#3a3a3a' }, g);
    const hit = el('circle', { cx, cy, r: 9, fill: 'transparent' }, g);
    hover(hit, () => {
      const st = measured.has(site.id) ? (fresh.has(site.id) && !s.isFull ? 'Measured this round' : 'Measured') : next.has(site.id) ? 'Queued next' : 'Not measured';
      const o = camp.obsAll[site.id];
      const cds = measured.has(site.id) ? `<br><span class="t2">CD ${o.map((v) => v.toFixed(1)).join(' / ')} nm</span>` : '';
      return `Die R${site.row + 1}C${site.col + 1} · ${sgn(site.focus)} nm · ${site.dose.toFixed(1)} mJ/cm²<br><span class="t2">${st}</span>${cds}`;
    });
  }
  // heatmap hover via overlay
  const ov = el('rect', { x: f.m.l, y: f.m.t, width: f.iw, height: f.ih, fill: 'transparent' }, g);
  g.insertBefore(ov, g.childNodes[1]);
  hover(ov, () => '');
  ov.addEventListener('mousemove', (e) => {
    const r = f.svg.getBoundingClientRect();
    const px = (e.clientX - r.left) / r.width * f.W, py = (e.clientY - r.top) / r.height * f.H;
    const F = x.inv(px), E = y.inv(py);
    const i = Math.round((F + 122) / 4), j = Math.round((E - 26) / 0.2);
    if (i < 0 || j < 0 || i >= FOCUS_AX.length || j >= DOSE_AX.length) return;
    tip.innerHTML = `${sgn(FOCUS_AX[i])} nm · ${DOSE_AX[j].toFixed(1)} mJ/cm²<br><span class="t2">P(all features in spec) ${(d.map[j][i] * 100).toFixed(0)}%</span>`;
  });
  document.getElementById('pwLegend').innerHTML =
    `<span>P(in spec)<span class="ramp"></span></span><span><span class="sw" style="background:#8a8a8a"></span>Measured</span>` +
    (s.picked && !s.isFull ? `<span><span class="sw" style="background:#fff"></span>This round</span>` : '') +
    (d.next ? `<span><span class="sw" style="border:1.4px dashed #fff;background:#000"></span>Queued</span>` : '') +
    (showTruth ? `<span><span class="ln" style="border-color:#a3a3a3"></span>Truth</span>` : '');
}

// ------------------------------------------------------------- wafer map
function renderWafer() {
  const s = stages[cur], d = derived(s);
  const f = frame(document.getElementById('wafer'), 360, 318, { l: 0, r: 0, t: 0, b: 0 });
  const g = f.svg, R = 150, sc = 1.0, cx = 180, cy = 166;
  const X = (mm) => cx + mm * sc, Y = (mm) => cy + mm * sc;
  el('circle', { cx, cy, r: R * sc, fill: '#0b0b0b', stroke: '#2f2f2f' }, g);
  el('circle', { cx, cy, r: (R - pw.DIE.edge) * sc, fill: 'none', stroke: '#1c1c1c', 'stroke-dasharray': '2 3' }, g);
  el('path', { d: `M ${cx - 6} ${cy + R * sc} a 6 6 0 0 1 12 0`, fill: '#000', stroke: '#2f2f2f' }, g);
  const measured = new Set(s.ids), fresh = new Set(s.picked || []), next = new Set(d.next ? d.next.ids : []);
  const order = d.next ? d.next.ids : [];
  for (const site of sites) {
    const x0 = X(site.x) + 1, y0 = Y(site.y) + 1, w = pw.DIE.w * sc - 2, h = pw.DIE.h * sc - 2;
    let fill = '#161616', stroke = 'none';
    if (measured.has(site.id)) fill = fresh.has(site.id) && !s.isFull ? '#ffffff' : '#7a7a7a';
    else if (next.has(site.id)) { fill = '#000'; stroke = '#fff'; }
    const r = el('rect', { x: x0, y: y0, width: w, height: h, rx: 2, fill, stroke, 'stroke-width': 1.2 }, g);
    if (next.has(site.id)) txt(g, x0 + w / 2, y0 + h / 2 + 4, String(order.indexOf(site.id) + 1), 'tick', 'middle', { fill: '#fff', 'font-size': 10.5 });
    hover(r, () => `Die R${site.row + 1}C${site.col + 1}<br><span class="t2">${sgn(site.focus)} nm · ${site.dose.toFixed(1)} mJ/cm²</span>`);
  }
  // axis cues
  txt(g, X(-136.5) + 10, 12, '−120', 'tick', 'middle');
  txt(g, X(136.5) - 10, 12, '+120 nm', 'tick', 'middle');
  txt(g, cx, 12, 'focus', 'tick', 'middle');
  document.getElementById('waferCount').textContent = `${s.ids.length} of ${sites.length} dies measured`;
}

// ------------------------------------------------------------- queue
function renderQueue() {
  const s = stages[cur], d = derived(s);
  const box = document.getElementById('queue');
  const title = document.getElementById('queueTitle'), hint = document.getElementById('queueHint');
  const btn = document.getElementById('btnNext');
  document.getElementById('navQueue').textContent = d.next ? d.next.ids.length : '0';
  if (s.isFull) {
    title.textContent = 'Full wafer'; hint.textContent = 'Reference: every exposed die measured';
    box.innerHTML = `<div class="done"><div class="ic"><i data-lucide="layers"></i></div><div><b>${sites.length} sites, ${sites.length * 3} CD readings</b>DOF ${fmt(s.sum.dofMed)} nm (${fmt(s.sum.dofLo)}–${fmt(s.sum.dofHi)}). The adaptive run reached ${fmt(lastAd.sum.dofMed)} nm (${fmt(lastAd.sum.dofLo)}–${fmt(lastAd.sum.dofHi)}) with ${lastAd.ids.length} sites.</div></div>`;
    btn.disabled = false; btn.innerHTML = '<i data-lucide="download"></i>Export window';
    return;
  }
  if (!d.next) {
    title.textContent = 'Window qualified'; hint.textContent = 'Stop rule: 90% interval within ±5 nm';
    box.innerHTML = `<div class="done"><div class="ic"><i data-lucide="check"></i></div><div><b>No more sites needed</b>DOF ${fmt(s.sum.dofMed)} nm, 90% interval ${fmt(s.sum.dofLo)}–${fmt(s.sum.dofHi)} nm after ${s.ids.length} of ${sites.length} sites. ${sites.length - s.ids.length} dies left unmeasured.</div></div>`;
    btn.disabled = false; btn.innerHTML = '<i data-lucide="download"></i>Export window';
    return;
  }
  title.textContent = 'Next sites'; hint.textContent = 'Ranked by how much they would narrow DOF';
  const rows = d.next.ids.map((id, i) => ({ site: sites[id], score: d.next.trace[i].score }));
  const max = Math.max(...rows.map((r) => r.score), 1e-9);
  box.innerHTML = `<table class="queue"><thead><tr><th>Die</th><th>Focus</th><th>Dose</th><th class="num">Info on DOF</th></tr></thead><tbody>` +
    rows.map((r, i) => `<tr><td><span class="qi">${i + 1}</span>R${r.site.row + 1}C${r.site.col + 1}</td><td>${sgn(r.site.focus)} nm</td><td>${r.site.dose.toFixed(1)}</td><td class="num">${(r.score * 100).toFixed(0)}%<span class="bar" style="width:${Math.max(3, 46 * r.score / max)}px"></span></td></tr>`).join('') +
    `</tbody></table>`;
  btn.disabled = false; btn.textContent = `Measure next ${rows.length} sites`;
}

// ------------------------------------------------------------- DOF histogram
function renderHist() {
  const s = stages[cur];
  const f = frame(document.getElementById('dofHist'), 400, 230, { l: 40, r: 10, t: 14, b: 44 });
  const vals = Array.from(s.post.dof);
  const lo = 20, hi = 120, bw = 2;
  const nb = (hi - lo) / bw, counts = new Array(nb).fill(0);
  for (const v of vals) { const b = Math.min(nb - 1, Math.max(0, Math.floor((v - lo) / bw))); counts[b]++; }
  const mx = Math.max(...counts);
  const x = scale(lo, hi, f.m.l, f.m.l + f.iw), y = scale(0, mx * 1.1, f.m.t + f.ih, f.m.t);
  axes(f.svg, f, x, y, [20, 40, 60, 80, 100, 120], [], 'DOF at 5% EL (nm)');
  counts.forEach((c, b) => {
    if (!c) return;
    const v0 = lo + b * bw;
    const inside = v0 + bw / 2 >= s.sum.dofLo && v0 + bw / 2 <= s.sum.dofHi;
    el('rect', { x: x(v0) + 0.6, y: y(c), width: x(v0 + bw) - x(v0) - 1.2, height: y(0) - y(c), fill: inside ? '#fff' : '#3a3a3a', rx: 1 }, f.svg);
  });
  const xm = x(s.sum.dofMed);
  txt(f.svg, xm, f.m.t + 2, `median ${fmt(s.sum.dofMed)}`, 'tick', 'middle', { fill: '#fff' });
  if (showTruth) {
    el('line', { x1: x(truth.dof), x2: x(truth.dof), y1: f.m.t + 8, y2: f.m.t + f.ih, stroke: '#a3a3a3', 'stroke-dasharray': '3 3' }, f.svg);
    txt(f.svg, x(truth.dof) + 4, f.m.t + 22, 'truth', 'tick', 'start');
  }
}

// ------------------------------------------------------------- ED curve
function renderED() {
  const s = stages[cur], d = derived(s);
  const f = frame(document.getElementById('edChart'), 400, 230, { l: 40, r: 10, t: 14, b: 44 });
  const x = scale(0, 15, f.m.l, f.m.l + f.iw), y = scale(0, 140, f.m.t + f.ih, f.m.t);
  axes(f.svg, f, x, y, [0, 5, 10, 15], [0, 40, 80, 120], 'Exposure latitude (%)', null, (v) => v + '%');
  const L = pw.EL_LEVELS.map((e) => e * 100);
  let area = `M ${x(L[0])} ${y(d.edBand[0].hi)}`;
  L.forEach((v, i) => { area += ` L ${x(v)} ${y(d.edBand[i].hi)}`; });
  for (let i = L.length - 1; i >= 0; i--) area += ` L ${x(L[i])} ${y(d.edBand[i].lo)}`;
  el('path', { d: area + ' Z', fill: '#fff', 'fill-opacity': 0.12 }, f.svg);
  el('path', { d: 'M ' + L.map((v, i) => `${x(v)} ${y(d.edBand[i].mid)}`).join(' L '), fill: 'none', stroke: '#fff', 'stroke-width': 1.6 }, f.svg);
  if (showTruth) el('path', { d: 'M ' + L.map((v, i) => `${x(v)} ${y(truthED[i])}`).join(' L '), fill: 'none', stroke: '#a3a3a3', 'stroke-width': 1.2, 'stroke-dasharray': '4 3' }, f.svg);
  el('line', { x1: x(5), x2: x(5), y1: f.m.t, y2: f.m.t + f.ih, stroke: '#2f2f2f', 'stroke-dasharray': '2 3' }, f.svg);
  el('circle', { cx: x(5), cy: y(s.sum.dofMed), r: 3.5, fill: '#fff' }, f.svg);
  txt(f.svg, x(5) + 7, y(s.sum.dofMed) - 7, `${fmt(s.sum.dofMed)} nm`, 'tick', 'start', { fill: '#fff' });
  txt(f.svg, f.m.l - 8, f.m.t - 2, 'nm', 'tick', 'end');
}

// ------------------------------------------------------------- Bossung
const DOSES = [29.3, 30.4, 31.5, 32.6, 33.7];
const GREYS = ['#4a4a4a', '#7a7a7a', '#ffffff', '#b5b5b5', '#5e5e5e'];
function renderBossung() {
  const s = stages[cur];
  const tabs = document.getElementById('featTabs');
  tabs.innerHTML = '';
  pw.FEATURES.forEach((ft, k) => {
    const b = document.createElement('button');
    b.type = 'button'; b.className = k === feat ? 'on' : ''; b.textContent = ['Dense', 'Semi', 'Iso'][k];
    b.onclick = () => { feat = k; renderBossung(); };
    tabs.appendChild(b);
  });
  document.getElementById('bossHint').textContent = pw.FEATURES[feat].detail;
  const f = frame(document.getElementById('bossung'), 500, 230, { l: 40, r: 12, t: 14, b: 44 });
  const x = scale(-130, 130, f.m.l, f.m.l + f.iw), y = scale(24, 50, f.m.t + f.ih, f.m.t);
  const cid = 'bclip';
  const cp = el('clipPath', { id: cid }, el('defs', {}, f.svg));
  el('rect', { x: f.m.l, y: f.m.t, width: f.iw, height: f.ih }, cp);
  // spec band
  el('rect', { x: f.m.l, y: y(44), width: f.iw, height: y(36) - y(44), fill: '#fff', 'fill-opacity': 0.045 }, f.svg);
  axes(f.svg, f, x, y, [-120, -60, 0, 60, 120], [28, 36, 40, 44, 48], 'Focus (nm)', null, (v) => (v > 0 ? '+' + v : v));
  const plot = el('g', { 'clip-path': `url(#${cid})` }, f.svg);
  txt(f.svg, f.m.l - 8, f.m.t - 2, 'CD nm', 'tick', 'end');
  const focusAx = Array.from({ length: 53 }, (_, i) => -130 + 5 * i);
  DOSES.forEach((E, di) => {
    const band = pw.bossungBand(s.post.draws, feat, E, focusAx);
    if (E === 31.5) {
      let a = 'M ' + band.map((b) => `${x(b.F)} ${y(b.hi)}`).join(' L ');
      a += ' L ' + band.slice().reverse().map((b) => `${x(b.F)} ${y(b.lo)}`).join(' L ');
      el('path', { d: a + ' Z', fill: '#fff', 'fill-opacity': 0.1 }, plot);
    }
    el('path', { d: 'M ' + band.map((b) => `${x(b.F)} ${y(b.mid)}`).join(' L '), fill: 'none', stroke: GREYS[di], 'stroke-width': E === 31.5 ? 1.7 : 1.2 }, plot);
    // measured points at this dose
    for (const id of s.ids) {
      const site = sites[id];
      if (Math.abs(site.dose - E) > 1e-6) continue;
      const v = camp.obsAll[id][feat];
      const pt = el('circle', { cx: x(site.focus), cy: y(v), r: 2.8, fill: GREYS[di], stroke: '#080808', 'stroke-width': 1 }, plot);
      hover(pt, () => `${pw.FEATURES[feat].name} · ${sgn(site.focus)} nm · ${E.toFixed(1)} mJ/cm²<br><span class="t2">measured CD ${v.toFixed(1)} nm</span>`);
    }
  });
  document.getElementById('bossLegend').innerHTML = DOSES.map((E, di) => `<span><span class="lnk" style="border-color:${GREYS[di]}"></span>${E.toFixed(1)}</span>`).join('') + '<span class="u">mJ/cm²</span>';
}

// ------------------------------------------------------------- validation
let BT = null;
function renderValidation() {
  if (!BT) return;
  const sm = BT.summary;
  document.getElementById('valHint').textContent = `${BT.N} simulated tools, same FEM and spec; error = |estimated − true DOF| at 5% EL`;
  const f = frame(document.getElementById('valChart'), 420, 230, { l: 48, r: 14, t: 14, b: 44 });
  const x = scale(14, 42, f.m.l, f.m.l + f.iw), y = scale(0, 6, f.m.t + f.ih, f.m.t);
  axes(f.svg, f, x, y, [16, 24, 32, 40], [0, 2, 4, 6], 'Sites measured', 'Mean DOF error (nm)');
  const fullKey = Object.keys(sm).find((k) => k.startsWith('full@'));
  const yf = y(sm[fullKey].meanErr);
  el('line', { x1: f.m.l, x2: f.m.l + f.iw, y1: yf, y2: yf, stroke: '#a3a3a3', 'stroke-dasharray': '4 4' }, f.svg);
  txt(f.svg, f.m.l + 4, yf + 15, `Full wafer, all ${sm[fullKey].sites} sites`, 'tick', 'start');
  const series = [['adaptive', '#ffffff', 2], ['grid', '#7a7a7a', 1.4], ['random', '#4a4a4a', 1.4]];
  for (const [st, col, wdt] of series) {
    const pts = [16, 24, 32, 40].map((b) => { const k = (b === 16 ? 'adaptive' : st) + '@' + b; return [b, sm[k].meanErr]; });
    el('path', { d: 'M ' + pts.map(([a, b]) => `${x(a)} ${y(b)}`).join(' L '), fill: 'none', stroke: col, 'stroke-width': wdt }, f.svg);
    pts.forEach(([a, b]) => el('circle', { cx: x(a), cy: y(b), r: st === 'adaptive' ? 3.2 : 2.4, fill: col }, f.svg));
  }
  const name = { adaptive: 'Adaptive', grid: 'Space-filling', random: 'Seed + random', full: 'Full wafer' };
  const col = { adaptive: '#fff', grid: '#7a7a7a', random: '#4a4a4a', full: '#a3a3a3' };
  const rowsK = ['adaptive@40', 'grid@40', 'random@40', fullKey];
  document.getElementById('valTable').innerHTML =
    `<thead><tr><th>Strategy</th><th>Sites</th><th>Mean error</th><th>90th pct error</th><th>90% interval covers truth</th></tr></thead><tbody>` +
    rowsK.map((k) => { const r = sm[k]; const st = r.strategy; return `<tr class="${st === 'adaptive' ? 'me' : ''}"><td><span class="key" style="border-color:${col[st]};${st === 'full' ? 'border-top-style:dashed' : ''}"></span>${name[st]}</td><td>${r.sites}</td><td>${r.meanErr.toFixed(1)} nm</td><td>${r.p90Err.toFixed(1)} nm</td><td>${Math.round(r.coverage90 * 100)}%</td></tr>`; }).join('') +
    `</tbody>`;
}

// ------------------------------------------------------------- wiring
function render() {
  renderRounds(); renderKPIs(); renderPW(); renderWafer(); renderQueue(); renderHist(); renderED(); renderBossung(); renderValidation();
  if (window.lucide) window.lucide.createIcons();
}
function exportWindow() {
  const s = stages[cur], d = derived(s), w = d.repWin;
  const report = {
    tool: 'LT-01', wafer: 'W-0418', stage: s.label, sitesMeasured: s.ids.length, sitesOnWafer: sites.length,
    spec: { targetNm: pw.SPEC.target, tol: pw.SPEC.tol }, requiredEL: EL_REQ,
    dofNm: { median: s.sum.dofMed, lo90: s.sum.dofLo, hi90: s.sum.dofHi },
    bestFocusNm: s.sum.bestFocus, bestDose: +s.sum.bestDose.toFixed(3),
    medianWindow: { focusNm: [w.f1, w.f2], dose: [+pw.fromU(w.ulo).toFixed(3), +pw.fromU(w.uhi).toFixed(3)], el: +w.el.toFixed(4) },
    measured: s.ids.map((id) => ({ die: `R${sites[id].row + 1}C${sites[id].col + 1}`, focusNm: sites[id].focus, dose: sites[id].dose, cdNm: camp.obsAll[id].map((v) => +v.toFixed(2)) })),
    note: 'Synthetic data. Unofficial concept demo by Mohith Gajjela; not affiliated with Source Foundry.',
  };
  const a = document.createElement('a');
  a.href = URL.createObjectURL(new Blob([JSON.stringify(report, null, 2)], { type: 'application/json' }));
  a.download = `process-window-LT-01-W-0418-${s.label.toLowerCase().replace(/\s+/g, '-')}.json`;
  a.click();
}
document.getElementById('btnNext').onclick = () => {
  const s = stages[cur];
  if (!s.isFull && cur < stages.length - 2) { cur++; render(); } else exportWindow();
};
document.getElementById('btnTruth').onclick = (e) => { showTruth = !showTruth; e.currentTarget.classList.toggle('on', showTruth); render(); };
fetch('data/backtest.json').then((r) => r.json()).then((j) => { BT = j; renderValidation(); }).catch((err) => console.error(err));
render();
