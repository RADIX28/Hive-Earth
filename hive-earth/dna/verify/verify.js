// Checks Helix Lab's engine against Biopython reference data and known statistics.
// Usage (from hive-earth/dna):
//   python -m http.server 8799 &
//   python verify/reference.py verify/reference.json        (needs: pip install biopython)
//   node verify/verify.js verify/reference.json [url]        (needs: npm i playwright)
const { chromium } = require('playwright');
const fs = require('fs');
const ref = JSON.parse(fs.readFileSync(process.argv[2] || 'verify/reference.json'));
const url = process.argv[3] || 'http://localhost:8799/index.html';
(async () => {
  const b = await chromium.launch(); const p = await b.newPage();
  await p.route('**/*', r => r.request().url().startsWith(new URL(url).origin) ? r.continue() : r.abort());
  await p.goto(url); await p.waitForTimeout(500);
  const checks = await p.evaluate(ref => {
    const C = [], ok = (name, pass, detail) => C.push({ name, pass: !!pass, detail });
    const map = { 1: 'std', 2: 'mito', 6: 'cil' };
    for (const k in map) ok('NCBI table ' + k + ' (' + map[k] + ')', G.CODES[map[k]].aa === ref.tables[k].aa, '');
    ok('PL0344 code: only TAA/TAG differ from standard', [...STD].every((c, i) => i === 10 || i === 11 ? true : G.CODES.pl.aa[i] === c) && G.tr('TAA', 'pl') === 'K' && G.tr('TAG', 'pl') === 'E' && G.tr('TGA', 'pl') === '*', '');
    let tr = 0, rc = 0, orf = 0, orf6 = 0;
    const key = o => o.map(x => [x.strand === '+' ? '+' : '-', x.frame, x.start, x.aa]).map(x => x.join()).sort().join(';');
    const rkey = o => o.map(x => x.join()).sort().join(';');
    for (const d of ref.dna) {
      if (G.revComp(d.s) !== d.rc) rc++;
      for (const k in map) for (let f = 0; f < 3; f++) if (G.translate(d.s, f, map[k], false).protein !== d.tr[k][f]) tr++;
      if (key(G.orfs(d.s, 'std', 10)) !== rkey(d.orfs)) orf++;
      if (key(G.orfs(d.s, 'cil', 8)) !== rkey(d.orfs6)) orf6++;
    }
    ok('Translation, 3 frames x 3 tables', tr === 0, tr + ' mismatches of ' + ref.dna.length * 9);
    ok('Reverse complement', rc === 0, rc + ' mismatches');
    ok('Six-frame ORF finder (standard, ciliate)', orf + orf6 === 0, orf + '+' + orf6 + ' mismatches');
    let mw = 0, gr = 0, q = 0, piBad = 0;
    for (const x of ref.prot) {
      const s = G.stats(x.p);
      mw = Math.max(mw, Math.abs(s.kDa * 1000 - x.mw) / x.mw); gr = Math.max(gr, Math.abs(s.gravy - x.gravy)); q = Math.max(q, Math.abs(s.charge - x.q7));
      // Biopython only searches pI in 4.05-12; outside that its answer is a clamp, so compare where it is a true root.
      if (x.pI > 4.06 && x.pI < 11.99 && Math.abs(s.pI - x.pI) > 0.01) piBad++;
    }
    ok('Molecular weight (rel. error < 0.01%)', mw < 1e-4, (mw * 100).toFixed(4) + '%');
    ok('GRAVY (Kyte-Doolittle)', gr < 1e-9, gr);
    ok('Net charge at pH 7 (ExPASy pKa set)', q < 1e-6, q);
    ok('Isoelectric point (within 0.01)', piBad === 0, piBad + ' mismatches');
    const c = { syn: 0, mis: 0, non: 0 };
    for (const a of 'TCAG') for (const b2 of 'TCAG') for (const d2 of 'TCAG') { const cod = a + b2 + d2, aa = G.tr(cod); if (aa === '*') continue;
      for (let i = 0; i < 3; i++) for (const nb of 'TCAG') if (nb !== cod[i]) { const m = G.tr(cod.slice(0, i) + nb + cod.slice(i + 1)); c[m === '*' ? 'non' : m === aa ? 'syn' : 'mis']++; } }
    ok('Single-base changes: 134 silent / 392 missense / 23 nonsense', c.syn === ref.subst.syn && c.mis === ref.subst.mis && c.non === ref.subst.non, JSON.stringify(c));
    const tm = k => G.tmSegments(PROT[k]).map(x => [x.s + 1, x.e]);
    ok('Sarcolipin TM helix = 8-26', JSON.stringify(tm('SLN')) === '[[8,26]]', JSON.stringify(tm('SLN')));
    ok('Phospholamban TM helix inside 28-52', tm('PLN').length === 1 && tm('PLN')[0][0] >= 28 && tm('PLN')[0][1] <= 52, JSON.stringify(tm('PLN')));
    ok('No TM helix in soluble HBB / INS / GFP', !tm('HBB').length && !tm('INS').length && !tm('GFP').length, '');
    const hn = PRESETS.find(x => x.id === 'humanin').seq(), mc = PRESETS.find(x => x.id === 'motsc').seq();
    ok('Humanin: 24 aa (standard) / 21 aa (mitochondrial)', G.cds(hn).protein === PROT.HN + '*' && G.cds(hn, 'mito').protein === 'MAPRGFSCLLLLTSEIDLPVK*', '');
    ok('MOTS-c ORF = MRWQEMGYIFYPRKLR (51 nt)', mc.length === 51 && G.cds(mc).protein === PROT.MOTSC + '*', '');
    const g = PROT.GFP, hb = PROT.HBB;
    ok('GFP chromophore S65-Y66-G67, T203', g.slice(64, 67) === 'SYG' && g[202] === 'T', '');
    ok('HBB Glu6 / Glu26 (clinical numbering)', hb[6] === 'E' && hb[26] === 'E', '');
    ok('Goldman 2013 rotation table', ['A', 'C', 'G', 'T'].map(x => Secret.rot(x).join('')).join() === 'CGT,GTA,TAC,ACG', '');
    ok('Secret codes round-trip (incl. emoji)', ['bin', 'rot'].every(s => ['Hi 🧬', 'Ünïcødé', ''].every(m => Secret.decodeWith(s, Secret.encodeWith(s, m)) === m)), '');
    // random-DNA ORF lengths follow the geometric model
    const runs = []; for (let n = 0; n < 200; n++) { let s = ''; for (let i = 0; i < 3000; i++) s += 'ACGT'[(Math.random() * 4) | 0];
      for (let i = s.indexOf('ATG'); i >= 0; i = s.indexOf('ATG', i + 1)) { let k = 1, j = i + 3; for (; j + 3 <= s.length; j += 3, k++) if (G.tr(s.substr(j, 3)) === '*') break; if (j + 3 <= s.length) runs.push(k); } }
    const pS = G.stopProb('', 'std'), s50 = runs.filter(x => x >= 50).length / runs.length;
    ok('Random ORFs >= 50 codons match (61/64)^49', Math.abs(s50 - G.orfLuck(50, pS)) < 0.015, 'sim ' + s50.toFixed(4) + ' vs ' + G.orfLuck(50, pS).toFixed(4));
    // pond mutation process
    const TS = { AG: 1, GA: 1, CT: 1, TC: 1 }, g0 = Pond.founderGenome(), L = GENE_KEYS.reduce((a, k) => a + g0[k].length, 0);
    let subs = 0, ts = 0, tv = 0; const N = 15000;
    for (let n = 0; n < N; n++) { const m = Pond.mutate(g0); for (const k of GENE_KEYS) if (m[k].length === g0[k].length) for (let i = 0; i < m[k].length; i++) if (m[k][i] !== g0[k][i]) { subs++; TS[g0[k][i] + m[k][i]] ? ts++ : tv++; } }
    const nominal = Pond.mut / 1000 * L, got = subs / N;
    ok('Pond substitution rate = UI rate (±10%)', Math.abs(got - nominal) / nominal < 0.1, got.toFixed(3) + ' vs ' + nominal.toFixed(3) + ' per baby');
    ok('Pond Ti/Tv ≈ 2 (human-like)', Math.abs(ts / tv - 2) < 0.2, (ts / tv).toFixed(2));
    const T = computeTraits(g0), silent = Object.assign({}, g0); silent.size = silent.size.replace(/GGC/, 'GGT');
    ok('Silent mutation leaves traits unchanged', JSON.stringify(computeTraits(silent)) === JSON.stringify(T), '');
    return C;
  }, ref);
  let fail = 0;
  for (const c of checks) { if (!c.pass) fail++; console.log((c.pass ? 'PASS ' : 'FAIL ') + c.name + (c.detail !== '' ? '  [' + c.detail + ']' : '')); }
  console.log(`\n${checks.length - fail}/${checks.length} checks passed`);
  await b.close(); process.exit(fail ? 1 : 0);
})();
