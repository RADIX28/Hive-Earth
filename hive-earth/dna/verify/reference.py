# Builds reference answers with Biopython (NCBI codon tables, ProtParam-style protein maths).
# Usage: python reference.py reference.json
import json, random, sys
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.SeqUtils.ProtParam import ProteinAnalysis
random.seed(2026)
B = 'TCAG'
out = {'tables': {}, 'dna': [], 'prot': [], 'subst': {}}
for tid in (1, 2, 6):
    t = CodonTable.unambiguous_dna_by_id[tid]
    s = ''.join(t.forward_table.get(a+b+c, '*') if a+b+c not in t.stop_codons else '*' for a in B for b in B for c in B)
    out['tables'][tid] = {'aa': s, 'starts': t.start_codons}
# single-substitution classes for the standard table
t = CodonTable.unambiguous_dna_by_id[1]
syn = mis = non = 0; ti = {'syn':0,'mis':0,'non':0}
TS = {('A','G'),('G','A'),('C','T'),('T','C')}
for a in B:
  for b in B:
    for c in B:
      cod = a+b+c
      if cod in t.stop_codons: continue
      aa = t.forward_table[cod]
      for i in range(3):
        for nb in B:
          if nb == cod[i]: continue
          m = cod[:i]+nb+cod[i+1:]
          k = 'non' if m in t.stop_codons else ('syn' if t.forward_table[m] == aa else 'mis')
          if k=='syn': syn+=1
          elif k=='mis': mis+=1
          else: non+=1
          if (cod[i], nb) in TS: ti[k]+=1
out['subst'] = {'syn': syn, 'mis': mis, 'non': non, 'ti': ti}
def ref_orfs(s, table, mn):
    res = []
    for strand, ss in (('+', s), ('-', str(Seq(s).reverse_complement()))):
        for f in range(3):
            prot = str(Seq(ss[f:f + (len(ss)-f)//3*3]).translate(table=table))
            i = 0; open_ = -1
            for k, aa in enumerate(prot):
                cod = ss[f+3*k:f+3*k+3]
                if open_ < 0:
                    if cod == 'ATG': open_ = k
                elif aa == '*':
                    n = k - open_
                    if n >= mn: res.append([strand, f, f+3*open_, n])
                    open_ = -1
    return sorted(res)
for n in range(250):
    L = random.randint(20, 420)
    s = ''.join(random.choice('ACGT') for _ in range(L))
    d = {'s': s, 'rc': str(Seq(s).reverse_complement()), 'tr': {}, 'orfs': ref_orfs(s, 1, 10), 'orfs6': ref_orfs(s, 6, 8)}
    for tid in (1, 2, 6):
        d['tr'][tid] = [str(Seq(s[f:f + (L-f)//3*3]).translate(table=tid)) for f in range(3)]
    out['dna'].append(d)
AAS = 'ACDEFGHIKLMNPQRSTVWY'
presets = {
 'HBB': 'MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH',
 'INS': 'MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN',
 'GFP': 'MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK',
 'HN': 'MAPRGFSCLLLLTSEIDLPVKRRA', 'MOTSC': 'MRWQEMGYIFYPRKLR',
 'SLN': 'MGINTRELFLNFTIVLITVILMWLLVRSYQY', 'PLN': 'MEKVQYLTRSAIRRASTIEMPQQARQKLQNLFINFCLILICLLLICIIVMLL'}
plist = list(presets.items()) + [('rand%d' % i, 'M' + ''.join(random.choice(AAS) for _ in range(random.randint(5, 300)))) for i in range(150)]
for name, p in plist:
    pa = ProteinAnalysis(p)
    out['prot'].append({'name': name, 'p': p, 'mw': pa.molecular_weight(), 'gravy': pa.gravy(), 'pI': pa.isoelectric_point(), 'q7': pa.charge_at_pH(7.0)})
json.dump(out, open(sys.argv[1], 'w'))
print('tables', {k: v['aa'] for k, v in out['tables'].items()})
print('subst', out['subst'])
