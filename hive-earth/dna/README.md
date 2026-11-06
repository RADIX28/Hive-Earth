# Helix Lab: a DNA game

A playful DNA simulator in one self-contained file (`index.html`). Double-click
**Open Helix Lab.bat** on Windows, or just open `index.html` in a browser.

| Room | What you do |
|---|---|
| 🧪 **DNA Lab** | Rewrite DNA letter by letter and watch DNA → RNA → protein update live. Every edit is classified (silent, missense, nonsense, frameshift, start/stop lost…). Switch between four genetic-code "dialects", scan all six reading frames for hidden genes and microproteins, fold the protein as a cartoon, or view it in real 3D. |
| 🐟 **Evolution Pond** | Creatures whose colour, size, speed and heat tolerance come from their translated genes. Change the floor, temperature, food, predators and mutation rate, trigger heat waves and bottlenecks, and watch traits evolve on the charts. Junk DNA can give birth to brand-new genes. Open any creature in the Lab, edit it, and release it back. |
| 🔐 **Secret Code** | Encode text into DNA (2-bit, Goldman-style no-repeat, or "protein words"), damage it, and decode it, with optional triple-copy error repair. |
| 📚 **Library** | 37 cited, beginner-friendly cards, from the basics through the dark proteome and microproteins to discoveries up to November 2026. Live lookups from UniProt / Ensembl / AlphaFold and a newest-papers feed from Europe PMC. |
| ⭐ **Discoveries** | 25 badges, each a real piece of biology. |

**Offline vs online.** Everything works offline except the live features: UniProt
and Europe PMC lookups, real gene DNA from Ensembl, AlphaFold and ESMFold 3D
structures (3Dmol.js loads from jsDelivr), and the Google Fonts.

**Honesty notes.** Offline presets spell real protein sequences with common human
codons; "Get the real DNA" fetches the true coding sequence from Ensembl. The Pond
uses simplified *game rules* for how proteins become traits, and the "Will it still
work?" meter is a playful estimate, not a prediction.

Progress is saved in your browser (localStorage) and stays on your device.

## Verification

`verify/` checks the simulator against independent reference data (25 checks, all passing):

- **Genetic code:** standard, vertebrate-mitochondrial and ciliate tables match NCBI (via Biopython)
  exactly; 2,250 translations, reverse complements and six-frame ORF scans of random DNA agree.
  The PL0344 dialect differs from standard only at TAA→Lys and TAG→Glu (McGowan et al. 2023).
- **Protein maths:** molecular weight is within 0.002% of Biopython (HBB 15,998 Da, insulin 11,981 Da,
  GFP 26,886 Da). GRAVY is exact. Charge and pI use the ExPASy/Bjellqvist pKa method and match
  Biopython (HBB pI 6.74, insulin 5.22, GFP 5.67).
- **Known biology:** the 549 single-base changes split 134 silent / 392 missense / 23 nonsense.
  Membrane helices are predicted at sarcolipin 8–26 and phospholamban 31–49 (known 8–26 and 31–52).
  The real humanin ORF gives 24 aa with the standard code and 21 aa with the mitochondrial code.
  MOTS-c gives MRWQEMGYIFYPRKLR.
- **Statistics:** random-DNA ORF lengths follow the geometric model (STOP every ~21 codons). The pond
  mutates at exactly the rate shown, with transitions twice as common as transversions
  (Ti/Tv ≈ 2, as in humans) and ~1 in 10 changes an insertion/deletion.

Run it from this folder:

```
python -m http.server 8799 &
pip install biopython && python verify/reference.py verify/reference.json
npm i playwright && node verify/verify.js verify/reference.json
```
