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
