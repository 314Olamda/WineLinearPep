# 🔗 Wine LinearPep — De novo Linear Peptide Identification (a/b/y + Immonium Scoring)

> *Step 3b of the Wine Peptidome series — mass-directed composition search, a/b/y/immonium fragment scoring, and target-decoy FDR validation for linear peptides (2–50 AA, optimized for 2–10 AA) from wine lees LC-MS/MS data.*

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue?logo=python&logoColor=white)](https://www.python.org/) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE) [![Pyteomics](https://img.shields.io/badge/Pyteomics-mass%20engine-8B0000)](https://pyteomics.readthedocs.io/) [![ORCID](https://img.shields.io/badge/ORCID-0000--0002--7720--3733-a6ce39?logo=orcid)](https://orcid.org/0000-0002-7720-3733)

---

## 🍷 Where this fits in the Wine Peptidome series

This is **Step 3b** of the programmatic pipeline, sitting alongside WineCycloPep as the linear counterpart to cyclic peptide detection:

| Step | Repository | What it does |
|------|------------|---------------|
| 1 | [WinePeptidome](https://github.com/314Olamda/WinePeptidome) | Retrieves *S. cerevisiae* & *V. vinifera* proteins (500 Da – 100 kDa) from UniProt REST API + Proteins API |
| 2 | [WineStructure](https://github.com/314Olamda/WineStructure) | AlphaFold 3D structures + per-residue pLDDT confidence |
| 3 | [WineCycloPep](https://github.com/314Olamda/WineCycloPep) | De novo **cyclic** peptide detection from Bruker `.d` files — bn-ion rotations, FDR, 3D conformers |
| **3b** | **WineLinearPep ← you are here** | De novo **linear** peptide identification — a/b/y/immonium scoring, target-decoy FDR |

---

## 🔬 Why a separate pipeline for linear peptides

Linear and cyclic peptides of identical composition are indistinguishable by precursor mass alone — the distinction is entirely in the fragmentation pattern. WineCycloPep's scoring logic (bn-ion rotations, y1 *absence*, residue-loss ions) is specific to head-to-tail cyclic topology and doesn't transfer directly. Linear peptides — the dominant form of dipeptides and short bioactive fragments annotated in wine lees LMW fractions (e.g. Ala-Phe, Gly-Leu, Ala-Ala) — need their own scoring logic built around the fragmentation physics that actually applies to them:

| | Cyclic (WineCycloPep) | Linear (this repo) |
|---|---|---|
| Fragmentation series | *n* ring-opening rotations → *n* overlapping b-ion ladders | One b-ion ladder, one y-ion ladder |
| y-ion presence | **Absence** of y1 is diagnostic (no free C-terminus) | **Presence** of y-ions is diagnostic (free C-terminus) |
| a-ions | Redundant with the rotated b-ladder — not scored separately | Independent evidence — its own scoring dimension |
| Sequence permutations | Rotational + reverse equivalents collapse to one candidate | Every permutation is a genuinely distinct molecule (IAA ≠ AAI ≠ AIA) |
| Immonium ions | Composition confirmation only, weight 0.10 | Same role, same weight |

---

## 🧬 Pipeline architecture

```
graph TD
    A[MS2 spectrum: mz_obs, int_obs] --> B[Precursor neutral mass\nlinear: sum residues + H2O]
    B --> C[compositions_for_linear_mass\nbranch-and-bound search]
    C --> D[all_permutations_linear\nevery ordering is distinct]

    D --> E[score_spectrum_vs_linear]
    E --> E1[b-ion coverage\nweight 0.30]
    E --> E2[y-ion coverage\nweight 0.30]
    E --> E3[a-ion coverage\nweight 0.15]
    E --> E4[intensity score\nweight 0.15]
    E --> E5[immonium score\nweight 0.10]

    E1 & E2 & E3 & E4 & E5 --> F[composite score\nconfidence tier]

    F --> G{score >= min_score?}
    G -->|No| H[discard]
    G -->|Yes| I[candidate retained]

    I --> J[generate_decoys_linear\npermutations excluding target\n+ isobaric I/L fallback]
    J --> K[global decoy pool\nempirical p-value]
    K --> L[BH-FDR correction\nStorey q-value\nHotelling T2 / F-stat]

    L --> M[(ranked candidates)]
    M --> N[sequence, composite,\np_value, q_storey, rejected_BH]
```

---

## ⚡ Quick start

```bash
# 1. Clone
git clone https://github.com/314Olamda/WineLinearPep.git
cd WineLinearPep

# 2. Install dependencies
pip install pyteomics numpy scipy

# 3. Run the self-test (synthetic spectrum, no data needed)
python linear_denovo.py
```

```python
# 4. Use on your own spectrum
from linear_denovo import identify_linear_peptide

candidates = identify_linear_peptide(
    mz_obs, int_obs,          # your MS2 spectrum arrays
    precursor_mz=237.1234,
    charge=1,
    n_min=2, n_max=10,         # 2-10 AA search window
    tol=0.02,
    min_score=0.10,
    n_decoys=100,
    fdr_alpha=0.05,
)

for c in candidates[:10]:
    print(c.sequence, c.composite, c.q_storey, c.rejected_BH)
```

Results are sorted by `q_storey` ascending (best first) — use it as your primary reporting metric in publications, same convention as WineCycloPep.

---

## 📦 Repository structure

| File | Role |
|------|------|
| `mass_utils.py` | Pyteomics-backed mass/ion layer — monoisotopic residue masses (20 canonical AA only, explicitly excluding Pyteomics' J/U/O ambiguity codes), neutral mass, b/loss/immonium ion calculators |
| `linear_peptide_scoring.py` | `score_spectrum_vs_linear()` — a/b/y/immonium coverage scoring for one candidate sequence against one spectrum |
| `linear_denovo.py` | Full pipeline — composition search, permutation enumeration, decoy generation, BH-FDR / Storey q-value / Hotelling T² |

---

## ⚙️ Configuration

Key parameters, all exposed as function arguments in `identify_linear_peptide()`:

```python
n_min          = 2      # minimum residues (2-10 AA range validated)
n_max          = 10     # maximum residues
tol            = 0.02   # fragment mass tolerance, Da
min_score      = 0.10   # minimum composite score to report
max_perms_per_composition = 200   # cap on permutations per composition
n_decoys       = 100    # decoy shuffles per candidate (FDR)
fdr_alpha      = 0.05   # FDR significance threshold
allowed_aa     = None   # optionally restrict the AA alphabet
```

Scoring weights (sum to 1.0), tuned for linear-specific fragmentation logic:

```
composite = 0.30 x b_coverage + 0.30 x y_coverage + 0.15 x a_coverage +
            0.15 x intensity_score + 0.10 x immonium_score
```

| Score | Tier | Interpretation |
|-------|------|-----------------|
| ≥ 0.60 | **HIGH** | Strong evidence across b, y, and a series |
| 0.35–0.59 | **MEDIUM** | Probable; validate with authentic standard |
| 0.15–0.34 | **LOW** | Possible; check for isobaric alternatives |
| < 0.15 | VERY LOW | Discard or review manually |

---

## 📊 Statistical validation

Identical framework to WineCycloPep — target-decoy competition adapted from shotgun proteomics (Elias & Gygi 2007):

- **Empirical p-value**: one-tailed, against a global decoy pool, with a Laplace pseudocount to avoid p = 0
- **Benjamini-Hochberg FDR** (1995): standard step-up correction
- **Storey q-value** (2003): π₀-adjusted, less conservative when many true signals are present — **use `q_storey` for publication reporting**
- **Hotelling T² / F-statistic**: tests whether the 5-dimensional score vector (b/y/a/intensity/immonium jointly) separates targets from decoys, with per-dimension Mann-Whitney U + BH correction identifying which criteria drive the separation

Decoys are generated by permuting the candidate's own composition, excluding the candidate itself — with an isobaric I↔L substitution fallback (shared monoisotopic mass 113.084 Da) for short or compositionally degenerate sequences where too few unique permutations exist.

---

## ⚠️ Known limitation: isobaric residues

I and L are isobaric — identical monoisotopic mass, indistinguishable by accurate-mass MS2 alone. When a composition contains I or L, expect tied top candidates (e.g. `IAA` and `LAA` scoring identically). This is reported correctly by the pipeline rather than arbitrarily resolved — disambiguation requires orthogonal evidence (retention time, ion mobility 1/K₀ from timsTOF PASEF, or chemical derivatization), consistent with the literature on short-peptide identification limits below ~5 residues.

---

## 🔗 Series & related resources

- **Step 3 (cyclic):** [WineCycloPep](https://github.com/314Olamda/WineCycloPep)
- [Pyteomics](https://pyteomics.readthedocs.io/) — mass calculation engine underlying this pipeline
- [GNPS molecular networking](https://gnps.ucsd.edu/) — downstream spectral annotation
- [reLees project](https://relees.uniwa.gr) — wine lees circular economy research

---

## 📄 Citation

```bibtex
@software{gimenez_gil_wine_linearpep_2025,
  author  = {Giménez-Gil, Pol},
  title   = {Wine LinearPep: De novo Linear Peptide Identification with a/b/y/Immonium Scoring},
  year    = {2025},
  url     = {https://github.com/314Olamda/WineLinearPep},
  orcid   = {0000-0002-7720-3733},
  note    = {Step 3b of the Wine Peptidome series. Cyclic counterpart: github.com/314Olamda/WineCycloPep}
}
```

---

## 👤 Author

**Pol Giménez-Gil**, PhD
Postdoctoral Researcher — ISVV, Université de Bordeaux
Scopus ID: 57219336109 · ORCID: [0000-0002-7720-3733](https://orcid.org/0000-0002-7720-3733)
ResearchGate: [Pol_Gimenez2](https://www.researchgate.net/profile/Pol_Gimenez2)

---

## 📜 License

MIT — see [LICENSE](LICENSE)
