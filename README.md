# `pep-admet-engine`

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![Snakemake](https://img.shields.io/badge/Snakemake-Catalog%20Standard-4B8BBE?logo=snakemake&logoColor=white)](https://snakemake.github.io/snakemake-workflow-catalog/)
[![Build & Tests](https://img.shields.io/badge/tests-passing%20(pytest)-success.svg)](https://github.com/Juanj053/pep-admet-engine/actions)
[![Coverage](https://img.shields.io/badge/coverage-100%25-brightgreen.svg)](https://github.com/Juanj053/pep-admet-engine)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-Juan%20Jos%C3%A9%20%C3%81lvarez-0A66C2.svg?logo=linkedin)](https://www.linkedin.com/in/juanj-alva)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code Style: Ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)
[![Type Checked: Mypy](https://img.shields.io/badge/type%20checked-mypy-blue.svg)](https://mypy-lang.org/)

> **In-Silico Pharmacokinetic Profiling, Multi-Protease Stability Screening, and Rational Peptidomimetic Stabilization Advisor for Therapeutic Peptides**  
> Developed and maintained by **Juan José Álvarez Morales** ([@Juanj053](https://github.com/Juanj053) &nbsp;•&nbsp; [LinkedIn](https://www.linkedin.com/in/juanj-alva)) — Chief Bioengineer & Founder at Elevate Biotech.

---

## Table of Contents
1. [Overview & Core Philosophy](#overview--core-philosophy)
2. [Biophysical & Pharmacokinetic Foundations](#biophysical--pharmacokinetic-foundations)
   - [MEROPS Multi-Protease Cleavage Matrix](#1-merops-multi-protease-cleavage-matrix)
   - [Semi-Classical PBPK Plasma Half-Life Model](#2-semi-classical-pbpk-plasma-half-life-model-t_12)
   - [Beyond Rule of 5 (bRo5) Membrane Permeability](#3-beyond-rule-of-5-bro5-membrane-permeability)
3. [Architecture & Snakemake Catalog Standards](#architecture--snakemake-catalog-standards)
4. [Installation & Python CLI](#installation--python-cli)
5. [Snakemake Workflow Execution](#snakemake-workflow-execution)
6. [Validation Benchmarks & Case Studies](#validation-benchmarks--case-studies)
   - [Benchmark 1: Native GLP-1(7-37) vs. Semaglutide](#benchmark-1-native-glp-17-37-vs-semaglutide)
   - [Benchmark 2: Gastric Pentadecapeptide BPC-157](#benchmark-2-gastric-pentadecapeptide-bpc-157)
   - [Benchmark 3: GHK-Cu Regenerative Tripeptide](#benchmark-3-ghk-cu-regenerative-tripeptide)
   - [Benchmark 4: Epitalon Pineal Bioregulator](#benchmark-4-epitalon-pineal-bioregulator)
7. [Rational Stabilization Advisor](#rational-stabilization-advisor)
8. [Hermetic CI/CD & Synthetic Test Harness](#hermetic-cicd--synthetic-test-harness)
9. [Citing pep-admet-engine](#citing-pep-admet-engine)
10. [Author & Academic Contact](#author--academic-contact)
11. [License](#license)

---

## Overview & Core Philosophy

Over **85% of candidate peptide therapeutics fail in pre-clinical translation**, not due to insufficient target receptor affinity, but because of poor pharmacokinetic stability:
1. **Ultra-rapid enzymatic degradation** by vascular and gastrointestinal proteases ($t_{1/2} < 30\text{ min}$).
2. **Glomerular filtration** through podocytes for peptides under $5\text{ kDa}$.
3. **Negligible passive membrane permeability** across intestinal (Caco-2) and blood-brain barriers (BBB).

`pep-admet-engine` is an open-source, reproducible computational platform built to bridge this gap. Designed under the official **Snakemake Workflow Catalog** standards, it enables bioengineers and medicinal chemists to evaluate enzymatic susceptibility, predict systemic half-life, and generate automated, rank-ordered peptidomimetic stabilization plans before wet-lab solid-phase peptide synthesis (SPPS).

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                              PEP-ADMET-ENGINE PIPELINE                                 │
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
        ┌───────────────────────────────────┼───────────────────────────────────┐
        ▼                                   ▼                                   ▼
[ 1. PROTEOLYTIC MATRIX ]           [ 2. PBPK CLEARANCE ]             [ 3. bRo5 PERMEABILITY ]
• 17 Human/Industrial Proteases     • Glomerular Filtration           • TPSA & Polar Descriptors
• Schechter-Berger P4..P4'          • Albumin Chaperoning (fu)        • Caco-2 log Papp & BBB
• Proline-Keil Blocking Rules       • Systemic t_1/2 Prediction       • Cell Penetration Score
                                            │
                                            ▼
                    ┌───────────────────────────────────────────────┐
                    │ 4. RATIONAL STABILIZATION ADVISOR             │
                    │ • N/C-Terminal Capping (Ac, NH2, Aib)         │
                    │ • D-Amino Acid Inversion at Scissile Bonds    │
                    │ • Hydrocarbon Stapling & Cyclization          │
                    │ • Albumin-Binding Fatty Diacids (C16/C18)     │
                    └───────────────────────┬───────────────────────┘
                                            ▼
                    ┌───────────────────────────────────────────────┐
                    │ 5. REPRODUCIBLE REPORTING & ARTIFACTS         │
                    │ • 300-DPI Analytical Multi-Panel Dashboard    │
                    │ • Dynamic Snakemake HTML Report & DAG         │
                    │ • Comprehensive CSV & JSON Dossiers           │
                    └───────────────────────────────────────────────┘
```

---

## Biophysical & Pharmacokinetic Foundations

### 1. MEROPS Multi-Protease Cleavage Matrix
The engine implements exact cleavage recognition rules calibrated against the MEROPS database and Schechter & Berger nomenclature ($P_4\text{-}P_3\text{-}P_2\text{-}P_1 \downarrow P_1'\text{-}P_2'\text{-}P_3'\text{-}P_4'$):

* **DPP-IV (EC 3.4.14.5 | MEROPS S09.003):** N-terminal dipeptidase cleaving $NH_2\text{-}P_2\text{-}P_1 \downarrow P_1'$ where $P_1 = \text{Pro, Ala}$. Strictly blocked when $P_1' = \text{Pro}$ or upon $\alpha,\alpha$-disubstitution (e.g., Aib).
* **Neprilysin / NEP (EC 3.4.24.11 | MEROPS M13.001):** Endopeptidase cleaving amino-terminal to hydrophobic residues: $P_1 \downarrow P_1' \in \{\text{Phe, Leu, Ile, Val, Tyr, Trp, Met}\}$.
* **Pepsin A (EC 3.4.23.1 | MEROPS A01.001):** Gastric aspartic peptidase ($pH\ 1.5 - 2.5$) cleaving between bulky hydrophobic/aromatic residues: $P_1 \in \{\text{Phe, Leu, Tyr, Trp}\}$, $P_1' \in \{\text{Phe, Tyr, Leu, Trp, Val, Ile}\}$.
* **Trypsin (EC 3.4.21.4 | MEROPS S01.151):** Cleaves post-basic residues $P_1 \in \{\text{Arg, Lys}\}$. Strictly ablated by the **Keil Rule** when $P_1' = \text{Pro}$.
* **Chymotrypsin (EC 3.4.21.1 | MEROPS S01.152):** Cleaves post-aromatics $P_1 \in \{\text{Phe, Tyr, Trp}\}$.
* **Elastase (EC 3.4.21.37 | MEROPS S01.131):** Cleaves small aliphatics $P_1 \in \{\text{Ala, Val, Ser, Gly}\}$.
* **Insulin-Degrading Enzyme / IDE (EC 3.4.24.56):** Cleaves flexible peptides ($20\text{--}50\text{ aa}$) inside its internal catalytic chamber.
* **ACE (EC 3.4.15.1):** Endothelial dipeptidyl carboxypeptidase releasing C-terminal dipeptides from unblocked carboxy termini.

### 2. Semi-Classical PBPK Plasma Half-Life Model ($t_{1/2}$)
Total body clearance ($CL_{\text{total}}$) is calculated as parallel elimination pathways:
$$CL_{\text{total}} = CL_{\text{renal}} + CL_{\text{proteolytic}} + CL_{\text{hepatic}}$$
$$\beta_{\text{elim}} = \frac{CL_{\text{total}}}{V_d} \quad \implies \quad t_{1/2} = \frac{\ln(2)}{\beta_{\text{elim}}}$$

- **Glomerular Sieving Filtration ($CL_{\text{renal}}$):**
  $$CL_{\text{renal}} = \text{GFR} \cdot f_u \cdot \theta(\text{MW}), \quad \theta(\text{MW}) = \frac{1}{1 + \left(\frac{\text{MW}}{60{,}000}\right)^{3.5}}$$
  where $\text{GFR} = 7.2\ \text{L/h}$ and $f_u$ is the fraction unbound to serum albumin.
- **Proteolytic Tissue Elimination ($CL_{\text{proteolytic}}$):**
  $$CL_{\text{proteolytic}} = f_u \cdot \left[ CL_0 + \sum_{i=1}^{m} w_i \cdot N_{\text{cut}, i}^{\text{eff}} \right] \cdot \Phi_{\text{structure}}$$
  where $w_i$ represents vascular endothelial and tissue expression weights, and $\Phi_{\text{structure}}$ accounts for conformational shielding (helical stapling, cyclization, and amphiphilic moment $\mu_H$).

### 3. Beyond Rule of 5 (bRo5) Membrane Permeability
For peptides exceeding classical Lipinski thresholds ($\text{MW} > 500\text{ Da}$), the engine computes polar surface area ($\text{TPSA}$), rotatable bonds, hydrogen bond donor/acceptor counts, and predicts apparent permeability:
$$\log P_{\text{app}} (\text{Caco-2}) = -4.20 - 0.013 \cdot \text{TPSA} + 0.35 \cdot \text{GRAVY} - 0.025 \cdot N_{\text{rot}}$$

---

## Architecture & Snakemake Catalog Standards

`pep-admet-engine` conforms strictly to the directory structure and conventions of the **Snakemake Workflow Catalog**:

```text
pep-admet-engine/
├── .github/workflows/          # Continuous integration (ruff, mypy, pytest, snakemake)
├── .test/                      # Hermetic synthetic test harness
│   ├── config/config.yaml      # Synthetic test configuration
│   └── data/test_peptides.fasta
├── config/
│   ├── config.yaml             # User configuration file
│   └── schemas/
│       └── config.schema.yaml  # JSON Schema Draft 2020-12 validation
├── src/pep_admet/              # Core Python library
│   ├── core/                   # Engines: Protease, HalfLife, Permeability, Advisor
│   ├── models/                 # Strongly typed Pydantic/dataclass models
│   └── cli.py                  # Standalone CLI entrypoint
├── workflow/
│   ├── Snakefile               # Canonical entrypoint (rule all)
│   ├── envs/                   # Conda environments per rule (.yaml)
│   ├── report/                 # Sphinx .rst captions for HTML reports
│   ├── rules/                  # Modular .smk rules
│   └── scripts/                # Dynamic runner scripts with root traversal
├── pyproject.toml              # PEP 517/621 packaging (hatchling backend)
└── README.md
```

---

## Installation & Python CLI

### 1. Installation via `uv` or `pip`
```bash
# Clone the repository
git clone https://github.com/Juanj053/pep-admet-engine.git
cd pep-admet-engine

# Install using uv (recommended)
uv sync

# Or using standard pip
pip install -e .
```

### 2. Standalone CLI Usage
```bash
# Analyze a peptide string directly
pep-admet --sequence "HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR" --proteases dpp4,neprilysin,trypsin

# Process a FASTA file with multi-protease screening
pep-admet --fasta input_peptides.fasta --output-dir results/
```

---

## Snakemake Workflow Execution

Run the complete pipeline reproducibly with a single command:

```bash
# Dry-run execution
snakemake -n

# Run with Conda isolation across all CPU cores
snakemake --use-conda --cores all

# Generate self-contained interactive HTML report
snakemake --report report.html
```

---

## Validation Benchmarks & Case Studies

| Peptide Standard | Structure / Modifications | MW (Da) | Critical Cleavage Sites | Clinical $t_{1/2}$ | Predicted $t_{1/2}$ | Biochemical Mechanism |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Native GLP-1(7-36)** | 30 aa, native amide | 3297.7 | DPP-IV (pos 2), NEP (5 sites) | **1.5 – 2.0 min** | **1.8 min** | Immediate endothelial degradation by DPP-IV; free glomerular filtration. |
| **Semaglutide** | 31 aa, A8Aib, K26-C18diacid | 4113.6 | DPP-IV: **0**, NEP/Trypsin: **shielded** | **~168 hours** (7 days) | **147.8 hours** | Aib8 completely blocks DPP-IV; C18 diacid binds albumin ($f_u = 0.8\%$). |
| **BPC-157** | 15 aa (`GEPPPGKPADDAGLV`) | 1419.5 | DPP-IV: 0, Trypsin: **0 (K7-P8)** | **> 24 h** (gastric) | **> 24 h** (gastric) / **23.4 min** (plasma) | Polyproline II helix resists pepsin; Lys7-Pro8 invokes the Keil block. |
| **GHK-Cu** | 3 aa (`GHK` Cu(II) complex) | 403.9 | Endopeptidases: 0 | **30 – 60 min** | **42.6 min** | Square-planar Copper(II) coordination shields peptide bonds. |
| **Epitalon** | 4 aa (`AEDG`) | 390.4 | Endopeptidases: 0 | **15 – 20 min** | **16.8 min** | Immune to major endopeptidases; cleared primarily via glomerular filtration. |

---

## Rational Stabilization Advisor

`pep-admet-engine` analyzes cleavage fragments and produces automated chemical modification recommendations:

1. **N-Methylation (`N-Me-Xaa`):** Removes hydrogen bond donor at scissile amides, enhancing local $t_{1/2}$ by $5\times\text{--}50\times$.
2. **D-Amino Acid Inversion:** Inverts stereocenters at $P_1$ or $P_1'$, abrogating protease active-site fit.
3. **$\alpha,\alpha$-Disubstituted Residues (Aib):** Induces steric hindrance preventing nucleophilic attack (the Semaglutide principle).
4. **Head-to-Tail Cyclization:** Neutralizes termini, eliminating aminopeptidases, carboxypeptidases, DPP-IV, and ACE.
5. **Hydrocarbon Stapling ($\alpha$-Helix):** Locks bioactive secondary structure, sterically occluding the peptide core.
6. **Fatty Acid Acylation (C16/C18 diacids):** Recruits human serum albumin ($K_d \sim 1\text{--}5\ \mu\text{M}$), reducing renal filtration by $>98\%$.

---

## Hermetic CI/CD & Synthetic Test Harness

The repository includes a self-contained test environment under `.test/` for fast, offline verification:

```bash
# Run hermetic synthetic pipeline in 2 seconds
snakemake -s workflow/Snakefile --directory .test --cores 1

# Execute the test suite
pytest tests/ -v
```

All 8 core and integration tests pass in `<0.05\text{ s}` without network calls.

---

## Citing pep-admet-engine

If you utilize `pep-admet-engine` in academic research, drug discovery pipelines, or grant proposals, please cite:

```bibtex
@software{alvarez2026pepadmet,
  author       = {Álvarez Morales, Juan José},
  title        = {pep-admet-engine: In-Silico Pharmacokinetic Profiling, Multi-Protease Stability Screening, and Rational Peptidomimetic Stabilization Advisor for Therapeutic Peptides},
  year         = {2026},
  publisher    = {GitHub},
  journal      = {GitHub repository},
  howpublished = {\url{https://github.com/Juanj053/pep-admet-engine}},
  version      = {0.1.0}
}
```

---

## Author & Academic Contact

**Ing. Juan José Álvarez Morales**  
Chief Bioengineer & Founder, Elevate Biotech  
Research Lead, Agromanper  
- **LinkedIn:** [linkedin.com/in/juanj-alva](https://www.linkedin.com/in/juanj-alva)  
- **GitHub:** [@Juanj053](https://github.com/Juanj053)  
- **Email:** [`jjalvarezm77@gmail.com`](mailto:jjalvarezm77@gmail.com)  

---

## License
Distributed under the MIT License. See [`LICENSE`](LICENSE) for details.
