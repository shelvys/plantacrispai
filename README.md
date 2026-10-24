<div align="center">

<img src="public/logo192.png" alt="PlantaCRISPAI" width="110" />

# PlantaCRISPAI

**An AI-assisted CRISPR genotyping platform for Sanger sequencing data from genome-edited plants.**

[![Live](https://img.shields.io/badge/live-plantacrispai.netlify.app-00D4AA?style=flat-square)](https://plantacrispai.netlify.app)
[![License](https://img.shields.io/badge/license-Proprietary-FF4D6D?style=flat-square)](LICENSE)
[![Copyright](https://img.shields.io/badge/copyright-application%20pending-F0C040?style=flat-square)](#license-and-copyright)
[![React](https://img.shields.io/badge/React-19-4D9FFF?style=flat-square&logo=react&logoColor=white)](https://react.dev)
[![ICAR-IARI](https://img.shields.io/badge/ICAR--IARI-Plant%20Physiology-2E7D32?style=flat-square)](https://iari.res.in)

[Live app](https://plantacrispai.netlify.app) · [Quick start](#quick-start) · [How it works](#how-it-works) · [Architecture](docs/ARCHITECTURE.md) · [Scientific methods](docs/SCIENTIFIC_METHODS.md) · [Citing](#citation)

</div>

---

## What this is

Screening a CRISPR plant population by Sanger sequencing is routine, tedious and
easy to get wrong. You get back a folder of `.ab1` traces, and for every plant you
need to answer the same questions: did the edit happen, how big is the indel, is
it in frame, and is this plant wild-type, mosaic, heterozygous or homozygous?

PlantaCRISPAI answers those questions from the raw traces. Drop in the `.ab1`
files, paste the guide RNA and the reference amplicon, pick the generation and
ploidy, and it returns a per-plant genotype call with the evidence behind it —
alignment, indel size and direction, frameshift status, R² either side of the cut
site, an ICE-equivalent score, and a confidence from an ensemble classifier. It
then tells you what that means for your breeding programme at T0, T1, T2 or T3,
designs a CAPS or dCAPS marker so you can screen the rest of the population by
gel instead of by sequencing, and writes the Methods paragraph for your
manuscript.

**Everything runs in the browser.** Sequences are parsed, aligned and classified
client-side. No file is uploaded, no account is needed, and the app works offline
once loaded — which matters when the data is unpublished.

> [!IMPORTANT]
> PlantaCRISPAI is a research tool. It is not a certified diagnostic device and
> is not validated for regulatory submission or seed certification. Confirm
> genotype calls by an independent laboratory method before acting on them.

---

## Why it exists

Built at the **Division of Plant Physiology, ICAR-Indian Agricultural Research
Institute**, under the PIU-GET project *"Enhancing Climate Resilience and
Ensuring Food Security with Genome Editing Tools"* (Code 12-197).

The existing tools are good but awkward for plant work. ICE and TIDE are
web services that want your sequences on someone else's server, assume a diploid
mammalian context, and say nothing about T-generation advancement or somaclonal
variation. CRISPResso targets amplicon NGS rather than the Sanger traces most
plant labs actually generate. PlantaCRISPAI was written to close that gap:
generation-aware, ploidy-aware, plant-first, and private by construction.

---

## Features

| | Capability | Detail |
|---|---|---|
| 🧬 | **Native `.ab1` parsing** | Reads the Applied Biosystems ABIF binary directly — base calls, per-base Phred scores, all four fluorescence channels, peak locations, sample name and run date. No conversion step. |
| 🤖 | **Ensemble zygosity calling** | Rule-based + ICE-equivalent + Random Forest over 45 features, with generation-aware thresholds. Returns wild-type, mosaic, heterozygous or homozygous with a confidence score. |
| 📐 | **Gap-aware alignment** | Smith-Waterman local alignment anchored to the reference, with automatic reverse-primer detection and strand correction. |
| ✂️ | **Multi-gRNA / dual-cut** | Analyses each guide independently and detects inter-cut deletions, with per-guide ICE, R² and zygosity. |
| 🌾 | **T0–T3 biology** | Expected Mendelian segregation ratios, T-DNA screening advice, somaclonal variation flags and homozygous-advancement recommendations per generation. |
| 🔬 | **Chromatogram viewer** | Canvas traces with gRNA and PAM highlighting, cut-site marker, per-base quality bars, wheel zoom, and an edited-vs-wild-type dual overlay. |
| 🧪 | **CAPS / dCAPS designer** | Scans 67 restriction enzymes for sites the indel disrupts, scores reliability, simulates the gel, and designs dCAPS primers when no natural marker exists. |
| 🎯 | **gRNA designer** | Scores candidate 20 bp + NGG targets on 28 features — position weights, nearest-neighbour thermodynamics, secondary structure, plant-specific corrections and CAPS compatibility. |
| 📄 | **Publication output** | Allele plots, indel distribution charts, per-sample chromatograms, CSV export and a pre-written Methods section. |
| 🔒 | **Private by construction** | No server, no upload, no telemetry, no signup. |

---

## Quick start

### Use it

Go to **[plantacrispai.netlify.app](https://plantacrispai.netlify.app)** and press
*Launch App*. A built-in *OsAAP3* rice T1 dataset is included under **Quick Demo**
if you want to see the pipeline before committing your own data.

### Run it locally

Requires **Node.js ≥ 18** (20 LTS recommended).

```bash
git clone https://github.com/shelvys/plantacrispai.git
cd plantacrispai
npm install
npm start           # http://localhost:3000
```

```bash
npm run build       # optimised production bundle in build/
npm run lint        # ESLint over src/
```

> Local use is permitted for non-commercial academic research under the terms in
> [`LICENSE`](LICENSE). It does not grant redistribution or modification rights.

---

## How it works

```
 .ab1 / FASTA / text          reference amplicon + gRNA(s)
         │                              │
         ▼                              ▼
 ┌─────────────────┐          ┌──────────────────────┐
 │  ABIF parser    │          │  cut-site derivation │
 │  traces·Q·calls │          │  strand-aware, N×gRNA│
 └────────┬────────┘          └──────────┬───────────┘
          │  quality trimming            │
          └──────────────┬───────────────┘
                         ▼
          ┌──────────────────────────────┐
          │ Smith-Waterman alignment     │
          │ gap-aware · auto-orienting   │
          └──────────────┬───────────────┘
                         ▼
          ┌──────────────────────────────┐
          │ indel size · direction       │
          │ frameshift · R²(up/down)     │
          │ downstream mismatch rate     │
          └──────────────┬───────────────┘
                         ▼
          ┌──────────────────────────────┐
          │ ICE-equivalent score         │
          │ 45 features → RF + rules     │
          │ generation-aware thresholds  │
          └──────────────┬───────────────┘
                         ▼
     zygosity call · confidence · flags · interpretation
                         │
          ┌──────────────┴───────────────┐
          ▼                              ▼
   CAPS / dCAPS marker            report · CSV · Methods
   + simulated gel
```

1. **Configure** — guide RNA sequence(s), reference amplicon, organism,
   generation (T0–T3) and ploidy. Multi-guide experiments supported.
2. **Upload** — drag in `.ab1` files, or paste sequences. Choose a trimming mode;
   the app warns if trimming would encroach on the cut site.
3. **Analyse** — alignment, indel detection, ICE scoring and ensemble
   classification run locally.
4. **Explore and publish** — dashboard, chromatograms, advancement strategy,
   CAPS design, report download.

Full algorithmic detail is in **[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md)**;
the published methods behind each score are cited in
**[docs/SCIENTIFIC_METHODS.md](docs/SCIENTIFIC_METHODS.md)**.

---

## Validation

The classifier was developed against an *OsAAP3* rice T1 panel genotyped
independently with **ICE** and **TIDE**, covering wild-type, mosaic,
heterozygous and homozygous classes.

---

## Project layout

```
plantacrispai/
├── public/                     static assets and HTML template
├── src/
│   ├── engine/                 science — no React, independently testable
│   │   ├── ab1Parser.js          ABIF binary reader + quality trimming
│   │   ├── analysisEngine.js     alignment, ICE, zygosity, reporting (v5.0)
│   │   ├── capsEngine.js         enzyme DB, CAPS/dCAPS design, gel bands
│   │   ├── grnaDesigner.js       candidate discovery and ranking
│   │   ├── grnaFeatures.js       position-specific feature extraction
│   │   ├── thermodynamics.js     nearest-neighbour ΔG, Tm, folding
│   │   └── demoData.js           OsAAP3 T1 demonstration panel
│   ├── components/             presentation
│   │   ├── AB1UploadPanel.js     upload + four trimming modes
│   │   ├── ChromatogramViewer.js canvas trace renderer
│   │   ├── CAPSDesignerPanel.js  CAPS workflow (3 modes)
│   │   ├── GRNACAPSDesignerPanel.js  pre-editing guide + marker design
│   │   ├── GelSimulator.js       agarose gel simulation
│   │   ├── MLQualityPanel.js     per-sample confidence and QC
│   │   └── Navbar.js
│   ├── pages/
│   │   ├── LandingPage.js        marketing, team, citation
│   │   └── AnalysisApp.js        eight-tab analysis workspace
│   ├── App.js  index.js  index.css
├── docs/                       architecture and scientific methods
├── .github/                    CI, issue and PR templates
├── LICENSE  NOTICE  THIRD_PARTY_LICENSES.md  CITATION.cff
└── netlify.toml
```

The `engine/` modules are plain JavaScript with no React dependency. That
separation is deliberate: the science can be unit-tested, ported or called from a
CLI without touching the UI.

---

## Technology

Deliberately minimal. The only runtime dependencies are **React 19** and
**React DOM** — every algorithm, chart, chromatogram renderer and gel simulation
is implemented from scratch in this repository. There is no charting library, no
bioinformatics package, no component framework and no backend.

| | |
|---|---|
| Runtime | React 19, React DOM |
| Build | Create React App (`react-scripts` 5.0.1), webpack, Babel |
| Rendering | HTML5 Canvas (chromatograms, gels), inline styles |
| Deployment | Netlify — static hosting, SPA redirect, immutable asset caching |

`react-scripts` is deprecated upstream. It is retained for now so this source
tree reproduces the deployed artifact exactly; a Vite migration is planned for
1.1.0.

---

## License and copyright

Copyright © 2026 **ICAR-Indian Agricultural Research Institute, Division of Plant
Physiology, New Delhi** and **Shelvy S.** All rights reserved.

A **copyright application covering this software is currently under processing**
through ICAR-IARI. The source is published here for transparency, citation and
academic evaluation — not as open source.

**You may** read the source, run it unmodified for non-commercial academic
research or teaching, and evaluate it. **You may not** redistribute it, modify
it, create derivative works, or use it commercially without prior written
permission. Full terms in **[`LICENSE`](LICENSE)**.

> [!WARNING]
> Because this repository is public, GitHub's Terms of Service grant other
> GitHub users a limited right to view and fork it within GitHub, regardless of
> the license above. If the institute requires that no copy exist outside its
> control while the application is processing, make the repository private —
> Settings → General → Danger Zone → Change visibility.

For commercial licensing, redistribution rights or collaboration, contact
<crispricar@gmail.com>.

Third-party components retain their own licenses — see
[`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md) and [`NOTICE`](NOTICE).

---

## Citation

If PlantaCRISPAI contributed to your research, please cite it. GitHub reads
[`CITATION.cff`](CITATION.cff), so the *Cite this repository* button on the
sidebar will generate these for you.

**APA (7th edition)**

> Shelvy, S., Archana, S., Ramesh, R., Sathee, L., & Viswanathan, C. (2026).
> *PlantaCRISPAI: An AI-powered universal CRISPR genotyping platform for Sanger
> sequencing data from genome-edited plants* [Software]. Division of Plant
> Physiology, ICAR-Indian Agricultural Research Institute.
> https://plantacrispai.netlify.app

**BibTeX**

```bibtex
@software{shelvy2026plantacrispai,
  author    = {Shelvy, S. and Archana, S. and Ramesh, R. and Sathee, L. and Viswanathan, C.},
  title     = {{PlantaCRISPAI}: An {AI}-powered universal {CRISPR} genotyping platform for {Sanger} sequencing data from genome-edited plants},
  year      = {2026},
  publisher = {Division of Plant Physiology, ICAR-Indian Agricultural Research Institute},
  url       = {https://plantacrispai.netlify.app}
}
```

*Preprint — citation to be updated following peer review.*

---

## Team

**Development** — [Shelvy S.](AUTHORS.md) (lead), Archana S., Ramesh R.
**Supervision** — Dr. Lekshmy Sathee, Dr. Viswanathan C.

Full roles and contributions in [`AUTHORS.md`](AUTHORS.md).

---

## Contact

**PIU-GET Project (Code 12-197)**
Division of Plant Physiology, ICAR-IARI
Pusa, New Delhi – 110012, India

📧 <crispricar@gmail.com> · 🌐 [genetools.iari.res.in](https://genetools.iari.res.in) · 🏛 [ICAR](https://icar.org.in/) · [IARI](https://iari.res.in/en/index.php)

Bug reports and feature requests: [open an issue](https://github.com/shelvys/plantacrispai/issues).
Security concerns: see [`SECURITY.md`](SECURITY.md).

---

<div align="center">
<sub>Built for plant scientists, by plant scientists · ICAR-IARI, New Delhi</sub>
</div>
