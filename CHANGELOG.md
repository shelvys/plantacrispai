# Changelog

All notable changes to PlantaCRISPAI are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] — 2026-10-25

First public source release, corresponding to the deployment live at
<https://plantacrispai.netlify.app>.

### Added — analysis engine (v5.0)

- Universal genotyping pipeline accepting `.ab1`, FASTA, multi-FASTA and plain
  text, across generations T0–T3 and haploid through hexaploid ploidy.
- Gap-aware Smith-Waterman local alignment with reference anchoring.
- Automatic read-strand detection and reverse-complement correction for
  reverse-primer reads.
- Guide-RNA-driven cut-site localisation, strand aware, with support for
  multiple guides and inter-cut deletion detection.
- Segmented R² upstream and downstream of the cut site, plus downstream
  mismatch rate.
- ICE-equivalent scoring and a generation-aware zygosity classifier
  (wild-type / mosaic / heterozygous / homozygous).
- 45-feature extraction feeding a rule-based + Random Forest ensemble.
- Tiered, generation-aware flagging and per-generation biological
  interpretation, including Mendelian segregation expectations and
  advancement recommendations.
- Dashboard statistics, AI-generated insights, report text and CSV export.

### Added — Sanger chromatogram handling

- Native ABIF (`.ab1`) binary parser reading `PBAS.2` base calls, `PCON.2`
  quality scores, `DATA.9–12` fluorescence channels, `PLOC.2` peak locations,
  `FWO_.1` channel-to-base mapping, sample name and run date.
- Four trimming modes — Standard (Phred ≥ 20, 10 bp window), Quality (user
  threshold), Manual (per-sample sliders) and No Trim — with cut-site proximity
  warnings.
- Canvas chromatogram viewer with Bezier-smoothed traces, gRNA and PAM
  highlighting, cut-site marker, per-base quality bars, wheel zoom, hover
  tooltips, raw/normalised modes, and a dual-trace edited-vs-wild-type overlay.

### Added — marker and guide design

- Restriction enzyme database of 67 entries with recognition sequence, cut
  position, overhang, incubation temperature, reliability score and supplier.
- CAPS predictor that tests whether a detected indel disrupts a restriction
  site, with an ML reliability score and gel band size calculation.
- dCAPS primer designer for targets with no natural CAPS marker.
- Agarose gel simulator rendering expected band patterns per sample.
- Guide RNA designer scanning both strands for 20 bp + NGG targets, ranked by a
  28-feature ensemble covering position-specific weights, nearest-neighbour
  thermodynamics, secondary structure, homopolymer and repeat penalties,
  plant-specific corrections and CAPS compatibility.
- Nearest-neighbour thermodynamics module for gRNA:DNA duplex ΔG, seed-region
  stability, self-folding energy and accurate melting temperature.

### Added — application

- Landing page covering features, workflow, institutional context, team and
  citation, with copyable APA and BibTeX records.
- Eight-tab analysis workspace: Setup, Dashboard, Explorer, AI Insights,
  ML Model, Report, Chromatogram and CAPS Designer.
- Built-in *OsAAP3* rice T1 demonstration dataset.
- Fully client-side execution — no sequence data is transmitted to any server.

### Repository metadata added in this release

These changes affect packaging only; no application logic was altered.

- `public/manifest.json` replaced the unmodified Create React App placeholder
  (`"Create React App Sample"`) with correct application name, description and
  brand theme colours, and the manifest is now linked from `index.html`.
- `apple-touch-icon` link added.
- `netlify.toml` added to pin the build command, Node version, SPA redirect and
  security/caching headers that were previously configured only in the Netlify
  dashboard.
- Licensing, citation, contribution, security and architecture documentation
  added.
- Continuous integration workflow added.

### Known issues

- `public/logo.png` is a 2.2 MB PNG used as both the favicon and the hero
  image. It should be downscaled and a separate small favicon generated.
- Reported classifier accuracy is inconsistent across the codebase: the landing
  page states 89.5 % on 23 rice T1 samples, `engine/analysisEngine.js` states
  18/19 (94.7 %), and `engine/demoData.js` states 17/20. These figures must be
  reconciled against the validation table before any manuscript submission. See
  [`docs/SCIENTIFIC_METHODS.md`](docs/SCIENTIFIC_METHODS.md).
- `react-scripts` 5.0.1 (Create React App) is deprecated upstream. A migration
  to Vite is planned; the toolchain is retained for now so this source tree
  reproduces the deployed artifact exactly.
- The contact form on the landing page is presentational and does not deliver
  messages. Use <crispricar@gmail.com>.
- No automated test suite yet.

[1.0.0]: https://github.com/shelvys/plantacrispai/releases/tag/v1.0.0
