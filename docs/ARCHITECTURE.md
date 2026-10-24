# Architecture

How PlantaCRISPAI is put together, and why. For the published science behind each
score, see [SCIENTIFIC_METHODS.md](SCIENTIFIC_METHODS.md).

## Design principles

**The science is separable from the UI.** Everything in `src/engine/` is plain
JavaScript with no React import. Each module takes data and returns data. This
means the pipeline can be unit-tested, called from Node, or ported to a CLI or
server without rewriting it, and it keeps numerical logic out of render paths.

**No computation leaves the browser.** There is no API, no backend and no
telemetry. Unpublished sequence data is handled client-side from parse to report.
This constrains the design — everything must run acceptably fast in a single tab
on a laptop — and rules out server-side tools like BLAST for off-target search.

**No heavy dependencies.** React and React DOM are the only runtime packages.
Alignment, thermodynamics, feature extraction, the classifier, chromatogram
rendering and gel simulation are all implemented here. The bundle is ~187 kB
gzipped.

**Faithful to the wet lab.** Thresholds, flags and interpretations are expressed
in terms a plant breeder uses — generation, ploidy, segregation ratio,
advancement decision — rather than in abstract statistics.

## Module map

```
src/
├── engine/          pure functions · no React
│   ├── ab1Parser.js
│   ├── analysisEngine.js
│   ├── capsEngine.js
│   ├── grnaDesigner.js
│   ├── grnaFeatures.js
│   ├── thermodynamics.js
│   └── demoData.js
├── components/      presentation · props in, callbacks out
└── pages/           composition and application state
```

Dependency direction is strictly one-way: `pages → components → engine`. No
engine module imports from `components/` or `pages/`.

Within the engine there is one internal dependency chain:

```
thermodynamics.js ──┐
                    ├──> grnaFeatures.js ──> grnaDesigner.js
capsEngine.js ──────┘                            │
      └──────────────────────────────────────────┘
analysisEngine.js ──> (consumed by capsEngine for CAPS prediction)
ab1Parser.js ───────> (feeds analysisEngine)
```

## State management

Application state is plain React `useState`, held at two levels:

- `App.js` holds only the active view (`landing` | `app`). Navigation is a state
  toggle, not a router — the app is a single page with no deep-linkable routes,
  so React Router would add a dependency for no benefit.
- `pages/AnalysisApp.js` owns the analysis session: configuration, loaded
  samples, results, active tab and selection. Children receive props and raise
  callbacks.

There is no global store. The state graph is shallow enough that Context or Redux
would add indirection without reducing coupling.

## Pipeline

### 1 · Parsing — `ab1Parser.js`

ABIF is Applied Biosystems' binary container, used by essentially every Sanger
facility. The parser validates the `ABIF` magic bytes, walks the directory of
tagged entries, and extracts:

| Tag | Content |
|-----|---------|
| `PBAS.2` | base calls — the sequence |
| `PCON.2` | per-base Phred quality |
| `DATA.9–12` | raw fluorescence, four channels |
| `PLOC.2` | peak locations in the trace |
| `FWO_.1` | filter wheel order — channel → base mapping |
| `SMPL.1` | sample name |
| `RUND.1` | run date |

`FWO_.1` matters: channel order is not fixed across instruments, so mapping
channels to bases by position rather than by this tag produces silently wrong
traces.

When only a text sequence is supplied, `generateSyntheticTraces()` synthesises
plausible traces so the chromatogram viewer stays usable, clearly derived rather
than measured.

**Trimming** supports four modes — Standard (Phred ≥ 20 over a 10 bp window),
Quality (user-set threshold and window), Manual (explicit per-sample boundaries)
and None. `applyTrim()` accepts the cut-site position and warns when a trim
boundary approaches it, because over-aggressive trimming that clips the cut site
is a common source of false wild-type calls.

### 2 · Cut-site localisation — `analysisEngine.js`

`deriveCutPosition()` locates the guide in the reference on both strands and
places the blunt Cas9 cut 3 bp 5′ of the PAM, respecting orientation.
`deriveAllCutPositions()` extends this to multiple guides, and
`analyseMultiGrna()` reports per-guide metrics plus inter-cut deletions for
dual-cut designs.

### 3 · Alignment

`alignPairwise()` is a Smith-Waterman local alignment returning gapped strings
for query and reference. Before aligning, `detectReadStrand()` and
`autoOrientQuery()` reverse-complement reverse-primer reads — without this,
reverse reads align poorly and are miscalled as heavily edited.

Cut-site coordinates are then mapped from reference space into alignment space by
`findCutInAlignment()`, so indels upstream do not shift the analysis window.

### 4 · Metrics

Around the alignment-space cut position:

- `windowedIndelSize()` — net indel within ±30 bp, with direction and frameshift
  status.
- `r2Segment()` — gap-aware R² computed separately upstream and downstream. A
  clean single-allele edit keeps high R² on both sides; a heterozygote shows the
  downstream trace degrading as two alleles diverge.
- `downstreamMismatchRate()` / `upstreamMismatchRate()` — mismatch fraction either
  side.
- `buildDiscordanceProfile()` — where disagreement begins, used for flagging.

### 5 · Scoring and classification

`computeICE()` produces an ICE-equivalent editing score from downstream mismatch
rate, indel size and downstream R².

`classifyZygosity()` applies generation-aware thresholds. The decisive signal for
a homozygous call is a *clean* post-cut trace — a single allele shifted by the
indel — expressed as `indelSize ≥ 1 && r2dn ≥ 0.97 && downMM < 0.03`. An earlier
version required `indelSize ≥ 3`, which systematically misclassified small
homozygous indels; this is documented in the file header.

`computeMLScore()` combines the rule-based call with Random Forest output over 45
extracted features, yielding a confidence. `flagSample()` raises tiered,
generation-aware warnings — low quality, ambiguous zygosity, cut site near a trim
boundary, discordance between scores.

`getBiologicalInterpretation()` converts the call into breeding guidance per
generation: chimerism expectations and T-DNA screening at T0, Mendelian
segregation at T1, homozygosity confirmation and advancement at T2–T3.

### 6 · CAPS and dCAPS — `capsEngine.js`

A curated database of **67 restriction enzymes** carries recognition sequence,
cut offset, overhang type, incubation temperature, a reliability score and
suppliers. `scanAllRESites()` finds sites (IUPAC-aware via `iupacToRegex()`);
`checkDisruption()` tests whether the detected indel destroys or creates a site;
`calculateCAPSScore()` scores practical reliability; `calcBandSizes()` predicts
gel bands, rendered by `GelSimulator`.

When no natural CAPS marker exists, `designDCAPS()` finds near-matches within two
mismatches and designs a primer carrying a deliberate mismatch that completes a
restriction site in one allele only.

### 7 · gRNA design — `grnaFeatures.js`, `grnaDesigner.js`

`findAllGRNACandidates()` scans both strands for 20 bp protospacers followed by
NGG. `extractGRNAFeatures()` builds a feature vector spanning global composition,
position-specific nucleotide identity, dinucleotide context, thermodynamics from
`thermodynamics.js`, plant-specific corrections, structural penalties
(homopolymers, repeats), PAM context and CAPS compatibility.

`predictGRNAEfficiency()` combines sub-models into a single score;
`rankGRNACandidates()` applies `findParetoOptimal()` so guides that are strong on
efficiency *and* CAPS compatibility surface together rather than being collapsed
into one weighted number.

`assessOffTargetRisk()` is **sequence-heuristic only** — seed-region and
composition based. There is no genome-wide search, because no genome index can be
shipped to a browser. Candidates must be checked against a genome with
CRISPOR, Cas-OFFinder or equivalent before ordering.

### 8 · Thermodynamics — `thermodynamics.js`

Nearest-neighbour model with RNA:DNA parameters (Sugimoto 1995) for guide:target
duplex stability and DNA:DNA parameters (SantaLucia 1998) for primer Tm.
Computes duplex ΔG, seed-region ΔG, self-folding penalty and salt-adjusted Tm at
a plant transformation temperature of 22 °C.

## Rendering

`ChromatogramViewer` draws to HTML5 Canvas rather than SVG or a charting library:
traces run to thousands of points across four channels, and Canvas redraws under
zoom and pan without the DOM node count SVG would require.

Two details are worth noting, both recorded in the file header as past bugs:
trace data is **pixel-indexed, not base-indexed**, so peak scanning must use
pixel coordinates; and the guide is located in the *reference* then mapped to
sample coordinates through the alignment, so highlighting stays correct in the
presence of indels.

## Build and deployment

Create React App 5 (webpack + Babel) produces a content-hashed bundle in
`build/`. `netlify.toml` pins the build command, Node 20, an SPA catch-all
redirect, security headers and immutable caching for `/static/*`.

The deployed CSS hash `main.5c4935f7.css` is reproduced byte-identically by
`npm run build` from this source tree.

## Known constraints

- **Single-threaded.** Alignment runs on the main thread; very large batches
  will block the UI. A Web Worker pool is the natural fix.
- **No persistence.** Reloading loses the session. Deliberate — nothing is
  written to disk — but IndexedDB with explicit opt-in would be reasonable.
- **Heuristic off-target scoring.** See above.
- **No test suite.** The engine's pure-function design makes it directly
  testable; this is the highest-value next contribution.
- **CRA is deprecated.** Vite migration planned for 1.1.0.
