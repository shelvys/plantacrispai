# Contributing to PlantaCRISPAI

Thank you for your interest. Please read this first — the licensing situation
here is unusual.

## Licensing — read before you write code

PlantaCRISPAI is **not open source**. It is published under a proprietary
license ([`LICENSE`](LICENSE)) while a copyright application filed through
ICAR-IARI is processing.

This means:

- **Bug reports, feature requests, scientific feedback and validation data are
  very welcome.** These are the most useful contributions right now.
- **Unsolicited pull requests cannot be merged.** Accepting external code during
  a pending copyright application complicates authorship and chain of title. If
  you have a change in mind, open an issue first so we can discuss whether and
  how it can be accepted.
- If a code contribution is agreed, you will be asked to sign a contributor
  agreement assigning copyright in the contribution to ICAR-IARI before it can be
  merged.

If you want to build on this work, email <crispricar@gmail.com>. Collaboration is
genuinely welcome; it just needs to go through the institute.

## Reporting a bug

Open an [issue](https://github.com/shelvys/plantacrispai/issues) using the bug
template. The most useful report includes:

- what you did, what you expected, what happened
- browser and operating system
- experiment context — organism, generation, ploidy, guide RNA length, single or
  multi-guide
- any browser console output (F12 → Console)

> [!CAUTION]
> **Do not attach `.ab1` files or sequences you cannot publish.** Issues are
> public. If the bug depends on your data, say so and we will arrange a private
> channel. Where possible, reproduce with the built-in *OsAAP3* demo dataset and
> report that instead.

## Reporting a scientific problem

Misclassifications matter more than crashes. If PlantaCRISPAI calls a genotype
that disagrees with your ICE, TIDE or cloning result, please tell us. Include:

- the call PlantaCRISPAI made and the call you believe is correct
- what established the correct answer — ICE, TIDE, cloning and sequencing,
  progeny segregation
- the reported indel size, R² upstream and downstream, and ICE-equivalent score
- generation and ploidy

These reports drive classifier improvement more than anything else.

## Development setup

Requires Node.js ≥ 18 (20 LTS recommended).

```bash
git clone https://github.com/shelvys/plantacrispai.git
cd plantacrispai
npm install
npm start
```

| Command | Purpose |
|---|---|
| `npm start` | dev server with hot reload, http://localhost:3000 |
| `npm run build` | production bundle in `build/` |
| `npm run lint` | ESLint over `src/` |
| `npm test` | test runner (no suite yet) |

## Code conventions

If a change is agreed and you are writing code, follow what is already there.

**Keep the science out of the UI.** Anything in `src/engine/` must stay free of
React imports — pure functions, data in and data out. This is what makes the
pipeline testable and portable, and it is the one architectural rule that should
not bend.

**Document algorithms in the file header.** Every engine module opens with a
block stating purpose, method, the published basis, and — importantly — what
changed and why in each version. Those headers record real debugging history
(the pixel-vs-base indexing bug in the chromatogram viewer, the homozygous
threshold correction in the classifier). Keep writing them.

**Cite your source.** A new scoring term, weight or threshold needs a citation or
an explicit note that it is empirical. Add it to
[`docs/SCIENTIFIC_METHODS.md`](docs/SCIENTIFIC_METHODS.md).

**Do not inflate claims.** If a number is heuristic, say so. If a model was
tuned on the same data it is evaluated on, say so. The documentation flags these
honestly and should continue to.

**Style:** 2-space indent, semicolons, single quotes, named exports from engine
modules. ESLint (`react-app` config) must pass clean. Run `npm run lint` before
proposing anything.

**Privacy is architectural.** Nothing may transmit user sequence data off the
client. No analytics, no error reporting services, no remote API calls on
sequence data. A change that breaks this will be rejected regardless of merit.

## Commit messages

Conventional Commits:

```
feat(engine): add tetraploid-aware zygosity thresholds
fix(ab1Parser): honour FWO_.1 channel order on 3730xl traces
docs(methods): reconcile reported validation accuracy
```

Types: `feat`, `fix`, `docs`, `refactor`, `test`, `chore`, `perf`.

## What would help most

In rough order of value:

1. **Validation data** — `.ab1` panels from other species, genes, generations or
   ploidy levels with independently confirmed genotypes. The classifier's biggest
   limitation is that it has only been tested on rice T1 *OsAAP3*.
2. **A test suite.** The engine modules are pure functions and directly
   testable; there is no suite yet.
3. **Reconciling the reported accuracy figures** — three inconsistent numbers
   are recorded in the codebase; see
   [`docs/SCIENTIFIC_METHODS.md`](docs/SCIENTIFIC_METHODS.md).
4. **Web Worker offloading** for alignment, so large batches do not block the UI.
5. **Documentation and usability feedback** from people screening real
   populations.

## Code of conduct

Participation is governed by [`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md).

## Questions

<crispricar@gmail.com> — Division of Plant Physiology, ICAR-IARI, New Delhi.
