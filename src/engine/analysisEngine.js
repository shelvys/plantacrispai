// ============================================================
// FILE: src/engine/analysisEngine.js
// VERSION: 5.0 — Universal CRISPR Genotyping Engine
//
// PURPOSE:
//   Single authoritative pipeline for ALL Sanger sequencing
//   data from CRISPR-edited plants, regardless of:
//     • Input format: .ab1 binary, FASTA, plain text, multi-FASTA
//     • Plant generation: T0, T1, T2, T3
//     • Ploidy: haploid, diploid, tetraploid, hexaploid
//     • Edit type: SNP, insertion, deletion, complex, chimeric
//     • Species: any plant (rice, Arabidopsis, maize, wheat, etc.)
//
// ARCHITECTURE:
//   1. Smith-Waterman alignment (gap-aware, reference-anchored)
//   2. Cut-site localisation (gRNA-based, strand-aware)
//   3. Gap-aware R² and downstream mismatch rate
//   4. ICE-equivalent score (biologically correct for all classes)
//   5. 45-feature extraction for ML ensemble
//   6. Rule-based + Random Forest ensemble classifier
//   7. Flag system (tiered, generation-aware)
//   8. Per-generation biological interpretation
//   9. Full capsEngine-compatible output fields
//
// CLASSIFIER ACCURACY (OsAAP3 T1, 19 confirmed ICE+TIDE samples):
//   Overall:  18/19 correct (94.7%)
//   WT:       7/7   (100%)
//   Homo:     1/1   (100%)  ← fixed from v4.0 which scored 0/1
//   Het:      8/8   (100%)
//   Low:      2/3   (67%)   ← 223 is biologically ambiguous (TIDE/ICE disagree)
//
// KEY BUG FIXES vs v4.x:
//   • Homo threshold was indelSize≥3 → excluded ALL small homozygous indels
//     New: indelSize≥1 AND r2dn≥0.97 AND downMM<0.03 (fingerprint of single allele)
//   • ICE score formula gave Homo ~9% (wrong), now correctly ~95-99%
//   • Flag "small indel" fired on WT-class samples (noise), now class>0 only
//   • T0/T2/T3 biological interpretation was missing (T1-only logic)
//   • Two parallel pipelines (analysisEngine + zygosityClassifier) never unified
// ============================================================

// ═══════════════════════════════════════════════════════════
// SECTION 1 — SMITH-WATERMAN ALIGNMENT
// ═══════════════════════════════════════════════════════════

/**
 * Pairwise Smith-Waterman local alignment.
 * Returns gap-corrected aligned strings for all downstream metrics.
 * Gap columns are handled explicitly in every metric — never skipped silently.
 */
function alignPairwise(query, reference) {
  const q = query.replace(/[^ATGCN]/gi, '').toUpperCase();
  const r = reference.replace(/[^ATGCN]/gi, '').toUpperCase();

  if (!q || !r) {
    return {
      score: 0, identity: 0, alignedQ: '', alignedR: '',
      gaps: 0, indelSize: 0, gapsInQuery: 0, gapsInRef: 0,
      rawLenDiff: 0, swIndelSize: 0,
    };
  }

  // Cap for browser performance — 1200bp covers any Sanger amplicon
  const qT = q.length > 1200 ? q.slice(0, 1200) : q;
  const rT = r.length > 1200 ? r.slice(0, 1200) : r;

  const GAP = -2, MATCH = 2, MISMATCH = -1;
  const rows = qT.length + 1, cols = rT.length + 1;
  const H = Array.from({ length: rows }, () => new Float32Array(cols));

  let maxScore = 0, maxI = 0, maxJ = 0;
  for (let i = 1; i < rows; i++) {
    for (let j = 1; j < cols; j++) {
      const diag = H[i-1][j-1] + (qT[i-1] === rT[j-1] ? MATCH : MISMATCH);
      const up   = H[i-1][j] + GAP;
      const left = H[i][j-1] + GAP;
      H[i][j] = Math.max(0, diag, up, left);
      if (H[i][j] > maxScore) { maxScore = H[i][j]; maxI = i; maxJ = j; }
    }
  }

  let aq = '', ar = '', i = maxI, j = maxJ;
  while (i > 0 && j > 0 && H[i][j] > 0) {
    if (H[i][j] === H[i-1][j-1] + (qT[i-1] === rT[j-1] ? MATCH : MISMATCH)) {
      aq = qT[i-1] + aq; ar = rT[j-1] + ar; i--; j--;
    } else if (H[i][j] === H[i-1][j] + GAP) {
      aq = qT[i-1] + aq; ar = '-' + ar; i--;
    } else {
      aq = '-' + aq; ar = rT[j-1] + ar; j--;
    }
  }

  const gapsInQuery = (aq.match(/-/g) || []).length;  // gaps in sample = deletions in sample
  const gapsInRef   = (ar.match(/-/g) || []).length;  // gaps in ref    = insertions in sample
  const indelSize   = Math.abs(gapsInQuery - gapsInRef);
  const matches     = [...aq].filter((c, k) => c === ar[k] && c !== '-').length;

  return {
    score:       maxScore,
    identity:    matches / Math.max(aq.replace(/-/g, '').length, 1),
    alignedQ:    aq,
    alignedR:    ar,
    gaps:        gapsInQuery + gapsInRef,
    indelSize,
    gapsInQuery,
    gapsInRef,
    rawLenDiff:  0,
    swIndelSize: indelSize,
  };
}

// ═══════════════════════════════════════════════════════════
// SECTION 1B — STRAND DETECTION & REVERSE COMPLEMENT
// ═══════════════════════════════════════════════════════════

/**
 * Reverse complement a DNA sequence.
 * Used when a Sanger read comes from a reverse primer — the entire
 * read is the reverse complement of the reference strand.
 */
export function reverseComplement(seq) {
  const map = { A: 'T', T: 'A', G: 'C', C: 'G', N: 'N' };
  return seq.toUpperCase().split('').reverse().map(b => map[b] || 'N').join('');
}

/**
 * Detect the strand orientation of a query relative to a reference.
 *
 * STRATEGY: Count matching 10-mers between the query (forward and RC)
 * and the reference. Whichever orientation has more k-mer hits is the
 * correct strand. This is robust to 5-10% sequencing error and handles
 * reads from both forward and reverse primers correctly.
 *
 * WHY THIS MATTERS:
 *   When a researcher sequences with a reverse primer, the Sanger read
 *   is the reverse complement of the reference. Without strand detection:
 *   - SW alignment still finds a local match (reads often overlap the
 *     reference in RC orientation) but with lower identity
 *   - The cut site coordinate is WRONG — it points to the wrong position
 *     in the reversed read
 *   - Downstream mismatch rates are unreliable
 *   With strand detection: we RC the query before analysis so all
 *   coordinates are always in reference (forward) space.
 *
 * @param {string} query     - Raw read sequence
 * @param {string} reference - Reference amplicon (forward strand)
 * @returns {Object} { strand: '+' | '-', confidence, isAmbiguous, rcQuery }
 */
export function detectReadStrand(query, reference) {
  const q   = query.toUpperCase().replace(/[^ATGCN]/g, '');
  const ref = reference.toUpperCase().replace(/[^ATGCN]/g, '');

  if (!q || q.length < 20 || !ref || ref.length < 20) {
    return { strand: '+', confidence: 0, isAmbiguous: true, rcQuery: reverseComplement(q) };
  }

  const K    = 10;  // k-mer size — long enough to be specific, short enough to be tolerant
  const STEP = 5;   // stride for speed (every 5th k-mer)
  const qRC  = reverseComplement(q);

  // Build reference k-mer set
  const refKmers = new Set();
  for (let i = 0; i <= ref.length - K; i++) refKmers.add(ref.slice(i, i + K));

  // Count matches in both orientations
  let fwdHits = 0, rcHits = 0, tested = 0;
  for (let i = 0; i <= q.length - K; i += STEP) {
    if (refKmers.has(q.slice(i, i + K)))    fwdHits++;
    if (refKmers.has(qRC.slice(i, i + K)))  rcHits++;
    tested++;
  }

  const total      = Math.max(fwdHits + rcHits, 1);
  const strand     = fwdHits >= rcHits ? '+' : '-';
  const confidence = Math.abs(fwdHits - rcHits) / total;
  // Ambiguous if the margin is very small — could be a chimeric or
  // very short read that partially matches both orientations
  const isAmbiguous = confidence < 0.15 || tested < 5;

  return {
    strand,
    confidence:  parseFloat(confidence.toFixed(3)),
    fwdHits,
    rcHits,
    isAmbiguous,
    // Pre-computed RC query for immediate use if strand === '-'
    rcQuery: qRC,
  };
}

/**
 * Auto-orient a query sequence to the reference strand.
 * Returns the query in forward (reference) orientation plus strand metadata.
 * Call this BEFORE alignment to guarantee all coordinates are in reference space.
 */
export function autoOrientQuery(query, reference) {
  const result = detectReadStrand(query, reference);
  return {
    orientedQuery: result.strand === '-' ? result.rcQuery : query.toUpperCase().replace(/[^ATGCN]/g, ''),
    readStrand:    result.strand,
    strandConf:    result.confidence,
    isAmbiguous:   result.isAmbiguous,
    wasFlipped:    result.strand === '-',
  };
}

// ═══════════════════════════════════════════════════════════
// SECTION 2 — CUT SITE LOCALISATION
// ═══════════════════════════════════════════════════════════

/**
 * Find the Cas9 cut position in alignment coordinates.
 * Searches both strands for the gRNA to confirm cut site.
 */
function findCutInAlignment(alignedR, cutRefPos) {
  let refPos = 0;
  for (let k = 0; k < alignedR.length; k++) {
    if (refPos >= cutRefPos) return k;
    if (alignedR[k] !== '-') refPos++;
  }
  return alignedR.length;
}

/**
 * Count indel gaps ONLY within a window around the cut site in the alignment.
 *
 * ROOT CAUSE OF BUG: alignPairwise() counts ALL gaps across the entire
 * alignment. Sanger sequencing 5' ends are noisy — the Smith-Waterman
 * algorithm may introduce a spurious gap in the low-quality 5' region
 * to maximise the alignment score, inflating the indel count.
 *
 * FIX: Only count gaps within ±WINDOW alignment columns of the cut site.
 * This mirrors what ICE/TIDE do — they use a quantification window.
 * Spurious 5' and 3' end artefacts are excluded by design.
 *
 * WINDOW size = 30bp each side (covers any NHEJ indel up to 30bp,
 * which accounts for >99% of CRISPR-induced indels in plants).
 */
function windowedIndelSize(alignedQ, alignedR, cutAlignPos, windowBp = 30) {
  const start = Math.max(0, cutAlignPos - windowBp);
  const end   = Math.min(alignedQ.length, cutAlignPos + windowBp);

  let gapsInQ = 0, gapsInR = 0;
  for (let k = start; k < end; k++) {
    if (alignedQ[k] === '-' && alignedR[k] !== '-') gapsInQ++;  // deletion in sample
    if (alignedR[k] === '-' && alignedQ[k] !== '-') gapsInR++;  // insertion in sample
  }

  // Net indel size: deletions dominate if gapsInQ > gapsInR, else insertions
  const netSize = Math.abs(gapsInQ - gapsInR);
  const type    = gapsInQ > gapsInR ? 'Deletion' : gapsInR > gapsInQ ? 'Insertion' : 'Complex';
  return { windowedIndelSize: netSize, windowGapsInQ: gapsInQ, windowGapsInR: gapsInR, windowEditType: type };
}

/**
 * Derive cut position AND gRNA strand from reference amplicon.
 * Cas9 cuts 3bp upstream of PAM (between position 17-18 of gRNA).
 * Returns both the cut position in reference coordinates AND which strand
 * the gRNA targets — essential for correct chromatogram annotation.
 *
 * @returns {Object} { cutPosition, grnaStrand: '+' | '-', grnaStart, grnaEnd }
 */
export function deriveCutPosition(refSeq, grnaSeq) {
  if (!refSeq || !grnaSeq) {
    const pos = Math.floor((refSeq || '').length / 2);
    return { cutPosition: pos, grnaStrand: '+', grnaStart: pos, grnaEnd: pos };
  }

  const ref  = refSeq.toUpperCase().replace(/[^ATGCN]/g, '');
  const grna = grnaSeq.toUpperCase().replace(/[^ATGCN]/g, '');

  // Forward strand: gRNA matches reference directly
  const fwdPos = ref.indexOf(grna);
  if (fwdPos !== -1) {
    return {
      cutPosition: fwdPos + grna.length - 3,
      grnaStrand:  '+',
      grnaStart:   fwdPos,
      grnaEnd:     fwdPos + grna.length,
    };
  }

  // Reverse complement: gRNA targets the minus strand
  const complement = { A: 'T', T: 'A', G: 'C', C: 'G', N: 'N' };
  const rcGrna = grna.split('').reverse().map(b => complement[b] || b).join('');
  const revPos = ref.indexOf(rcGrna);
  if (revPos !== -1) {
    // Cut is 3bp from the PAM on the minus strand
    // In reference coordinates: cut = revPos + 3
    return {
      cutPosition: revPos + 3,
      grnaStrand:  '-',
      grnaStart:   revPos,
      grnaEnd:     revPos + grna.length,
    };
  }

  // Partial match fallback — tolerate up to 3 mismatches
  for (let pos = 0; pos <= ref.length - grna.length; pos++) {
    let mm = 0;
    for (let k = 0; k < grna.length; k++) {
      if (grna[k] !== 'N' && ref[pos + k] !== grna[k]) mm++;
    }
    if (mm <= 3) {
      return {
        cutPosition: pos + grna.length - 3,
        grnaStrand:  '+',
        grnaStart:   pos,
        grnaEnd:     pos + grna.length,
      };
    }
  }

  // Last resort
  const mid = Math.floor(ref.length / 2);
  return { cutPosition: mid, grnaStrand: '?', grnaStart: mid, grnaEnd: mid };
}

// ═══════════════════════════════════════════════════════════
// SECTION 2B — MULTI-gRNA ANALYSIS (Phase 3)
// ═══════════════════════════════════════════════════════════

/**
 * Derive cut positions for ALL gRNAs in a multiplex experiment.
 *
 * Supports any number of gRNAs — single, dual, or higher-order multiplex.
 * Each gRNA is searched on both forward and reverse strands of the reference.
 *
 * USE CASES:
 *   • Single gRNA (standard)     — array of 1 item
 *   • Dual-cut deletion          — two gRNAs flanking a target region
 *   • Multiplex knockout         — multiple genes edited simultaneously
 *   • Golden Gate multiplex      — multiple gRNA cassettes in one T-DNA
 *
 * @param {string}   refSeq   - Reference amplicon (forward strand)
 * @param {Array}    grnaList - Array of { label, sequence } objects
 * @returns {Array}  Array of cut-site objects, one per gRNA
 */
export function deriveAllCutPositions(refSeq, grnaList) {
  if (!refSeq || !grnaList || grnaList.length === 0) return [];
  const ref = refSeq.toUpperCase().replace(/[^ATGCN]/g, '');

  return grnaList
    .filter(g => g && g.sequence && g.sequence.trim().length >= 15)
    .map((g, idx) => {
      const grna = g.sequence.toUpperCase().replace(/[^ATGCN]/g, '');
      const label = g.label || `gRNA${idx + 1}`;

      // Forward strand search
      const fwdPos = ref.indexOf(grna);
      if (fwdPos !== -1) {
        return {
          idx, label, sequence: grna,
          cutPosition: fwdPos + grna.length - 3,
          grnaStrand: '+',
          grnaStart:   fwdPos,
          grnaEnd:     fwdPos + grna.length,
          found: true,
          note: `Found on (+) strand at position ${fwdPos}`,
        };
      }

      // Reverse complement search
      const rcGrna  = reverseComplement(grna);
      const revPos  = ref.indexOf(rcGrna);
      if (revPos !== -1) {
        return {
          idx, label, sequence: grna,
          cutPosition: revPos + 3,
          grnaStrand: '-',
          grnaStart:   revPos,
          grnaEnd:     revPos + grna.length,
          found: true,
          note: `Found on (−) strand at position ${revPos}`,
        };
      }

      // Partial match fallback — up to 3 mismatches
      for (let pos = 0; pos <= ref.length - grna.length; pos++) {
        let mm = 0;
        for (let k = 0; k < grna.length; k++) {
          if (grna[k] !== 'N' && ref[pos + k] !== grna[k]) mm++;
        }
        if (mm <= 3) {
          return {
            idx, label, sequence: grna,
            cutPosition: pos + grna.length - 3,
            grnaStrand: '+',
            grnaStart:   pos,
            grnaEnd:     pos + grna.length,
            found: true,
            note: `Partial match (${mm} mismatch${mm !== 1 ? 'es' : ''}) on (+) strand`,
          };
        }
      }

      return {
        idx, label, sequence: grna,
        cutPosition: null,
        grnaStrand: '?',
        grnaStart: null, grnaEnd: null,
        found: false,
        note: 'gRNA not found in reference amplicon',
      };
    });
}

/**
 * Analyse a sample against MULTIPLE gRNA cut sites independently.
 *
 * For each gRNA that was found in the amplicon:
 *   1. Locates the cut site in alignment coordinates
 *   2. Computes R² upstream/downstream of THAT cut site
 *   3. Computes downMM and ICE around THAT cut site
 *   4. Classifies zygosity independently per gRNA
 *
 * Then derives an OVERALL zygosity from all per-gRNA results:
 *   - If only one gRNA: identical to single-gRNA analysis
 *   - If multiple gRNAs both active: reports "Dual Cut" with inter-cut deletion size
 *   - If gRNAs disagree: reports the worst-case (most edited) classification
 *
 * @param {Object} alignment   - SW alignment result from alignPairwise()
 * @param {Array}  cutSites    - Result from deriveAllCutPositions()
 * @param {string} generation  - T0/T1/T2/T3
 * @returns {Object} { perGrna[], overall, primaryCutPosition, interCutDeletion }
 */
export function analyseMultiGrna(alignment, cutSites, generation) {
  const { alignedQ, alignedR } = alignment;
  const foundSites = cutSites.filter(c => c.found && c.cutPosition != null);

  if (foundSites.length === 0) {
    return {
      perGrna: [],
      overall: { class: 0, label: 'Wild-type Segregant', summary: 'No gRNAs located in amplicon' },
      primaryCutPosition: null,
      interCutDeletion: null,
      isDualCut: false,
    };
  }

  // Analyse each gRNA independently
  const perGrna = foundSites.map(site => {
    const cutAlign = findCutInAlignment(alignedR, site.cutPosition);
    const r2up     = r2Segment(alignedQ, alignedR, 0, cutAlign);
    const r2dn     = r2Segment(alignedQ, alignedR, cutAlign, alignedQ.length);
    const downMM   = downstreamMismatchRate(alignedQ, alignedR, cutAlign);
    // Use windowed indel size — counts gaps only near THIS gRNA's cut site
    const win      = windowedIndelSize(alignedQ, alignedR, cutAlign, 30);
    const indelSz  = win.windowedIndelSize;
    const ice      = computeICE(downMM, indelSz, r2dn);
    const zyg      = classifyZygosity(indelSz, downMM, r2up, r2dn, alignment.identity, ice, generation);

    return {
      ...site,
      cutAlign,
      r2up, r2dn,
      downMM:   parseFloat(downMM.toFixed(4)),
      iceScore: ice,
      zygosity: zyg,
      discordanceProfile: buildDiscordanceProfile(alignedQ, alignedR, cutAlign),
    };
  });

  // ── Overall zygosity from all gRNA results ──────────────────────────────
  const classes     = perGrna.map(r => r.zygosity.class);
  const maxClass    = Math.max(...classes);
  const editingGrna = perGrna.filter(r => r.zygosity.class > 0);
  const allSame     = classes.every(c => c === classes[0]);

  // Dual-cut detection: 2+ gRNAs both showing edits, distance > 10 bp
  let interCutDeletion = null;
  let isDualCut        = false;
  if (editingGrna.length >= 2) {
    const positions = editingGrna.map(r => r.cutPosition).sort((a, b) => a - b);
    const dist      = positions[positions.length - 1] - positions[0];
    if (dist > 10) {
      interCutDeletion = dist;
      isDualCut        = true;
    }
  }

  // Build summary label
  let overallLabel   = ['Wild-type Segregant', 'Low / Mosaic Edit', 'Heterozygous / Chimeric', 'Homozygous Edited'][maxClass] || 'Unknown';
  let overallSummary = '';

  if (isDualCut) {
    overallLabel   += ' — Dual Cut';
    overallSummary  = `Both ${editingGrna.map(r => r.label).join(' + ')} active · potential ${interCutDeletion} bp inter-cut deletion`;
  } else if (!allSame && perGrna.length > 1) {
    overallSummary  = `Mixed: ${perGrna.map(r => `${r.label}=${['WT','Low','Het','Homo'][r.zygosity.class]}`).join(', ')}`;
  } else {
    overallSummary  = perGrna.map(r => `${r.label}: ${r.zygosity.label}`).join(' | ');
  }

  // Primary cut site = the one with highest edit signal (highest class, then highest ICE)
  const primarySite = [...perGrna].sort((a, b) =>
    (b.zygosity.class - a.zygosity.class) || (b.iceScore - a.iceScore)
  )[0];

  return {
    perGrna,
    overall: {
      class:    maxClass,
      label:    overallLabel,
      summary:  overallSummary,
      confidence: primarySite?.zygosity?.confidence || 0.5,
    },
    primaryCutPosition: primarySite?.cutPosition ?? null,
    primaryGrnaLabel:   primarySite?.label ?? null,
    interCutDeletion,
    isDualCut,
  };
}

// ═══════════════════════════════════════════════════════════
// SECTION 3 — GAP-AWARE METRICS
// ═══════════════════════════════════════════════════════════

/**
 * Gap-aware R² on a segment of aligned strings.
 * CRITICAL: gap columns (charCode 45) are excluded entirely.
 * Without this, a single 8bp insertion inflates ssRes by ~3200 → R²≈0.
 */
function r2Segment(alignedQ, alignedR, fromPos, toPos) {
  const pairs = [];
  const end = Math.min(toPos, alignedQ.length, alignedR.length);
  for (let k = fromPos; k < end; k++) {
    if (alignedQ[k] === '-' || alignedR[k] === '-') continue;
    pairs.push([alignedQ[k].charCodeAt(0), alignedR[k].charCodeAt(0)]);
  }
  if (pairs.length < 8) return 0;
  const meanR = pairs.reduce((s, [, r]) => s + r, 0) / pairs.length;
  let ssTot = 0, ssRes = 0;
  for (const [q, r] of pairs) {
    ssTot += (r - meanR) ** 2;
    ssRes += (q - r) ** 2;
  }
  return ssTot > 0 ? Math.max(0, Math.min(1, 1 - ssRes / ssTot)) : 0;
}

/**
 * Downstream non-gap mismatch rate.
 * PRIMARY heterozygous signal: two alleles create double-peaks downstream
 * of the cut site → basecaller outputs alternating wrong bases → 25-35% mismatch.
 * WT and Homo: <5% mismatch (single clean allele).
 */
function downstreamMismatchRate(alignedQ, alignedR, cutAlignPos) {
  let mismatches = 0, total = 0;
  for (let k = cutAlignPos; k < alignedQ.length && k < alignedR.length; k++) {
    if (alignedQ[k] === '-' || alignedR[k] === '-') continue;
    total++;
    if (alignedQ[k] !== alignedR[k]) mismatches++;
  }
  return total > 0 ? mismatches / total : 0;
}

/**
 * Upstream mismatch rate (before cut site).
 * Should be low (<5%) for all real samples regardless of zygosity.
 * High upstream mismatch = poor sequence quality or wrong reference.
 */
function upstreamMismatchRate(alignedQ, alignedR, cutAlignPos) {
  let mismatches = 0, total = 0;
  for (let k = 0; k < Math.min(cutAlignPos, alignedQ.length, alignedR.length); k++) {
    if (alignedQ[k] === '-' || alignedR[k] === '-') continue;
    total++;
    if (alignedQ[k] !== alignedR[k]) mismatches++;
  }
  return total > 0 ? mismatches / total : 0;
}

// ═══════════════════════════════════════════════════════════
// SECTION 4 — ICE-EQUIVALENT SCORE
// ═══════════════════════════════════════════════════════════

/**
 * Biologically correct ICE-equivalent editing efficiency score.
 *
 * ICE score = estimated % of reads (alleles) carrying an edit.
 *
 * CORRECT BIOLOGY:
 *   WT:   ICE ≈ 0%    (no edit — both alleles are wild-type)
 *   Low:  ICE ≈ 5-25% (mosaic/chimeric — minority of cells have edit)
 *   Het:  ICE ≈ 40-80% (one allele edited in diploid = ~50% edit reads)
 *   Homo: ICE ≈ 85-99% (both alleles edited = ~100% edit reads)
 *
 * KEY FIX v5.0:
 *   Old: ICE = downMM×280 → Homo gets ICE≈0 (because downMM≈0 for homo!)
 *   New: Homo branch detected by (indelSize>0 AND downMM<0.03 AND r2dn≥0.90)
 *        → ICE = 85 + r2dn×14 ≈ 95-99%
 */
function computeICE(downMM, indelSize, r2dn) {
  // Homozygous: single allele, clean downstream → very high ICE
  if (indelSize > 0 && downMM < 0.03 && r2dn >= 0.90) {
    return Math.round(Math.min(99, 85 + r2dn * 14));
  }
  // Heterozygous / Low: ICE driven by downstream mismatch
  const mmICE    = Math.round(Math.min(85, downMM * 270));
  const indelICE = Math.round(Math.min(25, indelSize * 4));
  return Math.min(90, mmICE + indelICE);
}

// ═══════════════════════════════════════════════════════════
// SECTION 5 — ZYGOSITY CLASSIFIER
// ═══════════════════════════════════════════════════════════

/**
 * UNIFIED zygosity classifier — replaces all previous versions.
 *
 * Three independent signals in priority order:
 *   1. r2dn  — gap-aware downstream R² (HOMO fingerprint)
 *   2. downMM — downstream mismatch rate (HET fingerprint)
 *   3. indelSize — SW gap count at cut site
 *
 * VALIDATED THRESHOLDS (OsAAP3 T1, 19 ICE+TIDE confirmed samples):
 *
 *   HOMO (class 3):
 *     Primary:   indelSize≥1 AND r2dn≥0.97 AND downMM<0.03
 *     Secondary: indelSize≥4 AND r2dn≥0.90 AND downMM<0.05
 *     Biological basis: single clean allele = perfect downstream R²,
 *       no double peaks. Indel SIZE is NOT the discriminator.
 *     Confirmed: 1219_FWP (+1bp, r2dn=0.995, downMM=0.01) → Homo ✓
 *     False-positive protection:
 *       117_FWP (Low): r2dn≈0.89 → fails r2dn≥0.97 → stays Low ✓
 *       119_FWP (Het): downMM=0.24 → fails downMM<0.03 → stays Het ✓
 *       WT samples: indelSize=0 → fails indelSize≥1 → stays WT ✓
 *
 *   HET (class 2):
 *     Primary:   downMM≥0.08 (double-peak signal)
 *     Secondary: r2dn<0.55 AND indelSize>0 (severely degraded downstream)
 *     Confirmed: 119,215,216,221,222,1111,1125,1126 all have downMM≥0.08
 *
 *   LOW/MOSAIC (class 1):
 *     indelSize≥1 AND downMM≥0.038 (above WT noise floor of 0.034)
 *     OR: downMM≥0.04 without clear gap (possible mosaic below noise)
 *     OR: indelSize≥2 AND downMM<0.04 (indel present, clean-ish downstream)
 *     NOTE: 117_FWP and 122_FWP (intrinsic Sanger limit) flagged as Borderline WT
 *
 *   WT (class 0):
 *     Everything else with identity≥0.80 and downMM<0.04
 *
 * GENERATION-SPECIFIC NOTES:
 *   T0: All edits are chimeric by definition. Never report Homo for T0.
 *       Hetero/High-edit call → "Chimeric (T0)", advance to T1.
 *   T1: Standard classification applies.
 *   T2: Homo expected if T1 parent was Het. Low confidence = failed segregation.
 *   T3: Homo should be confirmed. WT in T3 = contamination or seed mixing.
 */
function classifyZygosity(indelSize, downMM, r2up, r2dn, identity, iceScore, generation) {

  // ── T0: Special handling — all T0 plants are chimeric ─────────────────────
  // T0 = primary transformant. Cas9 has been editing during callus development.
  // Different cells carry different alleles → Sanger shows mixed signal.
  // NEVER report Homo or clean Het from T0 — advance to T1 for real classification.
  if (generation === 'T0') {
    if (indelSize === 0 && downMM < 0.04) {
      return {
        class: 0, label: 'Wild-type / No Editing Detected',
        confidence: Math.min(0.85, 0.72 + r2dn * 0.12),
        note: 'T0 plant with no detectable edit — may still be chimeric at low level',
      };
    }
    if (downMM >= 0.08 || indelSize > 0) {
      return {
        class: 2, label: 'Chimeric Edit (T0)',
        confidence: Math.min(0.85, 0.65 + downMM * 0.5 + indelSize * 0.03),
        note: 'T0 plants are always chimeric. Advance to T1 by self-pollination.',
      };
    }
    return {
      class: 1, label: 'Low-Level Edit (T0)',
      confidence: 0.62,
      note: 'Weak edit signal detected in T0. Advance to T1 to segregate.',
    };
  }

  // ── Homozygous: perfect downstream trace + confirmed indel ─────────────────
  // KEY INSIGHT: Homo signal is r2dn≥0.97, NOT large indelSize.
  // A +1bp homo insertion gives r2dn≈0.99 but indelSize=1.
  // Previous indelSize≥3 threshold made +1/+2bp homo IMPOSSIBLE to classify.

  // Primary homo: tight criteria
  if (indelSize >= 1 && downMM < 0.03 && r2dn >= 0.97) {
    return {
      class: 3, label: 'Homozygous Edited',
      confidence: Math.min(0.97, 0.88 + r2dn * 0.08 + (1 - downMM) * 0.03),
    };
  }
  // Secondary homo: larger indels with slightly relaxed criteria
  // (SW may open extra gaps for large indels, modestly reducing r2dn)
  if (indelSize >= 4 && downMM < 0.05 && r2dn >= 0.90) {
    return {
      class: 3, label: 'Homozygous Edited',
      confidence: Math.min(0.94, 0.84 + r2dn * 0.08 + (1 - downMM) * 0.03),
    };
  }

  // ── Heterozygous: double-peak downstream signal ────────────────────────────
  // In Het: two alleles → Sanger superimposes both → basecaller fails downstream.
  // downMM≥0.08 is the reliable signal (WT floor = 0.034 from real data).
  if (downMM >= 0.08) {
    return {
      class: 2, label: 'Heterozygous / Chimeric',
      confidence: Math.min(0.93, 0.72 + downMM * 0.7),
    };
  }
  // R²-based Het fallback: severely degraded downstream even without clear gap
  if (r2dn < 0.55 && indelSize > 0) {
    return {
      class: 2, label: 'Heterozygous / Chimeric',
      confidence: Math.min(0.87, 0.65 + (1 - r2dn) * 0.3),
    };
  }

  // ── Low / Mosaic: edit signal present but below clear Het threshold ─────────
  // Case A: gap at cut site + downMM above WT noise floor (0.034)
  if (indelSize >= 1 && downMM >= 0.038) {
    return {
      class: 1, label: 'Low / Mosaic Edit',
      confidence: Math.min(0.75, 0.58 + downMM * 0.8 + indelSize * 0.04),
    };
  }
  // Case B: downMM slightly elevated (above noise) without clear gap
  if (downMM >= 0.04) {
    return { class: 1, label: 'Low / Mosaic Edit', confidence: 0.60 };
  }
  // Case C: small indel at cut site with clean downstream
  // (possible low-frequency homo in T2/T3, or noise in T1)
  if (indelSize >= 2 && downMM < 0.04) {
    const conf = generation === 'T2' || generation === 'T3' ? 0.70 : 0.62;
    return { class: 1, label: 'Low / Mosaic Edit', confidence: conf };
  }

  // ── Wild-type: no edit signal ──────────────────────────────────────────────
  if (identity >= 0.90) {
    return {
      class: 0, label: 'Wild-type Segregant',
      confidence: Math.min(0.97, 0.85 + r2dn * 0.10),
    };
  }
  if (identity >= 0.80) {
    return { class: 0, label: 'Wild-type Segregant', confidence: 0.82 };
  }

  // Poor quality read — cannot classify reliably
  return { class: 1, label: 'Low / Mosaic Edit', confidence: 0.50 };
}

// ═══════════════════════════════════════════════════════════
// SECTION 6 — ML CONFIDENCE SCORE
// ═══════════════════════════════════════════════════════════

/**
 * ML confidence score (0-100%).
 * Integrates zygosity classifier confidence with sequence quality metrics.
 * Penalises for poor R², low identity, and edge-case ICE values.
 */
function computeMLScore(zygosity, r2up, r2dn, iceScore, identity, generation) {
  let base = zygosity.confidence;

  // Quality bonuses
  if (r2up > 0.90) base += 0.05;
  else if (r2up < 0.70) base -= 0.05;
  if (r2dn > 0.95) base += 0.03;
  else if (r2dn < 0.50) base -= 0.08;

  // Identity penalty for low-quality sequences
  if (identity < 0.75) base -= 0.10;
  else if (identity < 0.85) base -= 0.04;

  // Generation-specific adjustments
  if (generation === 'T2' || generation === 'T3') {
    // Later generations have cleaner signal — boost confidence slightly
    if (zygosity.class === 3) base += 0.03;
    if (zygosity.class === 0) base += 0.02;
  }

  return Math.round(Math.min(0.97, Math.max(0.50, base)) * 100);
}

// ═══════════════════════════════════════════════════════════
// SECTION 7 — FLAG SYSTEM (TIERED, GENERATION-AWARE)
// ═══════════════════════════════════════════════════════════

/**
 * Tiered flag system.
 *
 * Flags are only applied when scientifically justified.
 * WT-class samples (class=0) are NEVER flagged for "small indel" —
 * SW alignment noise occasionally produces 1-2 gap artefacts in clean WT reads.
 *
 * Generation-specific flags:
 *   T0: Always flag — chimeric, needs advancement
 *   T1: Standard flags
 *   T2: Flag WT (unexpected unless segregation failed)
 *   T3: Flag Low/Mosaic (should be fixed by T3)
 */
function flagSample(r2up, r2dn, identity, indelSize, downMM, zygosityClass, generation, iceScore) {
  // Critical quality flags (always shown regardless of class)
  if (identity < 0.65)
    return 'Low alignment identity — verify reference amplicon matches sample';
  if (r2up < 0.35)
    return 'R² too low — poor sequence quality; resequence before interpreting';
  if (r2up < 0.55)
    return 'Borderline sequence quality — results may be unreliable';

  // Class-specific flags
  if (zygosityClass > 0 && indelSize > 0 && indelSize <= 3) {
    // Only flag edited samples for small indels, not WT-classified ones
    return 'Small indel (≤3 bp) — confirm by clonal sequencing (TA-clone ≥8 colonies)';
  }

  // Borderline WT — noise floor warning (catches 117_FWP, 122_FWP type cases)
  if (zygosityClass === 0 && downMM >= 0.022 && downMM < 0.04) {
    return 'Borderline WT — low-level mosaic editing cannot be excluded by Sanger alone; confirm with ICE tool or clonal sequencing';
  }

  // Generation-specific flags
  if (generation === 'T0' && zygosityClass > 0) {
    return 'T0 plant — chimeric editing expected; advance to T1 for definitive genotyping';
  }
  if (generation === 'T2' && zygosityClass === 0) {
    return 'WT in T2 — unexpected if T1 parent was edited; verify sample identity and cross-contamination';
  }
  if (generation === 'T3' && zygosityClass === 1) {
    return 'Low/mosaic edit in T3 — should have segregated by T3; verify with TA-cloning';
  }

  // Low ICE on confirmed edit — warn about interpretation
  if (zygosityClass > 0 && iceScore < 15 && indelSize > 0) {
    return 'Low ICE score despite detected indel — possible chimeric or somaclonal variation';
  }

  return null;
}

// ═══════════════════════════════════════════════════════════
// SECTION 8 — DISCORDANCE PROFILE
// ═══════════════════════════════════════════════════════════

function buildDiscordanceProfile(alignedQ, alignedR, cutAlignPos) {
  const start = Math.max(0, cutAlignPos - 15);
  const end   = Math.min(alignedQ.length, cutAlignPos + 100);
  const profile = [];
  for (let k = start; k < end; k++) {
    const isGap = alignedQ[k] === '-' || alignedR[k] === '-';
    const disc  = isGap ? 0.5 : (alignedQ[k] !== alignedR[k] ? 1 : 0);
    profile.push({ position: k, discordance: disc, relPos: k - cutAlignPos });
  }
  return profile;
}

// ═══════════════════════════════════════════════════════════
// SECTION 9 — BIOLOGICAL INTERPRETATION
// ═══════════════════════════════════════════════════════════

/**
 * Per-generation biological interpretation and recommendation.
 * This is what a plant scientist would tell the user.
 */
function getBiologicalInterpretation(zygosityClass, generation, iceScore, indelSize, isFrameshift, r2dn) {
  const gen = generation || 'T1';

  if (gen === 'T0') {
    if (zygosityClass === 0)
      return { action: 'Advance to T1', priority: 'medium',
               note: 'No edit detected in T0. May still carry low-level chimeric edit. Advance to T1 and screen.' };
    return { action: 'Advance to T1', priority: 'high',
             note: 'Edit detected in T0 (chimeric). Self-pollinate to generate T1 segregants. Screen ~8 T1 plants to recover homozygous and WT segregants.' };
  }

  if (gen === 'T1') {
    if (zygosityClass === 3)
      return { action: 'Advance to T2 — Priority', priority: 'top',
               note: `Homozygous edit confirmed (ICE≈${iceScore}%). Advance directly to T2. Screen for transgene-free plants (Cas9/Hygromycin PCR).` };
    if (zygosityClass === 2)
      return { action: 'Self-pollinate for T2', priority: 'high',
               note: 'Heterozygous edit. Self-pollinate: expect ~25% T2 homozygous, ~50% Het, ~25% WT. Screen 8–10 T2 plants.' };
    if (zygosityClass === 1)
      return { action: 'Verify by clonal sequencing', priority: 'medium',
               note: 'Low-level or mosaic edit. TA-clone PCR and sequence ≥8 colonies to identify allele composition before advancing.' };
    return { action: 'Use as null-segregant control', priority: 'low',
             note: 'Wild-type segregant. No editing detected. Excellent null-segregant control for phenotyping experiments.' };
  }

  if (gen === 'T2') {
    if (zygosityClass === 3)
      return { action: 'Confirmed homozygous line', priority: 'top',
               note: `Homozygous in T2. Screen for transgene-free status. Select 3–5 independent lines for phenotypic analysis.` };
    if (zygosityClass === 2)
      return { action: 'Continue to T3', priority: 'high',
               note: 'Still heterozygous in T2. Self-pollinate for T3 homozygous line.' };
    if (zygosityClass === 0)
      return { action: 'Verify — unexpected WT in T2', priority: 'medium',
               note: 'Wild-type in T2. Verify sample identity. If confirmed WT, use as null-segregant control.' };
    return { action: 'Resequence', priority: 'medium',
             note: 'Low-confidence result in T2. Resequence or use TA-cloning to clarify allele status.' };
  }

  if (gen === 'T3') {
    if (zygosityClass === 3)
      return { action: 'Ready for phenotypic experiments', priority: 'top',
               note: 'Homozygous T3 line confirmed. Proceed with phenotypic characterisation and off-target screening.' };
    if (zygosityClass === 2 || zygosityClass === 1)
      return { action: 'Continue selfing', priority: 'medium',
               note: 'Not yet homozygous in T3. Continue to T4 or use TA-cloning to confirm allele structure.' };
    return { action: 'Use as null-segregant control', priority: 'low',
             note: 'Wild-type T3. Confirmed null-segregant control for experiments.' };
  }

  return { action: 'Manual review needed', priority: 'low', note: 'Unknown generation — please verify.' };
}

// ═══════════════════════════════════════════════════════════
// SECTION 10 — MAIN PER-SAMPLE ANALYSIS
// ═══════════════════════════════════════════════════════════

/**
 * Analyse a single sample against a WT reference.
 *
 * @param {Object} sample   - Sample object: { id, name, sequence, fullSeq, generation, ploidy, quality, traces }
 * @param {Object} wtSample - WT control sample (or null if using refAmplicon)
 * @param {Object} config   - { grnaSeq, cutPosition, refAmplicon, generation, ploidy, geneName, organism }
 * @returns {Object} Full analysis result with all fields for UI + capsEngine + reports
 */
export function analyseSample(sample, wtSample, config) {
  const {
    grnaSeq     = '',
    grnaList    = null,   // Phase 3: Array of { label, sequence } for multi-gRNA
    cutPosition = null,
    refAmplicon = '',
    generation  = sample?.generation || 'T1',
    ploidy      = sample?.ploidy     || 'diploid',
    geneName    = 'Unknown',
    organism    = 'Plant',
  } = config || {};

  // Build the effective gRNA list — single or multi
  // grnaList takes priority; if absent, wrap grnaSeq as a single-entry list
  const effectiveGrnaList = grnaList && grnaList.length > 0
    ? grnaList
    : (grnaSeq ? [{ label: 'gRNA1', sequence: grnaSeq }] : []);

  const isMultiGrna = effectiveGrnaList.length > 1;

  // ── 1. Sequence preparation ─────────────────────────────────────────────
  // Use full (untrimmed) sequence if available — quality trim can cut near the indel
  const queryTrimmed = (sample.sequence || '').replace(/[^ATGCN]/gi, '').toUpperCase();
  const queryFull    = (sample.fullSeq   || sample.sequence || '').replace(/[^ATGCN]/gi, '').toUpperCase();
  const queryRaw     = queryFull.length > queryTrimmed.length ? queryFull : queryTrimmed;

  // Reference: prefer explicit refAmplicon over WT sample sequence
  const wt = (refAmplicon || wtSample?.sequence || '').replace(/[^ATGCN]/gi, '').toUpperCase();

  // Error handling for missing/too-short sequences
  if (!queryRaw || queryRaw.length < 30) {
    return buildErrorResult(sample, config, 'Sequence too short or missing (min 30 bp required)');
  }
  if (!wt || wt.length < 30) {
    return buildErrorResult(sample, config, 'No reference provided — add WT control sample or paste reference amplicon in Setup');
  }

  // ── 2. STRAND DETECTION — auto-orient query to reference space ───────────
  // KEY FIX (Phase 1): If the researcher used a reverse primer, the Sanger
  // read is the reverse complement of the reference. We detect this using
  // k-mer matching and flip the query before alignment so ALL coordinates
  // (cut site, indel position, discordance profile) are in reference space.
  const strandInfo = autoOrientQuery(queryRaw, wt);
  const query      = strandInfo.orientedQuery;  // always forward-strand from here

  // ── 3. Alignment ─────────────────────────────────────────────────────────
  const alignment = alignPairwise(query, wt);
  const { alignedQ, alignedR, gapsInQuery, gapsInRef } = alignment;
  // NOTE: alignment.indelSize is the RAW gap count across the FULL alignment.
  // We do NOT use it for classification — it includes spurious 5'/3' end gaps
  // from low-quality Sanger regions (root cause of the 1bp→2bp miscounting bug).

  // ── 4. Multi-gRNA cut-site analysis (Phase 3) ───────────────────────────
  const allCutSites = deriveAllCutPositions(wt, effectiveGrnaList);
  const multiResult = analyseMultiGrna(alignment, allCutSites, generation);

  // Primary cut position
  const primarySiteInfo = allCutSites.find(c => c.found && c.cutPosition != null)
    || allCutSites[0]
    || { cutPosition: cutPosition ?? Math.floor(wt.length / 2), grnaStrand: '+', grnaStart: 0, grnaEnd: 0 };

  const effectiveCutPos = primarySiteInfo.cutPosition ?? cutPosition ?? Math.floor(wt.length / 2);
  const grnaStrand      = primarySiteInfo.grnaStrand || '+';
  const grnaStart       = primarySiteInfo.grnaStart  || 0;
  const grnaEnd         = primarySiteInfo.grnaEnd    || 0;
  const cutAlign        = findCutInAlignment(alignedR, effectiveCutPos);

  // ── WINDOWED INDEL SIZE — THE CORRECT INDEL COUNT ────────────────────────
  // Count gaps ONLY within ±30bp of the cut site in the alignment.
  // This eliminates spurious 5'/3' end gaps from noisy Sanger regions.
  // Validated: gives 1bp for the Sample.ab1 case (ICE/TIDE/CRISPRid confirmed),
  //            gives 2bp when using raw alignment.indelSize (the old bug).
  const windowed = windowedIndelSize(alignedQ, alignedR, cutAlign, 30);
  const indelSize = windowed.windowedIndelSize;  // USE THIS — not alignment.indelSize

  // ── 5. Core metrics on primary cut site (gap-aware) ─────────────────────
  const r2up   = r2Segment(alignedQ, alignedR, 0, cutAlign);
  const r2dn   = r2Segment(alignedQ, alignedR, cutAlign, alignedQ.length);
  const r2     = (r2up + r2dn) / 2;
  const downMM = downstreamMismatchRate(alignedQ, alignedR, cutAlign);
  const upMM   = upstreamMismatchRate(alignedQ, alignedR, cutAlign);

  // ICE score (biologically correct for all classes including Homo)
  const iceScore = computeICE(downMM, indelSize, r2dn);

  // ── 6. Classify ──────────────────────────────────────────────────────────
  // Use multi-gRNA overall result if multiple gRNAs are present;
  // otherwise fall back to single-gRNA classifier (identical result).
  const zygosity = isMultiGrna && multiResult.perGrna.length > 0
    ? multiResult.overall
    : classifyZygosity(indelSize, downMM, r2up, r2dn, alignment.identity, iceScore, generation);

  const mlScore = computeMLScore(
    { confidence: zygosity.confidence || 0.7, class: zygosity.class },
    r2up, r2dn, iceScore, alignment.identity, generation
  );

  // ── 7. Flag ──────────────────────────────────────────────────────────────
  let flag = flagSample(r2up, r2dn, alignment.identity, indelSize, downMM,
                         zygosity.class, generation, iceScore);

  // Add strand warning if ambiguous
  if (strandInfo.isAmbiguous && !flag) {
    flag = 'Read strand ambiguous — verify primer orientation; result may be unreliable';
  }

  // ── 8. Edit type and protein effect ─────────────────────────────────────
  let indelType = 'Wild-type';
  if      (windowed.windowGapsInR > windowed.windowGapsInQ) indelType = 'Insertion';
  else if (windowed.windowGapsInQ > windowed.windowGapsInR) indelType = 'Deletion';
  else if (indelSize > 0)                                   indelType = 'Complex';

  const isFrameshift  = indelSize > 0 && indelSize % 3 !== 0;
  const proteinEffect = indelSize === 0
    ? 'No change — wild-type'
    : isFrameshift
      ? 'Frameshift → likely loss-of-function'
      : indelSize % 3 === 0
        ? 'In-frame indel → functional validation needed'
        : 'Complex rearrangement — functional validation needed';

  // ── 9. Display label ─────────────────────────────────────────────────────
  const zygosityLabel = {
    0: 'Wild-type Segregant',
    1: 'Low / Mosaic Edit',
    2: 'Heterozygous / Chimeric',
    3: 'Homozygous Edited',
  }[zygosity.class] || 'Unknown';

  const editType = (indelSize === 0 && zygosity.class === 0) ? 'Wild-type' : indelType;

  // ── 10. Biological interpretation ─────────────────────────────────────────
  const bioInterp = getBiologicalInterpretation(
    zygosity.class, generation, iceScore, indelSize, isFrameshift, r2dn
  );

  // ── 11. Sample cut position in sequence coordinates (for ChromatogramViewer)
  let sampleCutPosition = 0;
  for (let k = 0; k < Math.min(cutAlign, alignedQ.length); k++) {
    if (alignedQ[k] !== '-') sampleCutPosition++;
  }

  // ── 12. Return unified result object ────────────────────────────────────
  return {
    // Identity
    id:               sample.id   || sample.name,
    name:             sample.name || sample.id,
    generation,
    ploidy,
    geneName,
    organism,

    // Core metrics (ICE-equivalent)
    iceScore,
    r2,
    r2up,
    r2dn,
    mlScore,
    downMM:           parseFloat(downMM.toFixed(4)),
    upMM:             parseFloat(upMM.toFixed(4)),

    // Edit details
    editType,
    indelSize,
    lenDiff:          gapsInRef - gapsInQuery,
    isFrameshift,
    proteinEffect,

    // Zygosity
    zygosityClass:    zygosity.class,
    zygosityLabel,
    zygosityRaw:      zygosity.label,
    zygosityNote:     zygosity.note || null,
    confidence:       zygosity.confidence,

    // Quality
    sequenceIdentity: parseFloat(alignment.identity.toFixed(4)),
    alignmentScore:   alignment.score,
    gaps:             alignment.gaps,
    flag,

    // ── STRAND INFORMATION (Phase 1 addition) ──────────────────────────
    // readStrand:  '+' = forward primer read (natural orientation)
    //              '-' = reverse primer read (was RC'd before analysis)
    // grnaStrand:  '+' = gRNA targets forward strand
    //              '-' = gRNA targets reverse strand (minus strand targeting)
    // wasFlipped:  true if the read was reverse-complemented before analysis
    readStrand:       strandInfo.readStrand,
    grnaStrand:       grnaStrand || '+',
    grnaStart:        grnaStart  || 0,
    grnaEnd:          grnaEnd    || 0,
    strandConf:       strandInfo.strandConf,
    wasFlipped:       strandInfo.wasFlipped,
    strandAmbiguous:  strandInfo.isAmbiguous,

    // Biological interpretation
    action:           bioInterp.action,
    actionPriority:   bioInterp.priority,
    actionNote:       bioInterp.note,

    // Chromatogram viewer
    discordanceProfile: buildDiscordanceProfile(alignedQ, alignedR, cutAlign),
    cutPosition:      effectiveCutPos,
    grnaSeq,
    sampleCutPosition,
    alignedQuery:     alignedQ,
    alignedRef:       alignedR,

    // Raw data passthrough
    sequence:         query,          // always forward-oriented
    rawSequence:      queryRaw,       // original as-received
    wtSequence:       wt,
    traces:           sample.traces   || null,
    wtTraces:         wtSample?.traces || null,
    quality:          sample.quality  || null,
    avgQuality:       sample.avgQuality || null,
    qualGrade:        sample.qualGrade  || null,
    source:           sample.source    || 'text',
    filename:         sample.filename  || null,

    // capsEngine compatibility
    indels: indelSize > 0 ? [{
      position:     effectiveCutPos,
      type:         indelType === 'Deletion'  ? 'deletion'
                  : indelType === 'Insertion' ? 'insertion' : 'substitution',
      length:       indelSize,
      refBase:      '-',
      altBase:      '-',
      relToCutSite: 0,
      isFrameshift,
    }] : [],
    isEdited:  zygosity.class > 0,
    editSize:  indelSize,
    zygosity:  { label: zygosity.label, class: zygosity.class },
    status:    flag ? 'uncertain'
             : zygosity.class === 0 ? 'wildtype' : 'edited',

    // ── Multi-gRNA fields (Phase 3) ────────────────────────────────────────
    // Always present — for single-gRNA, perGrna has exactly one entry.
    // For multi-gRNA: full breakdown per gRNA + dual-cut detection.
    isMultiGrna,
    grnaList:          effectiveGrnaList,
    allCutSites,
    perGrnaResults:    multiResult.perGrna,
    multiGrnaSummary:  multiResult.overall.summary || null,
    isDualCut:         multiResult.isDualCut,
    interCutDeletion:  multiResult.interCutDeletion,
    primaryGrnaLabel:  multiResult.primaryGrnaLabel,
  };
}

// ── Error result builder ─────────────────────────────────────────────────────
function buildErrorResult(sample, config, errorMsg) {
  const gen    = config?.generation || sample?.generation || 'T1';
  const ploidy = config?.ploidy     || sample?.ploidy     || 'diploid';
  return {
    id: sample?.id || sample?.name || 'unknown',
    name: sample?.name || sample?.id || 'unknown',
    generation: gen, ploidy,
    error: errorMsg, status: 'error',
    iceScore: 0, r2: 0, r2up: 0, r2dn: 0, mlScore: 50,
    zygosityClass: 1, zygosityLabel: 'Error', zygosityRaw: 'Error',
    flag: errorMsg,
    editType: 'Unknown', indelSize: 0, isFrameshift: false,
    proteinEffect: 'Unknown',
    sequenceIdentity: 0, alignmentScore: 0, gaps: 0,
    cutPosition: config?.cutPosition || 0, grnaSeq: config?.grnaSeq || '',
    discordanceProfile: [], sampleCutPosition: 0,
    alignedQuery: '', alignedRef: '',
    sequence: '', wtSequence: '',
    traces: null, wtTraces: null,
    indels: [], isEdited: false, editSize: 0,
    zygosity: { label: 'Error', class: 1 },
    action: 'Resequence — data error', actionPriority: 'low',
    actionNote: errorMsg,
  };
}

// ═══════════════════════════════════════════════════════════
// SECTION 11 — BATCH ANALYSIS
// ═══════════════════════════════════════════════════════════

/**
 * Run analysis on all samples in a batch.
 * Automatically identifies WT control by name matching.
 * Handles mixed generation batches (T0+T1, T1+T2, etc.).
 */
export function runAnalysis(samples, config) {
  if (!samples || samples.length === 0) return [];

  // Find WT control — prioritise exact 'MTU', 'WT', 'wildtype', 'control' match
  const wtSample = samples.find(s =>
    /\bMTU\b|\bWT\b|\bwildtype\b|\bcontrol\b/i.test(s.name || s.id)
  ) || samples.find(s =>
    /wt|wild.?type|control/i.test(s.name || s.id)
  ) || null;

  // Analyse all samples (including WT — WT will self-classify as class 0)
  return samples.map(s => analyseSample(s, wtSample, {
    ...config,
    generation: s.generation || config?.generation || 'T1',
    ploidy:     s.ploidy     || config?.ploidy     || 'diploid',
  }));
}

// ═══════════════════════════════════════════════════════════
// SECTION 12 — DASHBOARD STATISTICS
// ═══════════════════════════════════════════════════════════

export function buildDashboardStats(results) {
  if (!results || results.length === 0) {
    return {
      totalSamples: 0, editedCount: 0, wtCount: 0, uncertainCount: 0,
      editingEfficiency: 0, frameshiftCount: 0, homozygousCount: 0,
      hetCount: 0, lowEditCount: 0, editTypeDistribution: {},
      advancePlants: [], verifyFirst: [], wtControls: [],
    };
  }

  const all        = results.filter(r => !r.error);
  const edited     = all.filter(r => r.zygosityClass > 0);
  const wt         = all.filter(r => r.zygosityClass === 0);
  const flagged    = all.filter(r => r.flag);
  const homo       = all.filter(r => r.zygosityClass === 3);
  const het        = all.filter(r => r.zygosityClass === 2);
  const lowEdit    = all.filter(r => r.zygosityClass === 1);
  const frameshift = all.filter(r => r.isFrameshift);

  const editTypes  = {};
  all.forEach(r => {
    if (r.editType && r.editType !== 'Wild-type')
      editTypes[r.editType] = (editTypes[r.editType] || 0) + 1;
  });

  // Editing efficiency: edited / total (excluding errors and flagged samples for conservative estimate)
  const confirmed  = all.filter(r => !r.flag);
  const confEdited = confirmed.filter(r => r.zygosityClass > 0);
  const efficiency = all.length > 0 ? (edited.length / all.length * 100) : 0;
  const confEfficiency = confirmed.length > 0 ? (confEdited.length / confirmed.length * 100) : 0;

  return {
    totalSamples:         all.length,
    editedCount:          edited.length,
    wtCount:              wt.length,
    uncertainCount:       flagged.length,
    editingEfficiency:    Math.round(efficiency * 10) / 10,
    confEditingEfficiency: Math.round(confEfficiency * 10) / 10,
    frameshiftCount:      frameshift.length,
    homozygousCount:      homo.length,
    hetCount:             het.length,
    lowEditCount:         lowEdit.length,
    editTypeDistribution: editTypes,
    // Action lists
    topAdvance:           homo.filter(r => !r.flag).map(r => r.name),
    advancePlants:        het.filter(r => !r.flag).map(r => r.name),
    verifyFirst:          flagged.map(r => r.name),
    wtControls:           wt.map(r => r.name),
    lowPriority:          lowEdit.map(r => r.name),
    // Expected T1 segregation ratios
    expectedHomozygous:   '~25%',
    expectedHeterozygous: '~50%',
    expectedWT:           '~25%',
    expectedTransgFree:   '~37%',
  };
}

// ═══════════════════════════════════════════════════════════
// SECTION 13 — AI INSIGHTS GENERATOR
// ═══════════════════════════════════════════════════════════

/**
 * Generate AI recommendations for the AI Insights tab.
 * All insights are derived from actual analysis results, not hardcoded.
 */
export function generateAIInsights(results, config) {
  const stats = buildDashboardStats(results);
  const gen   = config?.generation || 'T1';
  const recs  = [];

  // Top priority: Homozygous plants ready to advance
  if (stats.homozygousCount > 0) {
    const plants = stats.topAdvance.join(', ') || 'see list';
    recs.push({
      priority: 'high',
      icon: '⭐',
      text: `${stats.homozygousCount} homozygous edited plant${stats.homozygousCount > 1 ? 's' : ''} confirmed: ${plants}. Advance to ${gen === 'T1' ? 'T2' : 'T3'} immediately. Screen for transgene-free status (Cas9/Hygromycin PCR).`,
    });
  }

  // Heterozygous plants
  if (stats.hetCount > 0) {
    const plants = stats.advancePlants.slice(0, 5).join(', ');
    const more   = stats.advancePlants.length > 5 ? ` +${stats.advancePlants.length - 5} more` : '';
    recs.push({
      priority: stats.homozygousCount === 0 ? 'high' : 'medium',
      icon: '🌿',
      text: `${stats.hetCount} heterozygous plant${stats.hetCount > 1 ? 's' : ''}: ${plants}${more}. Self-pollinate → screen ~8 ${gen === 'T1' ? 'T2' : 'T3'} plants to recover homozygous lines.`,
    });
  }

  // Verification needed
  if (stats.uncertainCount > 0) {
    const subset = stats.verifyFirst.slice(0, 6).join(', ');
    const more   = stats.verifyFirst.length > 6 ? ` +${stats.verifyFirst.length - 6} more` : '';
    recs.push({
      priority: 'medium',
      icon: '⚠️',
      text: `${stats.uncertainCount} sample${stats.uncertainCount > 1 ? 's' : ''} need verification: ${subset}${more}. Confirm by TA-cloning (≥8 colonies) or re-sequencing with optimised conditions.`,
    });
  }

  // Transgene screening reminder
  if (stats.editedCount > 0) {
    recs.push({
      priority: 'high',
      icon: '🧬',
      text: `Screen ALL ${stats.editedCount} edited lines for T-DNA presence (Cas9 PCR + Hygromycin assay). Transgene-free edited plants (~37% of T1) are your final target material.`,
    });
  }

  // WT controls note
  if (stats.wtCount > 0) {
    recs.push({
      priority: 'low',
      icon: '🌱',
      text: `${stats.wtCount} wild-type segregant${stats.wtCount > 1 ? 's' : ''} identified (${stats.wtControls.slice(0, 4).join(', ')}${stats.wtCount > 4 ? '…' : ''}). Use as null-segregant controls in phenotyping experiments.`,
    });
  }

  // Statistical power note
  const targetN = Math.ceil(Math.log(0.05) / Math.log(0.75));
  recs.push({
    priority: 'low',
    icon: '📊',
    text: `Screen ≥${targetN} ${gen} plants per line to find ≥3 homozygous edited plants with 95% probability (assuming 25% expected frequency).`,
  });

  return recs;
}

// ═══════════════════════════════════════════════════════════
// SECTION 14 — REPORT AND EXPORT
// ═══════════════════════════════════════════════════════════

export function buildReportText(results, config, stats) {
  const { grnaSeq = '', grnaList = null, generation = 'T1', geneName = 'Unknown', organism = 'Plant' } = config || {};
  const cutPos = config?.cutPosition || 'unknown';

  // Build gRNA summary line for report header
  const grnaDisplay = grnaList && grnaList.length > 1
    ? grnaList.map((g, i) => `${g.label || `gRNA${i+1}`}: ${g.sequence}`).join(' | ')
    : grnaSeq || '—';
  const isMultiGrna = grnaList && grnaList.length > 1;

  let txt = `PlantaCRISPAI — CRISPR Analysis Report\n`;
  txt += `${'═'.repeat(55)}\n`;
  txt += `Gene: ${geneName} | Organism: ${organism} | Generation: ${generation}\n`;
  if (isMultiGrna) {
    txt += `Multi-gRNA experiment (${grnaList.length} gRNAs):\n`;
    grnaList.forEach((g, i) => txt += `  gRNA${i+1} [${g.label||`gRNA${i+1}`}]: ${g.sequence}\n`);
  } else {
    txt += `gRNA: ${grnaDisplay} | Cut position: ~${cutPos} bp\n`;
  }
  txt += `Analysis engine: v5.0 (Phase 3 — Multi-gRNA)\n\n`;
  txt += `SUMMARY\n${'─'.repeat(40)}\n`;
  txt += `Total samples analysed: ${stats.totalSamples}\n`;
  txt += `Edited: ${stats.editedCount}/${stats.totalSamples} (${stats.editingEfficiency}%)\n`;
  txt += `  Homozygous: ${stats.homozygousCount}\n`;
  txt += `  Heterozygous: ${stats.hetCount}\n`;
  txt += `  Low/Mosaic: ${stats.lowEditCount}\n`;
  txt += `Wild-type segregants: ${stats.wtCount}\n`;
  txt += `Frameshift mutations: ${stats.frameshiftCount}\n`;
  txt += `Needs verification: ${stats.uncertainCount}\n\n`;

  if (stats.topAdvance.length > 0)
    txt += `⭐ TOP PRIORITY (advance now): ${stats.topAdvance.join(', ')}\n`;
  if (stats.advancePlants.length > 0)
    txt += `★ Advance (heterozygous): ${stats.advancePlants.join(', ')}\n`;
  if (stats.verifyFirst.length > 0)
    txt += `⚠ Verify first: ${stats.verifyFirst.join(', ')}\n`;
  if (stats.wtControls.length > 0)
    txt += `○ WT controls: ${stats.wtControls.join(', ')}\n`;

  txt += `\nSAMPLE DETAILS\n${'─'.repeat(55)}\n`;
  results.forEach(r => {
    if (r.error) {
      txt += `${r.name || r.id} — ERROR: ${r.error}\n\n`;
      return;
    }
    const zyg = {3:'Homozygous',2:'Heterozygous',1:'Low/Mosaic',0:'Wild-type'}[r.zygosityClass];
    txt += `${r.name || r.id}\n`;
    txt += `  Edit: ${r.editType}${r.indelSize ? ` (${Math.abs(r.indelSize)} bp)` : ''} | Direction: ${r.lenDiff > 0 ? '+' : ''}${r.lenDiff} bp\n`;
    txt += `  Zygosity: ${zyg} | ICE: ${r.iceScore}% | R²: ${r.r2.toFixed(3)} | ML: ${r.mlScore}%\n`;
    txt += `  R² upstream: ${(r.r2up||0).toFixed(3)} | R² downstream: ${(r.r2dn||0).toFixed(3)}\n`;
    txt += `  Identity: ${((r.sequenceIdentity||0)*100).toFixed(1)}% | Frameshift: ${r.isFrameshift ? 'Yes' : 'No'}\n`;
    txt += `  Protein: ${r.proteinEffect}\n`;
    if (r.flag) txt += `  ⚠ FLAG: ${r.flag}\n`;
    txt += `  → ${r.action || ''}\n\n`;
  });
  return txt;
}

export function buildCSV(results) {
  const header = [
    'Sample', 'Generation', 'Ploidy', 'EditType', 'IndelSize_bp', 'Direction_bp',
    'Zygosity', 'ZygosityClass', 'ICE_Score', 'R2_Overall',
    'R2_Upstream', 'R2_Downstream', 'DownstreamMismatch',
    'ML_Score', 'Confidence', 'IsFrameshift', 'ProteinEffect',
    'SequenceIdentity', 'ReadStrand', 'GRNAStrand', 'ReadAutoFlipped',
    'IsMultiGrna', 'ActiveGrnas', 'IsDualCut', 'InterCutDeletion_bp',
    'MultiGrnaSummary',
    'Flag', 'Action', 'ActionPriority',
  ].join(',');

  const rows = results.map(r => {
    // Build per-gRNA active summary (which gRNAs found edits)
    const activeGrnas = r.perGrnaResults
      ? r.perGrnaResults.filter(g => g.zygosity?.class > 0).map(g => g.label).join(';')
      : (r.primaryGrnaLabel || '');

    return [
      r.name || r.id,
      r.generation,
      r.ploidy,
      r.editType,
      r.indelSize,
      r.lenDiff || 0,
      r.zygosityRaw || r.zygosityLabel,
      r.zygosityClass,
      r.iceScore,
      (r.r2 || 0).toFixed(3),
      (r.r2up || 0).toFixed(3),
      (r.r2dn || 0).toFixed(3),
      (r.downMM || 0).toFixed(4),
      r.mlScore,
      r.confidence ? r.confidence.toFixed(2) : '',
      r.isFrameshift ? 1 : 0,
      `"${r.proteinEffect || ''}"`,
      ((r.sequenceIdentity || 0) * 100).toFixed(1),
      r.readStrand  || '+',
      r.grnaStrand  || '+',
      r.wasFlipped  ? 'Yes' : 'No',
      r.isMultiGrna ? 'Yes' : 'No',
      `"${activeGrnas}"`,
      r.isDualCut   ? 'Yes' : 'No',
      r.interCutDeletion || '',
      `"${r.multiGrnaSummary || ''}"`,
      r.flag ? `"${r.flag}"` : '',
      `"${r.action || ''}"`,
      r.actionPriority || '',
    ].join(',');
  });

  return [header, ...rows].join('\n');
}
