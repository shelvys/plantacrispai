/* eslint-disable no-unused-vars */
// ============================================================
// FILE: src/engine/grnaDesigner.js  (v2 — ML-powered)
// PURPOSE: CAPS-aware gRNA designer with proper ML scoring
//
// UPGRADE from v1 (rule-based) to v2 (feature-engineered ML):
//   v1: if GC > 40 → +15 points (naive rules)
//   v2: 28-feature vector → 4 sub-models → ensemble prediction
//       with thermodynamic stability, position-specific weights,
//       structural penalties, and plant-specific corrections
//
// Scientific basis:
//   Doench et al. 2016 (Nature Biotechnology) — PWM model
//   Xu et al. 2015 (Plant Cell) — rice-specific features
//   Endo et al. 2015 (Plant Cell) — arabidopsis features
//   Tang et al. 2018 (Nature Plants) — wheat features
//   Sugimoto et al. 1995 — RNA:DNA thermodynamic parameters
// ============================================================

import { findRESites, scoreGRNAForCAPS,
         reverseComplement, RE_DATABASE } from './capsEngine';
import { predictGRNAEfficiency }          from './grnaFeatures';

// ══════════════════════════════════════════════════════════
// SECTION 1 — gRNA CANDIDATE FINDER
// ══════════════════════════════════════════════════════════

/**
 * Find all valid Cas9 target sites (20bp + NGG PAM)
 * on both strands of the input sequence
 */
export function findAllGRNACandidates(sequence,
                                       targetStart = 0,
                                       targetEnd = null) {
  const seq    = sequence.toUpperCase().replace(/[^ATGCN]/g, '');
  const end    = targetEnd || seq.length;
  const LEN    = 20;
  const cands  = [];
  const seen   = new Set();

  const addCandidate = (grna, pam, strand, pos, cutPos) => {
    const key = `${grna}:${strand}:${pos}`;
    if (seen.has(key)) return;
    seen.add(key);
    // Skip sequences with too many N
    if ((grna.match(/N/g)||[]).length > 2) return;
    // Skip if grna contains only one base type
    const uniq = new Set(grna.split('')).size;
    if (uniq < 2) return;
    cands.push({ sequence: grna, pam, strand, position: pos, cutPosition: cutPos });
  };

  // Forward strand — find NGG
  for (let i = targetStart; i <= end - LEN - 3; i++) {
    const pam = seq.slice(i + LEN, i + LEN + 3);
    if (pam[1] === 'G' && pam[2] === 'G') {
      addCandidate(seq.slice(i, i + LEN), pam, '+', i, i + LEN - 3);
    }
  }

  // Reverse strand — find CCN (= NGG on RC)
  for (let i = targetStart; i <= end - LEN - 3; i++) {
    const pamRC = seq.slice(i, i + 3);
    if (pamRC[0] === 'C' && pamRC[1] === 'C') {
      const grnaFwd = reverseComplement(seq.slice(i + 3, i + 3 + LEN));
      const pam     = reverseComplement(pamRC);
      const cutPos  = seq.length - i - 3;
      addCandidate(grnaFwd, pam, '-', seq.length - i - LEN - 3, cutPos);
    }
  }

  return cands.sort((a,b) => a.position - b.position);
}

// ══════════════════════════════════════════════════════════
// SECTION 2 — OFF-TARGET RISK ASSESSMENT
// More sophisticated than v1
// ══════════════════════════════════════════════════════════

/**
 * Assess off-target risk using sequence-based features
 * Full off-target prediction requires genome alignment
 * (BLAST/Bowtie) — this is an estimate from sequence alone
 */
function assessOffTargetRisk(grnaSeq, species = 'rice') {
  const seq  = grnaSeq.toUpperCase();
  let score  = 0; // 0 = low risk, 100 = high risk

  // 1. Seed region uniqueness proxy
  // Low complexity seed = more off-targets
  const seed    = seq.slice(0, 12);
  const seedUniq = new Set(seed.split('')).size;
  if (seedUniq <= 2) score += 30;
  else if (seedUniq <= 3) score += 10;

  // 2. Homopolymer runs indicate repeat regions
  for (const base of 'ATGC') {
    const maxRun = (seq.match(new RegExp(`${base}+`, 'g'))||[])
      .reduce((a,b) => Math.max(a, b.length), 0);
    if (maxRun >= 5) score += 20;
    else if (maxRun >= 4) score += 10;
    else if (maxRun >= 3) score += 4;
  }

  // 3. GC extremes increase off-targets
  const gc = (seq.match(/[GC]/g)||[]).length / seq.length * 100;
  if (gc < 20 || gc > 85) score += 20;
  else if (gc < 30 || gc > 75) score += 8;

  // 4. Known problematic sequences
  if (seq.includes('TTTT')) score += 15;
  if (seq.includes('AAAA')) score += 8;

  // 5. Species genome complexity adjustment
  // Rice has ~370Mbp diploid — smaller than wheat (17Gbp hexaploid)
  const genomeSizePenalty = {
    'rice': 0, 'arabidopsis': 0, 'maize': 5,
    'tomato': 3, 'wheat': 15, 'soybean': 8,
    'potato': 8, 'other': 5
  }[species] || 5;
  score += genomeSizePenalty;

  score = Math.min(100, Math.max(0, score));

  return {
    score,
    level: score <= 25 ? 'Low'
         : score <= 50 ? 'Medium'
         : score <= 75 ? 'High' : 'Very High',
    color: score <= 25 ? '#00D4AA'
         : score <= 50 ? '#4D9FFF'
         : score <= 75 ? '#F0C040' : '#FF4D6D',
    note: score <= 25
      ? 'Low off-target risk — good candidate'
      : score <= 50
      ? 'Moderate off-target risk — validate by sequencing'
      : 'High off-target risk — recommend genome-wide off-target check',
  };
}

// ══════════════════════════════════════════════════════════
// SECTION 3 — PARETO FRONTIER ANALYSIS
// Find optimal gRNAs that balance efficiency AND CAPS
// ══════════════════════════════════════════════════════════

/**
 * Find Pareto-optimal gRNAs
 * A gRNA is Pareto-optimal if no other gRNA beats it
 * on BOTH efficiency AND CAPS score simultaneously
 *
 * Example:
 *   gRNA-A: efficiency=85, CAPS=60  → Pareto optimal
 *   gRNA-B: efficiency=70, CAPS=80  → Pareto optimal
 *   gRNA-C: efficiency=65, CAPS=55  → NOT optimal (A beats on both)
 */
function findParetoOptimal(candidates) {
  return candidates.map((cand, i) => {
    const dominated = candidates.some((other, j) => {
      if (i === j) return false;
      return other.scores.efficiency     >= cand.scores.efficiency &&
             other.scores.capsCompat     >= cand.scores.capsCompat &&
             (other.scores.efficiency > cand.scores.efficiency ||
              other.scores.capsCompat > cand.scores.capsCompat);
    });
    return { ...cand, isPareto: !dominated };
  });
}

// ══════════════════════════════════════════════════════════
// SECTION 4 — MAIN SCORING + RANKING PIPELINE
// ══════════════════════════════════════════════════════════

/**
 * Score and rank all gRNA candidates using ML features
 * @param {Array} candidates - from findAllGRNACandidates
 * @param {string} refSequence - reference amplicon
 * @param {Object} options
 * @returns {Array} ranked candidates with full ML scores
 */
export function rankGRNACandidates(candidates, refSequence, options = {}) {
  const {
    efficiencyWeight = 0.50,
    capsWeight       = 0.35,
    offTargetWeight  = 0.15,
    minEfficiency    = 30,
    species          = 'rice',
  } = options;

  const scored = [];

  for (const cand of candidates) {
    // ── ML efficiency prediction ─────────────────────
    const mlResult = predictGRNAEfficiency(
      cand.sequence,
      cand.pam,
      refSequence,
      cand.cutPosition,
      species
    );

    if (mlResult.score < minEfficiency) continue;

    // ── CAPS compatibility ───────────────────────────
    let capsData;
    try {
      capsData = scoreGRNAForCAPS(
        cand.sequence, refSequence, cand.cutPosition
      );
    } catch (e) {
      capsData = { capsCompatibilityScore: 0, hasCAPSCompatibility: false,
                   capsCategory: 'None', bestEnzyme: null };
    }

    // ── Off-target risk ──────────────────────────────
    const offTarget = assessOffTargetRisk(cand.sequence, species);
    const offTargetScore = 100 - offTarget.score; // invert: high = good

    // ── Combined score ───────────────────────────────
    const combined = Math.round(
      mlResult.score                   * efficiencyWeight +
      capsData.capsCompatibilityScore  * capsWeight +
      offTargetScore                   * offTargetWeight
    );

    // ── Recommendation label ─────────────────────────
    const { text: recText, color: recColor } =
      combined >= 80 ? { text: '⭐ Excellent — use this gRNA',         color: '#00D4AA' }
    : combined >= 65 ? { text: '✓ Good — highly recommended',          color: '#4D9FFF' }
    : combined >= 50 ? { text: '△ Acceptable — verify off-targets',    color: '#F0C040' }
    : combined >= 35 ? { text: '▽ Poor — consider alternatives',       color: '#FF8C42' }
    :                  { text: '✗ Very poor — avoid',                  color: '#FF4D6D' };

    scored.push({
      ...cand,
      scores: {
        efficiency:   mlResult.score,
        capsCompat:   capsData.capsCompatibilityScore,
        offTarget:    offTargetScore,
        combined,
        breakdown:    mlResult.breakdown,
        ci:           mlResult.ci,
      },
      mlResult,
      capsData,
      offTarget,
      gcContent:      Math.round(
        (cand.sequence.match(/[GC]/g)||[]).length / cand.sequence.length * 100
      ),
      recommendation:      recText,
      recommendationColor: recColor,
      isPareto: false, // set by Pareto analysis below
    });
  }

  // Pareto analysis
  const withPareto = findParetoOptimal(scored);

  // Sort by combined score
  return withPareto.sort((a,b) => b.scores.combined - a.scores.combined);
}

// ══════════════════════════════════════════════════════════
// SECTION 5 — FULL DESIGN PIPELINE
// ══════════════════════════════════════════════════════════

/**
 * Complete gRNA + CAPS design pipeline
 * @param {string} sequence - target gene/amplicon
 * @param {Object} options
 * @returns {Object} complete design report
 */
export function designGRNAsWithCAPS(sequence, options = {}) {
  const {
    targetStart  = 0,
    targetEnd    = null,
    maxResults   = 20,
    effWeight    = 0.50,
    capsWeight   = 0.35,
    offTgtWeight = 0.15,
    species      = 'rice',
    minEfficiency = 30,
  } = options;

  if (!sequence || sequence.length < 30) {
    return {
      success: false,
      error: 'Sequence too short. Minimum 30bp required.',
      candidates: [],
    };
  }

  // Find all candidates
  const allCandidates = findAllGRNACandidates(sequence, targetStart, targetEnd);

  if (allCandidates.length === 0) {
    return {
      success: false,
      error: 'No valid NGG PAM sites found. Check sequence contains ATGC bases.',
      candidates: [],
    };
  }

  // Score and rank using ML
  const ranked = rankGRNACandidates(allCandidates, sequence, {
    efficiencyWeight: effWeight,
    capsWeight,
    offTargetWeight:  offTgtWeight,
    minEfficiency,
    species,
  });

  const top = ranked.slice(0, maxResults);

  // Segment by performance
  const excellent  = top.filter(c => c.scores.combined >= 75);
  const good       = top.filter(c => c.scores.combined >= 60 && c.scores.combined < 75);
  const withCAPS   = top.filter(c => c.capsData.hasCAPSCompatibility);
  const pareto     = top.filter(c => c.isPareto);

  // Best options
  const topOverall = top[0] || null;
  const bestCAPS   = [...top].sort((a,b) =>
    b.capsData.capsCompatibilityScore - a.capsData.capsCompatibilityScore)[0];
  const bestEff    = [...top].sort((a,b) =>
    b.scores.efficiency - a.scores.efficiency)[0];

  return {
    success:      true,
    species,
    totalFound:   allCandidates.length,
    totalScored:  ranked.length,
    candidates:   top,
    summary: {
      totalCandidates:  allCandidates.length,
      scoredCandidates: ranked.length,
      excellentGRNAs:   excellent.length,
      goodGRNAs:        good.length,
      withCAPSCompat:   withCAPS.length,
      paretoOptimal:    pareto.length,
      topGRNA:          topOverall,
      bestCAPSGRNA:     bestCAPS,
      bestEffGRNA:      bestEff,
    },
    recommendations: buildRecommendations(top, species),
    methodology: buildMethodsText(species),
  };
}

// ══════════════════════════════════════════════════════════
// SECTION 6 — RECOMMENDATIONS + METHODS TEXT
// ══════════════════════════════════════════════════════════

function buildRecommendations(ranked, species) {
  const recs = [];
  const top  = ranked[0];

  if (!top) return recs;

  // Primary recommendation
  recs.push({
    priority: 'high',
    icon:     '⭐',
    text: `Best overall gRNA: ${top.sequence} — Combined score ${top.scores.combined}% ` +
          `(Efficiency: ${top.scores.efficiency}%, CAPS: ${top.scores.capsCompat}%). ` +
          (top.capsData.hasCAPSCompatibility
            ? `Natural CAPS marker available with ${top.capsData.bestEnzyme}.`
            : `No natural CAPS — dCAPS design recommended after editing.`),
  });

  // CAPS-specific recommendation
  const bestCaps = ranked.find(c => c.capsData.capsCompatibilityScore >= 65);
  if (bestCaps && bestCaps !== top) {
    recs.push({
      priority: 'medium',
      icon:     '🔍',
      text: `Best CAPS compatibility: ${bestCaps.sequence.slice(0,10)}... ` +
            `with ${bestCaps.capsData.bestEnzyme} (CAPS score: ${bestCaps.capsData.capsCompatibilityScore}%). ` +
            `If easy T1 screening is the priority, consider this gRNA.`,
    });
  }

  // Thermodynamic note
  const topThermo = top.mlResult?.thermo;
  if (topThermo?.interpretation) {
    recs.push({
      priority: 'info',
      icon:     '🌡',
      text: `Thermodynamics: ${topThermo.interpretation.text} ` +
            `(ΔG = ${topThermo.fullDuplex?.dG} kcal/mol, seed ΔG = ${topThermo.seedRegion?.dG} kcal/mol).`,
    });
  }

  // Species-specific note
  const speciesNotes = {
    rice:        'Rice Cas9 editing works best at 25-28°C. Callus transformation efficiency is highest with high-GC gRNAs targeting coding exons.',
    wheat:       'Hexaploid wheat has 3 homeologous copies — verify gRNA specificity across A, B, D genomes. Low GC gRNAs may have higher homeolog specificity.',
    arabidopsis: 'Arabidopsis floral dip transformation is efficient. T1 seeds need kanamycin selection. gRNAs targeting early exons give more reliable knockouts.',
    maize:       'Maize requires biolistics or Agrobacterium. Target conserved exons for reliable knockout across inbred lines.',
    tomato:      'Tomato Agrobacterium transformation is standard. Target exon 1-3 for efficient knockout.',
  };
  if (speciesNotes[species]) {
    recs.push({
      priority: 'info',
      icon:     '🌾',
      text: speciesNotes[species],
    });
  }

  return recs;
}

/**
 * Generate methods text suitable for publication
 * This is a real contribution — helps researchers write their Methods section
 */
function buildMethodsText(species) {
  return `Guide RNA candidates were identified and scored using PlantaCRISPAI ` +
    `(v2.0), which implements a feature-engineered ensemble ML model for ` +
    `${species} CRISPR-Cas9 efficiency prediction. The model incorporates: ` +
    `(1) a position-weight matrix derived from Doench et al. (2016) with ` +
    `plant-specific corrections from Xu et al. (2015); ` +
    `(2) nearest-neighbor thermodynamic parameters for gRNA:DNA duplex stability ` +
    `(Sugimoto et al. 1995), including seed region (positions 1–12) and ` +
    `PAM-proximal (positions 1–8) free energy calculations; ` +
    `(3) structural penalty scores for self-complementarity, homopolymer runs, ` +
    `and RNA Pol III terminator sequences; and ` +
    `(4) ${species}-specific corrections for GC content optimum and CpG ` +
    `methylation context. CAPS marker compatibility was scored for each candidate ` +
    `based on restriction enzyme site proximity to the predicted cut site, ` +
    `band size distinguishability, and enzyme commercial reliability. ` +
    `Final ranking used a weighted ensemble: efficiency (50%), CAPS compatibility (35%), ` +
    `and off-target risk (15%). Pareto-optimal candidates were identified as those ` +
    `not dominated on both efficiency and CAPS dimensions simultaneously.`;
}
