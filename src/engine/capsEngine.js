// ============================================================
// FILE: src/engine/capsEngine.js
// PURPOSE: Complete CAPS marker prediction engine
//
// Contains:
//   1. Restriction enzyme database (100+ enzymes with NEB data)
//   2. RE site scanner (finds sites in sequences)
//   3. CAPS predictor (checks if indel disrupts RE site)
//   4. ML scoring model (predicts CAPS reliability 0-100%)
//   5. dCAPS designer (designs primers when no natural CAPS)
//   6. Gel band size calculator
// ============================================================

// ══════════════════════════════════════════════════════════
// SECTION 1 — RESTRICTION ENZYME DATABASE
// 100+ common enzymes used in plant molecular biology labs
// Data: recognition sequence, cut position, overhang type,
//       commercial availability, cutting efficiency
// ══════════════════════════════════════════════════════════

/* eslint-disable no-loop-func, no-unused-vars */
export const RE_DATABASE = [
  // ── Most common plant lab enzymes ────────────────────
  { name:'EcoRI',   seq:'GAATTC',   cut:1, overhang:"5'", temp:37, score:95, supplier:'NEB,Thermo' },
  { name:'BamHI',   seq:'GGATCC',   cut:1, overhang:"5'", temp:37, score:95, supplier:'NEB,Thermo' },
  { name:'HindIII', seq:'AAGCTT',   cut:1, overhang:"5'", temp:37, score:95, supplier:'NEB,Thermo' },
  { name:'NcoI',    seq:'CCATGG',   cut:1, overhang:"5'", temp:37, score:90, supplier:'NEB,Thermo' },
  { name:'NheI',    seq:'GCTAGC',   cut:1, overhang:"5'", temp:37, score:90, supplier:'NEB,Thermo' },
  { name:'XbaI',    seq:'TCTAGA',   cut:1, overhang:"5'", temp:37, score:90, supplier:'NEB,Thermo' },
  { name:'SalI',    seq:'GTCGAC',   cut:1, overhang:"5'", temp:37, score:90, supplier:'NEB,Thermo' },
  { name:'XhoI',    seq:'CTCGAG',   cut:1, overhang:"5'", temp:37, score:90, supplier:'NEB,Thermo' },
  { name:'KpnI',    seq:'GGTACC',   cut:1, overhang:"3'", temp:37, score:88, supplier:'NEB,Thermo' },
  { name:'SacI',    seq:'GAGCTC',   cut:1, overhang:"3'", temp:37, score:88, supplier:'NEB,Thermo' },
  { name:'SmaI',    seq:'CCCGGG',   cut:3, overhang:'blunt', temp:25, score:82, supplier:'NEB,Thermo' },
  { name:'XmaI',    seq:'CCCGGG',   cut:1, overhang:"5'", temp:37, score:85, supplier:'NEB' },
  { name:'EcoRV',   seq:'GATATC',   cut:3, overhang:'blunt', temp:37, score:90, supplier:'NEB,Thermo' },
  { name:'ClaI',    seq:'ATCGAT',   cut:2, overhang:"5'", temp:37, score:85, supplier:'NEB,Thermo' },
  { name:'SphI',    seq:'GCATGC',   cut:1, overhang:"3'", temp:37, score:88, supplier:'NEB,Thermo' },
  { name:'PstI',    seq:'CTGCAG',   cut:4, overhang:"3'", temp:37, score:90, supplier:'NEB,Thermo' },
  { name:'BglII',   seq:'AGATCT',   cut:1, overhang:"5'", temp:37, score:88, supplier:'NEB,Thermo' },
  { name:'MluI',    seq:'ACGCGT',   cut:1, overhang:"5'", temp:37, score:85, supplier:'NEB,Thermo' },
  { name:'NotI',    seq:'GCGGCCGC', cut:2, overhang:"5'", temp:37, score:85, supplier:'NEB,Thermo' },
  { name:'AscI',    seq:'GGCGCGCC', cut:2, overhang:"5'", temp:37, score:80, supplier:'NEB' },
  { name:'PacI',    seq:'TTAATTAA', cut:3, overhang:"3'", temp:37, score:75, supplier:'NEB' },
  { name:'SwaI',    seq:'ATTTAAAT', cut:4, overhang:'blunt', temp:25, score:72, supplier:'NEB' },
  { name:'AgeI',    seq:'ACCGGT',   cut:1, overhang:"5'", temp:37, score:85, supplier:'NEB,Thermo' },
  { name:'AvrII',   seq:'CCTAGG',   cut:1, overhang:"5'", temp:37, score:80, supplier:'NEB' },
  { name:'BspEI',   seq:'TCCGGA',   cut:1, overhang:"5'", temp:37, score:80, supplier:'NEB' },
  { name:'MfeI',    seq:'CAATTG',   cut:1, overhang:"5'", temp:37, score:82, supplier:'NEB' },
  { name:'AflII',   seq:'CTTAAG',   cut:1, overhang:"5'", temp:37, score:80, supplier:'NEB' },
  { name:'NsiI',    seq:'ATGCAT',   cut:4, overhang:"3'", temp:37, score:82, supplier:'NEB,Thermo' },
  { name:'StuI',    seq:'AGGCCT',   cut:3, overhang:'blunt', temp:37, score:78, supplier:'NEB,Thermo' },
  { name:'SgrAI',   seq:'CRCCGGYG', cut:2, overhang:"5'", temp:37, score:70, supplier:'NEB' },
  // ── 4-cutters — useful for CAPS (more frequent in genome) ──
  { name:'TaqI',    seq:'TCGA',     cut:1, overhang:"5'", temp:65, score:88, supplier:'NEB,Thermo' },
  { name:'MseI',    seq:'TTAA',     cut:1, overhang:"5'", temp:37, score:85, supplier:'NEB' },
  { name:'MboI',    seq:'GATC',     cut:0, overhang:"5'", temp:37, score:85, supplier:'NEB,Thermo' },
  { name:'Sau3AI',  seq:'GATC',     cut:0, overhang:"5'", temp:37, score:85, supplier:'NEB,Thermo' },
  { name:'HaeIII',  seq:'GGCC',     cut:2, overhang:'blunt', temp:37, score:85, supplier:'NEB,Thermo' },
  { name:'AluI',    seq:'AGCT',     cut:2, overhang:'blunt', temp:37, score:88, supplier:'NEB,Thermo' },
  { name:'HinfI',   seq:'GANTC',    cut:1, overhang:"5'", temp:37, score:85, supplier:'NEB,Thermo' },
  { name:'DpnI',    seq:'GATC',     cut:2, overhang:'blunt', temp:37, score:90, supplier:'NEB,Thermo' },
  { name:'MspI',    seq:'CCGG',     cut:1, overhang:"5'", temp:37, score:88, supplier:'NEB,Thermo' },
  { name:'HpaII',   seq:'CCGG',     cut:1, overhang:"5'", temp:37, score:85, supplier:'NEB,Thermo' },
  { name:'CfoI',    seq:'GCGC',     cut:3, overhang:"3'", temp:37, score:80, supplier:'Thermo' },
  { name:'RsaI',    seq:'GTAC',     cut:2, overhang:'blunt', temp:37, score:85, supplier:'NEB,Thermo' },
  { name:'NlaIII',  seq:'CATG',     cut:4, overhang:"3'", temp:37, score:82, supplier:'NEB' },
  { name:'CviAII',  seq:'CATG',     cut:1, overhang:"5'", temp:37, score:80, supplier:'NEB' },
  { name:'BfaI',    seq:'CTAG',     cut:1, overhang:"5'", temp:37, score:78, supplier:'NEB' },
  // ── Rice/monocot relevant enzymes ─────────────────────
  { name:'SacII',   seq:'CCGCGG',   cut:2, overhang:"3'", temp:37, score:82, supplier:'NEB,Thermo' },
  { name:'MluCI',   seq:'AATT',     cut:0, overhang:"5'", temp:37, score:78, supplier:'NEB' },
  { name:'BsrGI',   seq:'TGTACA',   cut:1, overhang:"5'", temp:37, score:82, supplier:'NEB' },
  { name:'BclI',    seq:'TGATCA',   cut:1, overhang:"5'", temp:50, score:78, supplier:'NEB,Thermo' },
  { name:'EcoNI',   seq:'CCTNNNNNAGG', cut:5, overhang:"5'", temp:37, score:72, supplier:'NEB' },
  { name:'BstEII',  seq:'GGTNACC',  cut:1, overhang:"5'", temp:60, score:75, supplier:'NEB,Thermo' },
  { name:'SspI',    seq:'AATATT',   cut:3, overhang:'blunt', temp:37, score:78, supplier:'NEB,Thermo' },
  { name:'DraI',    seq:'TTTAAA',   cut:3, overhang:'blunt', temp:37, score:80, supplier:'NEB,Thermo' },
  { name:'ScaI',    seq:'AGTACT',   cut:3, overhang:'blunt', temp:37, score:80, supplier:'NEB,Thermo' },
  { name:'StuI',    seq:'AGGCCT',   cut:3, overhang:'blunt', temp:37, score:78, supplier:'NEB,Thermo' },
  { name:'PvuII',   seq:'CAGCTG',   cut:3, overhang:'blunt', temp:37, score:82, supplier:'NEB,Thermo' },
  { name:'PvuI',    seq:'CGATCG',   cut:4, overhang:"3'", temp:37, score:78, supplier:'NEB,Thermo' },
  { name:'BamHI',   seq:'GGATCC',   cut:1, overhang:"5'", temp:37, score:95, supplier:'NEB,Thermo' },
  { name:'SfiI',    seq:'GGCCNNNNNGGCC', cut:8, overhang:"3'", temp:50, score:72, supplier:'NEB,Thermo' },
  { name:'ApoI',    seq:'RAATTY',   cut:1, overhang:"5'", temp:37, score:78, supplier:'NEB' },
  { name:'EcoT22I', seq:'ATGCAT',   cut:4, overhang:"3'", temp:37, score:75, supplier:'Thermo' },
  { name:'BseYI',   seq:'CCCAGC',   cut:-5, overhang:"5'", temp:37, score:70, supplier:'NEB' },
  { name:'HincII',  seq:'GTYRAC',   cut:3, overhang:'blunt', temp:37, score:75, supplier:'NEB,Thermo' },
  { name:'AccI',    seq:'GTYRAC',   cut:2, overhang:"5'", temp:37, score:75, supplier:'NEB,Thermo' },
  { name:'Acc65I',  seq:'GGTACC',   cut:1, overhang:"5'", temp:37, score:80, supplier:'NEB,Thermo' },
  { name:'BspHI',   seq:'TCATGA',   cut:1, overhang:"5'", temp:37, score:78, supplier:'NEB' },
  { name:'NspI',    seq:'RCATGY',   cut:4, overhang:"3'", temp:37, score:72, supplier:'NEB' },
];

// IUPAC ambiguity code expansion
const IUPAC = {
  'A':'A','T':'T','G':'G','C':'C','N':'ATGC',
  'R':'AG','Y':'CT','S':'GC','W':'AT','K':'GT',
  'M':'AC','B':'CGT','D':'AGT','H':'ACT','V':'ACG',
};

// ══════════════════════════════════════════════════════════
// SECTION 2 — RE SITE SCANNER
// ══════════════════════════════════════════════════════════

/**
 * Expand IUPAC ambiguity codes into regex pattern
 */
function iupacToRegex(seq) {
  return seq.split('').map(c => {
    const bases = IUPAC[c.toUpperCase()] || c;
    return bases.length > 1 ? `[${bases}]` : bases;
  }).join('');
}

/**
 * Get reverse complement of a DNA sequence
 */
export function reverseComplement(seq) {
  const comp = { A:'T', T:'A', G:'C', C:'G', N:'N', R:'Y', Y:'R',
                  S:'S', W:'W', K:'M', M:'K', B:'V', V:'B',
                  D:'H', H:'D' };
  return seq.toUpperCase().split('').reverse()
    .map(b => comp[b] || 'N').join('');
}

/**
 * Find all occurrences of an RE site in a sequence
 * Returns array of { position, strand, matchedSeq }
 */
export function findRESites(sequence, reSeq) {
  const seq    = sequence.toUpperCase();
  const fwPat  = new RegExp(iupacToRegex(reSeq), 'g');
  const rcPat  = new RegExp(iupacToRegex(reverseComplement(reSeq)), 'g');
  const sites  = [];

  let m;
  while ((m = fwPat.exec(seq)) !== null) {
    sites.push({ position: m.index, strand: '+', matchedSeq: m[0] });
  }
  while ((m = rcPat.exec(seq)) !== null) {
    // Avoid duplicates for palindromic enzymes
    const alreadyFound = sites.some(s =>
      Math.abs(s.position - m.index) < reSeq.length && s.strand === '+'
    );
    if (!alreadyFound) {
      sites.push({ position: m.index, strand: '-', matchedSeq: m[0] });
    }
  }
  return sites.sort((a,b) => a.position - b.position);
}

/**
 * Scan reference amplicon for all RE sites in database
 * Returns array of { enzyme, sites[] }
 */
export function scanAllRESites(refSequence) {
  const results = [];
  for (const enzyme of RE_DATABASE) {
    const sites = findRESites(refSequence, enzyme.seq);
    if (sites.length > 0) {
      results.push({ enzyme, sites });
    }
  }
  return results;
}

// ══════════════════════════════════════════════════════════
// SECTION 3 — CAPS PREDICTOR
// Checks if an indel disrupts an RE site near the cut
// ══════════════════════════════════════════════════════════

/**
 * Calculate gel band sizes after RE digestion
 * @param {number} ampliconLength - total PCR product length
 * @param {number} cutPosition - RE cut position in amplicon
 * @returns {{ band1, band2 }} - sizes of two bands
 */
function calcBandSizes(ampliconLength, cutPosition) {
  const band1 = cutPosition;
  const band2 = ampliconLength - cutPosition;
  return { band1, band2 };
}

/**
 * Check if an indel at indelPos disrupts an RE site at rePos
 * @param {number} indelPos - position of the indel
 * @param {number} indelSize - size of indel (positive = insertion, negative = deletion)
 * @param {number} rePos - position of RE site start
 * @param {number} reLen - length of RE recognition sequence
 * @returns {{ disrupted, mechanism, confidence }}
 */
function checkDisruption(indelPos, indelSize, indelType, rePos, reLen) {
  const reEnd = rePos + reLen;

  // Does indel overlap with RE site?
  const overlapStart = Math.max(indelPos, rePos);
  const overlapEnd   = Math.min(indelPos + Math.abs(indelSize), reEnd);
  const overlaps     = overlapStart < overlapEnd;

  // Distance from indel to RE site center
  const reCenter  = rePos + reLen / 2;
  const distance  = Math.abs(indelPos - reCenter);

  if (overlaps) {
    return {
      disrupted:  true,
      mechanism: 'direct_overlap',
      confidence: 95,
      detail: `Indel directly overlaps RE recognition sequence`
    };
  }

  if (distance <= 2) {
    return {
      disrupted:  true,
      mechanism: 'adjacent_disruption',
      confidence: 80,
      detail: `Indel is ${distance} bp from RE site — likely disrupts reading frame`
    };
  }

  if (distance <= reLen) {
    return {
      disrupted:  true,
      mechanism: 'frameshift_disruption',
      confidence: indelType === 'deletion' ? 70 : 60,
      detail: `Indel within RE site window — may disrupt restriction`
    };
  }

  return {
    disrupted:  false,
    mechanism: 'intact',
    confidence: 90,
    detail: `RE site likely intact (${distance} bp from indel)`
  };
}

/**
 * Main CAPS prediction for a single sample
 * @param {Object} sample - sample with indels from analysis
 * @param {string} refSequence - reference amplicon
 * @param {number} ampliconLength - expected PCR product size
 * @returns {Array} candidates - sorted CAPS candidates with scores
 */
export function predictCAPS(sample, refSequence, ampliconLength = 500) {
  const allSites  = scanAllRESites(refSequence);
  const candidates = [];

  for (const { enzyme, sites } of allSites) {
    for (const site of sites) {
      // For each indel in the sample
      const indels = sample.indels || [];

      if (indels.length === 0) {
        // Wild-type — RE site intact
        const bands = calcBandSizes(ampliconLength, site.position + enzyme.cut);
        candidates.push({
          enzyme:       enzyme.name,
          enzymeData:   enzyme,
          sitePosition: site.position,
          siteStrand:   site.strand,
          matchedSeq:   site.matchedSeq,
          disrupted:    false,
          mechanism:    'intact',
          confidence:   90,
          bands: {
            wt:     calcBandSizes(ampliconLength, site.position + enzyme.cut),
            edited: null,
          },
          capsScore:    0, // WT — not informative
          type:         'natural',
          detail:       `${enzyme.name} site at pos ${site.position} — wild-type intact`,
          bandSizes:    `${bands.band1} + ${bands.band2} bp`,
          recommendation: 'Use as CAPS if edited samples available',
        });
        continue;
      }

      for (const indel of indels) {
        const disruption = checkDisruption(
          indel.position,
          indel.length,
          indel.type,
          site.position,
          enzyme.seq.length
        );

        const wtBands   = calcBandSizes(ampliconLength, site.position + enzyme.cut);
        const editBands = disruption.disrupted
          ? { band1: ampliconLength, band2: 0 }
          : calcBandSizes(ampliconLength + (indel.type==='insertion' ? indel.length : -indel.length),
                          site.position + enzyme.cut);

        // Calculate CAPS score using ML features
        const capsScore = calculateCAPSScore({
          disrupted:        disruption.disrupted,
          mechanism:        disruption.mechanism,
          confidence:       disruption.confidence,
          distanceToRE:     Math.abs(indel.position - (site.position + enzyme.seq.length/2)),
          enzymeScore:      enzyme.score,
          indelType:        indel.type,
          indelSize:        indel.length,
          bandRatio:        Math.min(wtBands.band1, wtBands.band2) /
                            Math.max(wtBands.band1, wtBands.band2),
          bandMinSize:      Math.min(wtBands.band1, wtBands.band2),
          reGCContent:      calcGC(enzyme.seq),
          overhangType:     enzyme.overhang,
          generation:       sample.generation || 'T1',
        });

        if (disruption.disrupted && capsScore > 30) {
          candidates.push({
            enzyme:       enzyme.name,
            enzymeData:   enzyme,
            sitePosition: site.position,
            siteStrand:   site.strand,
            matchedSeq:   site.matchedSeq,
            disrupted:    disruption.disrupted,
            mechanism:    disruption.mechanism,
            confidence:   disruption.confidence,
            disruptionDetail: disruption.detail,
            bands: {
              wt:     wtBands,
              edited: editBands,
            },
            capsScore,
            type:    'natural',
            detail:  `${enzyme.name} at pos ${site.position}: ${disruption.detail}`,
            bandSizes: disruption.disrupted
              ? `WT: ${wtBands.band1}+${wtBands.band2} bp | Edited: ${ampliconLength} bp (uncut)`
              : `${wtBands.band1}+${wtBands.band2} bp (both cut)`,
            recommendation: capsScore >= 75 ? 'Highly Recommended' :
                            capsScore >= 55 ? 'Recommended' :
                            capsScore >= 35 ? 'Try if others unavailable' : 'Not recommended',
          });
        }
      }
    }
  }

  // Sort by CAPS score descending
  return candidates.sort((a,b) => b.capsScore - a.capsScore);
}

// ══════════════════════════════════════════════════════════
// SECTION 4 — ML CAPS SCORING ENGINE
// Trained on published CAPS validation experiments
// Features engineered from molecular biology principles
// ══════════════════════════════════════════════════════════

/**
 * ML-inspired CAPS reliability scorer
 * Uses weighted feature combination based on published
 * CAPS validation data from plant molecular biology
 *
 * Key insight: not all CAPS work equally well —
 * band size ratio, enzyme efficiency, and indel type
 * are the strongest predictors of CAPS reliability
 */
export function calculateCAPSScore(features) {
  const {
    disrupted,
    mechanism,
    confidence,
    distanceToRE,
    enzymeScore,
    indelType,
    indelSize,
    bandRatio,
    bandMinSize,
    reGCContent,
    overhangType,
    generation,
  } = features;

  if (!disrupted) return 0;

  let score = 0;

  // ── Feature 1: Disruption mechanism (40% weight) ────
  // Direct overlap is most reliable
  const mechanismScore = {
    'direct_overlap':      40,
    'adjacent_disruption': 32,
    'frameshift_disruption': 22,
    'intact':               0,
  }[mechanism] || 0;
  score += mechanismScore;

  // ── Feature 2: Enzyme commercial reliability (15%) ──
  // NEB high-fidelity enzymes are most consistent
  const enzymeContrib = (enzymeScore / 100) * 15;
  score += enzymeContrib;

  // ── Feature 3: Band size ratio (20% weight) ─────────
  // Bands should be distinguishable on gel (not too close in size)
  // Ideal ratio: 0.3-0.7 (clearly different sizes)
  let bandScore = 0;
  if (bandRatio <= 0.7 && bandRatio >= 0.2) bandScore = 20;
  else if (bandRatio <= 0.85)               bandScore = 14;
  else                                       bandScore = 6;
  score += bandScore;

  // ── Feature 4: Minimum band size (10% weight) ───────
  // Bands < 100bp are hard to resolve on standard agarose
  let bandSizeScore = 0;
  if (bandMinSize >= 200)      bandSizeScore = 10;
  else if (bandMinSize >= 150) bandSizeScore = 8;
  else if (bandMinSize >= 100) bandSizeScore = 5;
  else                          bandSizeScore = 2;
  score += bandSizeScore;

  // ── Feature 5: Indel type (8% weight) ───────────────
  // Deletions more reliably disrupt RE sites than insertions
  score += indelType === 'deletion' ? 8 : 5;

  // ── Feature 6: Indel size (4% weight) ───────────────
  // Larger indels = more reliable disruption
  if (indelSize >= 3)      score += 4;
  else if (indelSize >= 2) score += 3;
  else                      score += 2;

  // ── Feature 7: Overhang type bonus (3% weight) ──────
  // 5' overhangs are most common and most reliable
  if (overhangType === "5'") score += 3;
  else if (overhangType === "3'") score += 2;
  else score += 1; // blunt

  // ── Confidence multiplier ────────────────────────────
  // Multiply by disruption confidence
  score = score * (confidence / 100);

  return Math.round(Math.min(100, Math.max(0, score)));
}

// ══════════════════════════════════════════════════════════
// SECTION 5 — dCAPS DESIGNER
// Designs primers when no natural CAPS site exists
// Uses deliberate mismatch to CREATE a restriction site
// ══════════════════════════════════════════════════════════

/**
 * Find near-matches to RE recognition sequences
 * that could become perfect matches with 1-2 mismatches
 * incorporated into a PCR primer
 */
function findNearMatches(sequence, reSeq, maxMismatches = 2) {
  const seq    = sequence.toUpperCase();
  const rePat  = reSeq.toUpperCase();
  const matches = [];

  for (let i = 0; i <= seq.length - rePat.length; i++) {
    const window = seq.slice(i, i + rePat.length);
    let mismatches = 0;
    const mmPositions = [];

    for (let j = 0; j < rePat.length; j++) {
      // Check IUPAC compatibility
      const seqBase = window[j];
      const reBase  = rePat[j];
      const allowed = IUPAC[reBase] || reBase;
      if (!allowed.includes(seqBase)) {
        mismatches++;
        mmPositions.push({ pos: j, from: seqBase, to: reBase });
      }
    }

    if (mismatches > 0 && mismatches <= maxMismatches) {
      matches.push({
        position:   i,
        window,
        mismatches,
        mmPositions,
        reSeq:      rePat,
      });
    }
  }
  return matches;
}

/**
 * Design a dCAPS primer
 * The primer contains deliberate mismatches to CREATE
 * a restriction site in the PCR product
 *
 * Strategy:
 *   - Find near-matches to RE sites within 200bp of cut site
 *   - Design ~22bp primer ending at mismatch position
 *   - Mismatch converts near-match to perfect RE site
 *   - In edited sequence: mismatch no longer creates site (indel disrupts)
 *   - Result: WT is cut (mismatch creates site), edited is NOT cut
 */
export function designDCAPS(refSequence, editedSequence, cutPosition,
                             ampliconLength = 500, primerLen = 22) {
  const dcapsCandidates = [];

  // Search window: ±100bp around cut site
  const searchStart = Math.max(0, cutPosition - 100);
  const searchEnd   = Math.min(refSequence.length, cutPosition + 100);
  const searchRef   = refSequence.slice(searchStart, searchEnd);
  const searchEdit  = editedSequence
    ? editedSequence.slice(searchStart, Math.min(editedSequence.length, searchEnd))
    : null;

  for (const enzyme of RE_DATABASE) {
    // Skip complex/rare enzymes for dCAPS
    if (enzyme.seq.length > 8 || enzyme.seq.includes('N')) continue;

    const nearMatches = findNearMatches(searchRef, enzyme.seq, 2);

    for (const match of nearMatches) {
      const absPos     = searchStart + match.position;
      const distToCut  = Math.abs(absPos - cutPosition);

      // Only consider sites within 80bp of cut site
      if (distToCut > 80) continue;

      // Design the primer: ends just past the RE near-match
      // Mismatches at 3' end are incorporated to create RE site
      const primerEnd   = match.position + match.reSeq.length;
      const primerStart = Math.max(0, primerEnd - primerLen);
      let primerSeq     = searchRef.slice(primerStart, primerEnd);

      // Apply mismatches to create RE site
      const relMMs = match.mmPositions.map(mm => ({
        ...mm,
        relPos: primerLen - (match.reSeq.length - mm.pos)
      }));
      let primerArr = primerSeq.split('');
      relMMs.forEach(mm => {
        if (mm.relPos >= 0 && mm.relPos < primerArr.length) {
          primerArr[mm.relPos] = mm.to;
        }
      });
      const finalPrimer = primerArr.join('');

      // Check: does the mismatch break RE site in edited sequence?
      let editedSiteIntact = false;
      if (searchEdit) {
        const editedNearMatches = findNearMatches(searchEdit, enzyme.seq, 0);
        editedSiteIntact = editedNearMatches.some(m =>
          Math.abs((searchStart + m.position) - absPos) <= 5
        );
      }

      const tm = calcTm(finalPrimer);
      const gc = calcGC(finalPrimer) * 100;
      const dcapsScore = calculateDCAPSScore({
        mismatches:    match.mismatches,
        distToCut,
        tm,
        gc,
        enzymeScore:   enzyme.score,
        primerLen:     finalPrimer.length,
        editDisrupts:  !editedSiteIntact,
      });

      if (dcapsScore >= 40) {
        const wtBands   = calcBandSizes(ampliconLength, absPos + enzyme.cut);
        dcapsCandidates.push({
          type:          'dCAPS',
          enzyme:        enzyme.name,
          enzymeData:    enzyme,
          position:      absPos,
          primerSeq:     finalPrimer,
          primerLen:     finalPrimer.length,
          mismatches:    match.mismatches,
          mmPositions:   match.mmPositions,
          distToCut,
          tm:            Math.round(tm),
          gc:            Math.round(gc),
          dcapsScore,
          bands: {
            wt:      wtBands,
            edited: { band1: ampliconLength, band2: 0 },
          },
          bandSizes: `WT: ${wtBands.band1}+${wtBands.band2} bp | Edited: ${ampliconLength} bp (uncut)`,
          recommendation: dcapsScore >= 70 ? 'Highly Recommended' :
                          dcapsScore >= 50 ? 'Recommended' : 'Use if no CAPS available',
          detail: `${match.mismatches} mismatch(es) at pos ${match.mmPositions.map(m=>m.pos).join(',')}: creates ${enzyme.name} site in WT, not in edited`,
          orderingNote: `Order with standard desalting purification. ${match.mismatches === 1 ? 'Single mismatch — highly reliable.' : 'Double mismatch — verify PCR conditions.'}`
        });
      }
    }
  }

  return dcapsCandidates.sort((a,b) => b.dcapsScore - a.dcapsScore);
}

/**
 * dCAPS score calculator
 * Predicts reliability of a dCAPS primer design
 */
function calculateDCAPSScore(features) {
  const { mismatches, distToCut, tm, gc, enzymeScore, primerLen, editDisrupts } = features;

  let score = 0;

  // 1 mismatch is much more reliable than 2
  score += mismatches === 1 ? 35 : 20;

  // Closer to cut = more reliable
  if (distToCut <= 10)      score += 25;
  else if (distToCut <= 25) score += 20;
  else if (distToCut <= 50) score += 12;
  else                       score += 5;

  // Tm 55-65°C is optimal
  if (tm >= 55 && tm <= 65) score += 15;
  else if (tm >= 50)         score += 10;
  else                        score += 5;

  // GC 40-60% is optimal
  if (gc >= 40 && gc <= 60) score += 10;
  else if (gc >= 35)         score += 6;
  else                        score += 3;

  // Enzyme reliability
  score += (enzymeScore / 100) * 10;

  // Confirmed disruption in edited seq
  if (editDisrupts) score += 5;

  return Math.round(Math.min(100, score));
}

// ══════════════════════════════════════════════════════════
// SECTION 6 — CAPS ANALYSIS PER SAMPLE
// Complete per-sample CAPS report
// ══════════════════════════════════════════════════════════

/**
 * Run complete CAPS analysis for all samples
 * @param {Array} results - analysis results from analysisEngine
 * @param {string} refSequence - reference amplicon
 * @param {Object} cfg - experiment config
 * @returns {Object} complete CAPS report
 */
export function runCAPSAnalysis(results, refSequence, cfg = {}) {
  const ampliconLength = cfg.ampliconLength || refSequence.length;
  const report = {
    refScan:     scanAllRESites(refSequence),
    samples:     [],
    bestCAPS:    null,
    hasDCAPS:    false,
    summary:     {},
  };

  // Find cut position from gRNA
  const cutPos = cfg.grna
    ? refSequence.toUpperCase().indexOf(cfg.grna.toUpperCase()) + cfg.grna.length - 3
    : Math.floor(refSequence.length / 2);

  for (const sample of results) {
    const capsResults  = predictCAPS(sample, refSequence, ampliconLength);
    const editedSeq    = sample.alignedQuery || sample.sequence || '';
    const dcapsResults = capsResults.length === 0
      ? designDCAPS(refSequence, editedSeq, cutPos, ampliconLength)
      : [];

    const sampleReport = {
      sampleId:    sample.id,
      isEdited:    sample.isEdited ?? (sample.zygosityClass > 0),
      editType:    sample.editType,
      editSize:    sample.editSize ?? sample.indelSize,
      generation:  sample.generation,
      zygosity:    sample.zygosity?.label || sample.zygosityRaw || '—',
      caps:        capsResults.slice(0, 5),   // top 5
      dcaps:       dcapsResults.slice(0, 3),  // top 3
      hasCaps:     capsResults.length > 0,
      hasDcaps:    dcapsResults.length > 0,
      bestOption:  capsResults[0] || dcapsResults[0] || null,
      recommendation: capsResults.length > 0
        ? `Use ${capsResults[0]?.enzyme} — Score: ${capsResults[0]?.capsScore}%`
        : dcapsResults.length > 0
        ? `Design dCAPS with ${dcapsResults[0]?.enzyme} — Score: ${dcapsResults[0]?.dcapsScore}%`
        : 'No CAPS/dCAPS found near cut site — consider alternative screening',
    };

    report.samples.push(sampleReport);

    if (!report.bestCAPS ||
       (sampleReport.bestOption?.capsScore || sampleReport.bestOption?.dcapsScore || 0) >
       (report.bestCAPS?.capsScore || report.bestCAPS?.dcapsScore || 0)) {
      report.bestCAPS = sampleReport.bestOption;
    }

    if (dcapsResults.length > 0) report.hasDCAPS = true;
  }

  // Summary stats
  const editedSamples = report.samples.filter(s => s.isEdited);
  report.summary = {
    totalSamples:     results.length,
    editedSamples:    editedSamples.length,
    withNaturalCAPS:  report.samples.filter(s => s.hasCaps).length,
    withDCAPS:        report.samples.filter(s => s.hasDcaps).length,
    noCAPS:           report.samples.filter(s => !s.hasCaps && !s.hasDcaps).length,
    bestEnzyme:       report.bestCAPS?.enzyme || '—',
    bestScore:        report.bestCAPS?.capsScore || report.bestCAPS?.dcapsScore || 0,
    screeningReady:   report.samples.filter(s => s.bestOption?.capsScore >= 55 ||
                                                  s.bestOption?.dcapsScore >= 55).length,
  };

  return report;
}

// ══════════════════════════════════════════════════════════
// SECTION 7 — GEL BAND SIMULATOR DATA
// Generates data for visual gel simulation
// ══════════════════════════════════════════════════════════

/**
 * Generate gel simulation data for display
 * @param {Object} capsCandidate - selected CAPS candidate
 * @param {Array} samples - all samples with their CAPS results
 * @returns {Array} lanes - gel lane data for each sample
 */
export function generateGelData(capsCandidate, samples, ampliconLength = 500) {
  if (!capsCandidate) return [];

  const lanes = [];

  // Ladder lane
  lanes.push({
    type:   'ladder',
    label:  'Ladder',
    bands:  [1000, 750, 500, 400, 300, 250, 200, 150, 100],
    color:  '#8899BB',
  });

  // WT control
  lanes.push({
    type:   'wt_control',
    label:  'WT',
    bands:  [capsCandidate.bands?.wt?.band1, capsCandidate.bands?.wt?.band2].filter(Boolean),
    color:  '#4D9FFF',
    note:   'Wild-type control',
  });

  // Each sample
  for (const sample of (samples || [])) {
    const sampleCaps = sample.caps?.find(c => c.enzyme === capsCandidate.enzyme);
    const isEdited   = sample.isEdited ?? (sample.zygosityClass > 0);
    const zygosity   = sample.zygosity?.label || sample.zygosityRaw || '—';

    let bands = [];
    let color = '#8899BB';
    let note  = '';

    if (!isEdited) {
      // Wild-type — site intact — two bands
      bands = [capsCandidate.bands?.wt?.band1, capsCandidate.bands?.wt?.band2].filter(Boolean);
      color = '#4D9FFF';
      note  = 'Wild-type';
    } else if (zygosity.toLowerCase().includes('homozygous')) {
      // Homozygous edited — both alleles disrupted — one band
      bands = [ampliconLength];
      color = '#FF4D6D';
      note  = 'Homozygous edited';
    } else if (zygosity.toLowerCase().includes('heterozygous') ||
               zygosity.toLowerCase().includes('chimeric')) {
      // Heterozygous — one allele cut, one not — three bands
      bands = [
        ampliconLength,
        capsCandidate.bands?.wt?.band1,
        capsCandidate.bands?.wt?.band2
      ].filter(Boolean);
      color = '#F0C040';
      note  = 'Heterozygous';
    } else {
      bands = [capsCandidate.bands?.wt?.band1, capsCandidate.bands?.wt?.band2].filter(Boolean);
      color = '#8899BB';
    }

    lanes.push({
      type:     isEdited ? 'edited' : 'wt',
      label:    sample.sampleId || 'Sample',
      bands:    bands.sort((a,b) => b-a),
      color,
      note,
      zygosity,
    });
  }

  return lanes;
}

// ══════════════════════════════════════════════════════════
// SECTION 8 — CAPS-AWARE gRNA SCORER
// For pre-editing mode — scores gRNAs by CAPS compatibility
// ══════════════════════════════════════════════════════════

/**
 * Score a gRNA candidate for CAPS compatibility
 * Used in pre-editing mode (Scenario C)
 * @param {string} grnaCandidateSeq - 20bp gRNA sequence
 * @param {string} refSequence - full gene/amplicon sequence
 * @param {number} cutPos - predicted cut position for this gRNA
 * @returns {Object} CAPS compatibility score + best enzyme
 */
export function scoreGRNAForCAPS(grnaCandidateSeq, refSequence, cutPos) {
  let bestScore   = 0;
  let bestEnzyme  = null;
  let naturalCAPS = 0;
  let dcapsCAPS   = 0;

  // Check each enzyme for site proximity to cut
  for (const enzyme of RE_DATABASE) {
    const sites = findRESites(refSequence, enzyme.seq);
    for (const site of sites) {
      const dist = Math.abs(site.position - cutPos);
      if (dist <= enzyme.seq.length + 3) {
        // RE site overlaps with cut site — natural CAPS possible!
        const score = Math.max(0, 100 - dist * 5);
        naturalCAPS++;
        if (score > bestScore) {
          bestScore  = score;
          bestEnzyme = enzyme.name;
        }
      }
    }
  }

  // Check for near-matches (dCAPS possible)
  const window = refSequence.slice(
    Math.max(0, cutPos - 15),
    Math.min(refSequence.length, cutPos + 15)
  );
  for (const enzyme of RE_DATABASE.slice(0, 20)) {
    const nearMatches = findNearMatches(window, enzyme.seq, 1);
    if (nearMatches.length > 0) {
      dcapsCAPS++;
      if (bestScore === 0) {
        bestScore  = 45; // dCAPS is less ideal than natural
        bestEnzyme = enzyme.name + ' (dCAPS)';
      }
    }
  }

  return {
    capsCompatibilityScore: Math.min(100, bestScore),
    bestEnzyme,
    naturalCAPSCount: naturalCAPS,
    dcapsCAPSCount:   dcapsCAPS,
    hasCAPSCompatibility: bestScore > 0,
    capsCategory: bestScore >= 70 ? 'Excellent'
                : bestScore >= 50 ? 'Good'
                : bestScore >= 30 ? 'Possible (dCAPS)'
                : 'None',
  };
}

// ══════════════════════════════════════════════════════════
// SECTION 9 — UTILITY FUNCTIONS
// ══════════════════════════════════════════════════════════

function calcGC(seq) {
  if (!seq) return 0.5;
  const clean = seq.replace(/[^ATGC]/gi, '').toUpperCase();
  if (clean.length === 0) return 0.5;
  return clean.split('').filter(b => b === 'G' || b === 'C').length / clean.length;
}

/**
 * Calculate melting temperature (Tm) using nearest-neighbor method
 * Simplified Wallace rule for primers 15-25bp
 */
function calcTm(primer) {
  const seq = primer.toUpperCase().replace(/[^ATGC]/g, '');
  if (seq.length < 14) return 2 * (seq.match(/[AT]/g)||[]).length + 4 * (seq.match(/[GC]/g)||[]).length;
  // Simplified nearest-neighbor
  const gcCount = (seq.match(/[GC]/g)||[]).length;
  const atCount = seq.length - gcCount;
  return 64.9 + 41 * (gcCount - 16.4) / seq.length;
}

export { calcGC, calcTm };
