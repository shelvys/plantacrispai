/* eslint-disable no-unused-vars */
// ============================================================
// FILE: src/engine/thermodynamics.js
// PURPOSE: Thermodynamic calculations for gRNA efficiency
//
// Implements nearest-neighbor thermodynamic model for:
//   1. gRNA:DNA duplex free energy (ΔG)
//   2. gRNA secondary structure (self-folding ΔG)
//   3. Seed region stability
//   4. PAM-proximal melting temperature
//
// Based on:
//   Sugimoto et al. 1995 (RNA:DNA nearest-neighbor parameters)
//   Freier et al. 1986 (DNA:DNA parameters)
//   SantaLucia 1998 (unified nearest-neighbor parameters)
// ============================================================

// ── RNA:DNA nearest-neighbor parameters (ΔH kcal/mol, ΔS cal/mol·K) ──
// From Sugimoto et al. 1995 — used for gRNA:DNA hybrid stability
// Key for CRISPR: gRNA is RNA, target is DNA
const RNA_DNA_NN = {
  'AA': { dH: -7.8,  dS: -21.9 },
  'AC': { dH: -5.9,  dS: -12.3 },
  'AG': { dH: -9.1,  dS: -23.5 },
  'AT': { dH: -8.3,  dS: -23.9 },
  'CA': { dH: -9.0,  dS: -26.1 },
  'CC': { dH: -9.3,  dS: -23.2 },
  'CG': { dH: -16.3, dS: -47.1 },
  'CT': { dH: -7.0,  dS: -19.7 },
  'GA': { dH: -5.5,  dS: -13.5 },
  'GC': { dH: -8.0,  dS: -17.1 },
  'GG': { dH: -12.8, dS: -31.9 },
  'GT': { dH: -7.8,  dS: -21.6 },
  'TA': { dH: -7.8,  dS: -23.2 },
  'TC': { dH: -5.5,  dS: -13.5 },
  'TG': { dH: -9.4,  dS: -26.7 },
  'TT': { dH: -11.5, dS: -36.4 },
};

// ── DNA:DNA nearest-neighbor parameters ──
// SantaLucia 1998 — used for primer Tm calculation
const DNA_DNA_NN = {
  'AA': { dH: -7.9,  dS: -22.2 },
  'AC': { dH: -7.8,  dS: -21.0 },
  'AG': { dH: -7.2,  dS: -20.4 },
  'AT': { dH: -7.2,  dS: -21.3 },
  'CA': { dH: -8.5,  dS: -22.7 },
  'CC': { dH: -8.0,  dS: -19.9 },
  'CG': { dH: -10.6, dS: -27.2 },
  'CT': { dH: -7.8,  dS: -21.0 },
  'GA': { dH: -8.2,  dS: -22.2 },
  'GC': { dH: -9.8,  dS: -24.4 },
  'GG': { dH: -8.0,  dS: -19.9 },
  'GT': { dH: -8.4,  dS: -22.4 },
  'TA': { dH: -7.2,  dS: -20.4 },
  'TC': { dH: -7.2,  dS: -21.3 },
  'TG': { dH: -8.5,  dS: -22.7 },
  'TT': { dH: -7.9,  dS: -22.2 },
};

const R    = 1.987;   // gas constant cal/mol·K
const T_K  = 310.15;  // 37°C in Kelvin (Cas9 reaction temp)
const T_PLANT = 295.15; // 22°C — plant transformation temp

/**
 * Calculate nearest-neighbor free energy for a sequence
 * @param {string} seq - DNA/RNA sequence
 * @param {Object} params - NN parameter table
 * @param {number} temp - temperature in Kelvin
 * @returns {{ dG, dH, dS }} thermodynamic values
 */
function calcNearestNeighbor(seq, params, temp = T_K) {
  let totalDH = 0;
  let totalDS = 0;

  // Initiation parameters
  const initDH = 0.1;
  const initDS = -2.8;
  totalDH += initDH;
  totalDS += initDS;

  // Sum nearest-neighbor contributions
  for (let i = 0; i < seq.length - 1; i++) {
    const pair = seq[i] + seq[i+1];
    const nn   = params[pair];
    if (nn) {
      totalDH += nn.dH;
      totalDS += nn.dS / 1000; // convert cal to kcal
    }
  }

  // Terminal AT penalty
  if (seq[0] === 'A' || seq[0] === 'T') totalDS -= 0.0014;
  if (seq[seq.length-1] === 'A' || seq[seq.length-1] === 'T') totalDS -= 0.0014;

  const dG = totalDH - temp * totalDS;

  return { dG: parseFloat(dG.toFixed(2)), dH: totalDH, dS: totalDS };
}

/**
 * Calculate gRNA:DNA duplex free energy
 * This is the PRIMARY thermodynamic predictor of CRISPR efficiency
 *
 * More negative ΔG = more stable duplex = better cutting
 * But TOO stable can reduce strand release after cutting
 * Optimal range: -18 to -28 kcal/mol (from published data)
 *
 * @param {string} grnaSeq - 20bp gRNA sequence (DNA notation)
 * @returns {Object} thermodynamic profile
 */
export function calcGRNADuplexThermo(grnaSeq) {
  const seq = grnaSeq.toUpperCase().replace(/U/g, 'T');

  // Full duplex ΔG
  const full   = calcNearestNeighbor(seq, RNA_DNA_NN, T_K);

  // Seed region (positions 1-12, PAM-proximal) — most critical
  const seed   = seq.slice(0, 12);
  const seedNN = calcNearestNeighbor(seed, RNA_DNA_NN, T_K);

  // PAM-proximal region (positions 1-8) — Cas9 initial recognition
  const pamProx    = seq.slice(0, 8);
  const pamProxNN  = calcNearestNeighbor(pamProx, RNA_DNA_NN, T_K);

  // 5' region (positions 13-20) — tolerance for mismatches
  const fivePrime  = seq.slice(12, 20);
  const fivePrimeNN = calcNearestNeighbor(fivePrime, RNA_DNA_NN, T_K);

  // Melting temperature (simplified)
  const gc    = (seq.match(/[GC]/g)||[]).length / seq.length;
  const tm    = -full.dH / (-full.dS + R * Math.log(1e-6)) - 273.15;

  // Efficiency prediction from thermodynamics alone
  // Based on regression from Xu et al. 2015 plant CRISPR data
  let thermoScore = 50; // baseline
  if (full.dG >= -28 && full.dG <= -18)  thermoScore += 25; // optimal range
  else if (full.dG < -28)                 thermoScore += 10; // too stable
  else if (full.dG > -18 && full.dG > -12) thermoScore += 5;  // unstable
  else                                     thermoScore -= 10; // very unstable

  // Seed region stability bonus
  if (seedNN.dG <= -10)  thermoScore += 15;
  else if (seedNN.dG <= -8) thermoScore += 8;
  else                    thermoScore -= 5;

  // PAM-proximal should be moderately stable
  if (pamProxNN.dG >= -8 && pamProxNN.dG <= -5) thermoScore += 10;

  // Self-folding risk
  const selfFold = calcSelfFolding(seq);
  thermoScore -= selfFold.penalty;

  return {
    fullDuplex:    full,
    seedRegion:    seedNN,
    pamProximal:   pamProxNN,
    fivePrime:     fivePrimeNN,
    meltingTemp:   Math.round(tm * 10) / 10,
    thermoScore:   Math.round(Math.min(100, Math.max(0, thermoScore))),
    selfFolding:   selfFold,
    interpretation: interpretThermo(full.dG, seedNN.dG),
  };
}

/**
 * Estimate self-folding risk of gRNA
 * gRNAs that fold on themselves cannot bind target DNA
 * Simplified stem-loop detection
 */
function calcSelfFolding(seq) {
  let maxStemLen = 0;
  const minStemLen = 4;
  const rc = reverseComplement(seq);

  // Check for internal complementarity
  for (let i = 0; i < seq.length - minStemLen; i++) {
    for (let j = i + minStemLen; j < seq.length; j++) {
      const window1 = seq.slice(i, i + minStemLen);
      const window2 = reverseComplement(seq.slice(j, j + minStemLen));
      if (window1 === window2) {
        maxStemLen = Math.max(maxStemLen, minStemLen);
      }
    }
  }

  const penalty = maxStemLen >= 6 ? 15
                : maxStemLen >= 5 ? 8
                : maxStemLen >= 4 ? 3 : 0;

  return {
    maxStemLen,
    penalty,
    risk: maxStemLen >= 6 ? 'High' : maxStemLen >= 4 ? 'Medium' : 'Low',
    color: maxStemLen >= 6 ? '#FF4D6D' : maxStemLen >= 4 ? '#F0C040' : '#00D4AA',
  };
}

function interpretThermo(dG, seedDG) {
  if (dG < -30)
    return { text: 'Very stable duplex — may reduce strand release', color: '#F0C040' };
  if (dG >= -28 && dG <= -18)
    return { text: 'Optimal duplex stability for Cas9 activity', color: '#00D4AA' };
  if (dG > -18 && dG >= -12)
    return { text: 'Moderate stability — reduced efficiency expected', color: '#F0C040' };
  return { text: 'Unstable duplex — low efficiency predicted', color: '#FF4D6D' };
}

/**
 * Calculate melting temperature using nearest-neighbor
 * More accurate than Wallace rule for primers
 * @param {string} seq - primer or short DNA sequence
 * @param {number} primerConc - primer concentration (nM)
 * @returns {number} Tm in °C
 */
export function calcAccurateTm(seq, primerConc = 250) {
  const clean = seq.toUpperCase().replace(/[^ATGC]/g, '');
  if (clean.length < 8) {
    // Wallace rule for very short sequences
    return 2 * (clean.match(/[AT]/g)||[]).length +
           4 * (clean.match(/[GC]/g)||[]).length;
  }

  const { dH, dS } = calcNearestNeighbor(clean, DNA_DNA_NN, T_K);
  const concM = primerConc * 1e-9;
  const tm    = (dH * 1000) / (dS + R * Math.log(concM / 4)) - 273.15;

  return Math.round(tm * 10) / 10;
}

/**
 * Calculate folding energy for a short DNA/RNA window
 * Used to assess secondary structure formation risk
 */
export function calcFoldingEnergy(seq) {
  const clean = seq.toUpperCase().replace(/U/g, 'T');
  const { dG } = calcNearestNeighbor(clean, DNA_DNA_NN, T_PLANT);
  return dG;
}

// Reverse complement helper
function reverseComplement(seq) {
  const comp = { A:'T', T:'A', G:'C', C:'G', N:'N' };
  return seq.toUpperCase().split('').reverse()
    .map(b => comp[b] || 'N').join('');
}

export { calcNearestNeighbor, RNA_DNA_NN, DNA_DNA_NN };
