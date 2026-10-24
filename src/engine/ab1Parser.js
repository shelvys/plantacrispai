// ============================================================
// FILE: src/engine/ab1Parser.js
// PURPOSE: Parse real Sanger sequencing .ab1 binary files
//
// .ab1 is Applied Biosystems binary format used by ALL
// major sequencing facilities worldwide. Every Sanger
// result from your sequencing facility is in this format.
//
// KEY DATA BLOCKS EXTRACTED:
//   PBAS.2   — base calls (the DNA sequence)
//   PCON.2   — quality/confidence scores per base
//   DATA.9-12 — raw fluorescence traces (4 channels)
//   PLOC.2   — peak locations in trace
//   FWO_.1   — filter wheel order (channel → base mapping)
//   SMPL.1   — sample name
//   RUND.1   — run date
// ============================================================

const ABIF_MAGIC = [0x41, 0x42, 0x49, 0x46]; // "ABIF"

// ── Main parser ───────────────────────────────────────────

export function parseAB1(buffer, filename = 'unknown.ab1') {
  try {
    const view  = new DataView(buffer);
    const bytes = new Uint8Array(buffer);

    // Verify magic number "ABIF"
    for (let i = 0; i < 4; i++) {
      if (bytes[i] !== ABIF_MAGIC[i]) {
        return {
          success: false,
          error: `Not a valid .ab1 file — magic bytes mismatch. Got: "${
            String.fromCharCode(...bytes.slice(0,4))
          }", expected: "ABIF"`
        };
      }
    }

    const version   = view.getInt16(4, false);
    const rootEntry = readDirEntry(view, 6);

    // Read all directory entries into a lookup map
    const entries = {};
    for (let i = 0; i < rootEntry.numElements; i++) {
      const offset = rootEntry.dataOffset + i * 28;
      if (offset + 28 > buffer.byteLength) break;
      const entry = readDirEntry(view, offset);
      entries[entry.name + '.' + entry.number] = entry;
    }

    // ── Extract data blocks ───────────────────────────────
    const result = { success: true, filename, version, format: 'AB1' };

    // Sample name
    result.sampleName = readString(view, bytes, entries['SMPL.1'])
                     || filename.replace(/\.ab1$/i, '');

    // DNA sequence (base calls)
    const seq = readSequence(view, bytes, entries['PBAS.2']);
    if (!seq) return { success: false, error: 'Could not read base calls from file' };
    result.sequence = seq;

    // Quality scores (Phred)
    result.quality = readQuality(view, bytes, entries['PCON.2']);

    // Peak locations
    result.peakLocations = readShortArray(view, entries['PLOC.2']);

    // Filter wheel order → maps channel index to base letter
    const fwo = readString(view, bytes, entries['FWO_.1']) || 'GATC';
    result.filterWheelOrder = fwo;

    // Raw fluorescence traces — channels DATA.9 through DATA.12
    const rawTraces = {};
    [9, 10, 11, 12].forEach((tag, idx) => {
      const entry = entries[`DATA.${tag}`];
      if (entry) rawTraces[idx] = readShortArray(view, entry);
    });

    // Map channel index → base letter
    result.traces = {};
    for (let i = 0; i < fwo.length && i < 4; i++) {
      if (rawTraces[i] && rawTraces[i].length > 0) {
        result.traces[fwo[i]] = rawTraces[i];
      }
    }

    // Fall back to synthetic traces if real ones unavailable
    if (Object.keys(result.traces).length < 2) {
      result.traces = generateSyntheticTraces(result.sequence, result.quality);
      result.syntheticTraces = true;
    }

    result.traces = normalizeTraces(result.traces);

    // Run date
    result.runDate = readDate(view, entries['RUND.1']);

    // Quality stats
    result.seqLength  = result.sequence.length;
    result.avgQuality = result.quality.length > 0
      ? Math.round(result.quality.reduce((a,b) => a+b, 0) / result.quality.length)
      : 0;

    // Quality trimming (Phred ≥ 20 window method)
    const { trimStart, trimEnd } = findTrimPoints(result.quality);
    result.trimStart     = trimStart;
    result.trimEnd       = trimEnd;
    result.trimmedSeq    = result.sequence.slice(trimStart, trimEnd);
    result.trimmedLength = result.trimmedSeq.length;
    result.qualityGrade  = assessQuality(result.avgQuality, result.trimmedLength);

    return result;

  } catch (err) {
    return {
      success: false,
      error: `Parse error: ${err.message}. File may be corrupted or use unsupported AB1 version.`
    };
  }
}

// ── Directory entry reader ────────────────────────────────

function readDirEntry(view, offset) {
  try {
    return {
      name:        String.fromCharCode(
                     view.getUint8(offset),   view.getUint8(offset+1),
                     view.getUint8(offset+2), view.getUint8(offset+3)),
      number:      view.getInt32(offset+4,  false),
      elementType: view.getInt16(offset+8,  false),
      elementSize: view.getInt16(offset+10, false),
      numElements: view.getInt32(offset+12, false),
      dataSize:    view.getInt32(offset+16, false),
      dataOffset:  view.getInt32(offset+20, false),
    };
  } catch {
    return { name:'', number:0, numElements:0, dataOffset:0, dataSize:0 };
  }
}

// ── Data readers ──────────────────────────────────────────

function readSequence(view, bytes, entry) {
  if (!entry || entry.numElements === 0) return null;
  try {
    const offset = entry.dataSize <= 4 ? 20 : entry.dataOffset;
    let seq = '';
    for (let i = 0; i < entry.numElements; i++) {
      const b = bytes[offset + i];
      if (b >= 32 && b < 128) seq += String.fromCharCode(b);
    }
    return seq.replace(/[^ATGCN]/gi, 'N').toUpperCase() || null;
  } catch { return null; }
}

function readQuality(view, bytes, entry) {
  if (!entry || entry.numElements === 0) return [];
  try {
    const offset = entry.dataSize <= 4 ? 20 : entry.dataOffset;
    return Array.from({ length: entry.numElements }, (_, i) => bytes[offset + i]);
  } catch { return []; }
}

function readShortArray(view, entry) {
  if (!entry || entry.numElements === 0) return [];
  try {
    const offset = entry.dataSize <= 4 ? 20 : entry.dataOffset;
    return Array.from({ length: entry.numElements }, (_, i) =>
      view.getInt16(offset + i * 2, false));
  } catch { return []; }
}

function readString(view, bytes, entry) {
  if (!entry || entry.numElements === 0) return '';
  try {
    let offset = entry.dataSize <= 4 ? 20 : entry.dataOffset;
    let str = '';
    const isPascal = entry.elementType === 18 || entry.elementType === 2;
    const start    = isPascal ? bytes[offset++] : 0;
    const len      = isPascal ? start : entry.numElements;
    for (let i = 0; i < len; i++) {
      const ch = bytes[offset + i];
      if (ch === 0) break;
      if (ch >= 32) str += String.fromCharCode(ch);
    }
    return str.trim();
  } catch { return ''; }
}

function readDate(view, entry) {
  if (!entry) return null;
  try {
    const o = entry.dataOffset;
    const y = view.getInt16(o, false);
    const m = view.getUint8(o+2);
    const d = view.getUint8(o+3);
    return `${y}-${String(m).padStart(2,'0')}-${String(d).padStart(2,'0')}`;
  } catch { return null; }
}

// ── Quality processing ────────────────────────────────────

/**
 * Four trimming modes — exported so UI and engine can both call them.
 *
 * MODE 'none'     — raw untrimmed sequence (user trusts read quality)
 * MODE 'standard' — Phred ≥ 20, sliding window 10 bp (lab standard default)
 * MODE 'quality'  — user-set Phred threshold + window (Q25/Q30 for high-qual data)
 * MODE 'manual'   — user drags exact start/end position sliders
 *
 * WHY THIS MATTERS:
 *   The Cas9 cut site is often within the first 100-150 bp of a Sanger read,
 *   which is the lowest-quality region. Aggressive auto-trimming can remove
 *   the entire cut-site region, causing the engine to see only a perfect
 *   downstream trace → misclassified as WT or Homo. This function gives the
 *   researcher full control and always warns when trim points are dangerously
 *   close to the predicted cut site.
 *
 * @param {string}   sequence  - Full raw sequence string
 * @param {number[]} quality   - Phred score array (same length as sequence)
 * @param {string}   mode      - 'none' | 'standard' | 'quality' | 'manual'
 * @param {Object}   options   - { phredThreshold, window, manualStart, manualEnd }
 * @param {number}   [cutSite] - Predicted cut site position for warning generation
 * @returns {Object} { trimStart, trimEnd, trimmedSeq, mode, warning }
 */
export function applyTrim(sequence, quality, mode = 'standard', options = {}, cutSite = null) {
  const n = (sequence || '').length;
  if (n === 0) return { trimStart: 0, trimEnd: 0, trimmedSeq: '', mode, warning: null };

  let trimStart = 0;
  let trimEnd   = n;

  if (mode === 'none') {
    // No trimming — use everything
    trimStart = 0;
    trimEnd   = n;

  } else if (mode === 'standard') {
    // Classic Phred ≥ 20, 10-base sliding window (NCBI / Ewing et al. 1998)
    const PHRED = 20, WIN = 10;
    if (quality && quality.length >= WIN) {
      for (let i = 0; i <= n - WIN; i++) {
        if (quality.slice(i, i + WIN).reduce((a, b) => a + b, 0) / WIN >= PHRED) {
          trimStart = i; break;
        }
      }
      for (let i = n - WIN; i >= trimStart; i--) {
        if (quality.slice(i, i + WIN).reduce((a, b) => a + b, 0) / WIN >= PHRED) {
          trimEnd = i + WIN; break;
        }
      }
    }

  } else if (mode === 'quality') {
    // User-defined Phred threshold + window
    const PHRED = Math.max(1,  options.phredThreshold || 20);
    const WIN   = Math.max(3,  options.window         || 10);
    if (quality && quality.length >= WIN) {
      for (let i = 0; i <= n - WIN; i++) {
        if (quality.slice(i, i + WIN).reduce((a, b) => a + b, 0) / WIN >= PHRED) {
          trimStart = i; break;
        }
      }
      for (let i = n - WIN; i >= trimStart; i--) {
        if (quality.slice(i, i + WIN).reduce((a, b) => a + b, 0) / WIN >= PHRED) {
          trimEnd = i + WIN; break;
        }
      }
    }

  } else if (mode === 'manual') {
    // User-dragged slider positions
    trimStart = Math.max(0,          Math.min(options.manualStart || 0, n - 1));
    trimEnd   = Math.max(trimStart + 1, Math.min(options.manualEnd   || n, n));
  }

  // Clamp to sequence bounds
  trimStart = Math.max(0, trimStart);
  trimEnd   = Math.min(n, trimEnd);
  if (trimEnd <= trimStart) { trimStart = 0; trimEnd = n; }

  const trimmedSeq = (sequence || '').slice(trimStart, trimEnd);

  // ── Cut-site proximity warning ─────────────────────────────────────────
  // If the cut site is within 20 bp of either trim edge, warn the user.
  // This is the most common cause of wrong classification from trimming.
  let warning = null;
  if (cutSite != null) {
    const margin = Math.min(cutSite - trimStart, trimEnd - cutSite);
    if (cutSite < trimStart || cutSite > trimEnd) {
      warning = `⚠ Cut site (~${cutSite} bp) is OUTSIDE the trimmed region (${trimStart}–${trimEnd}). Analysis will be unreliable — use "No trimming" or adjust trim points.`;
    } else if (margin < 20) {
      warning = `⚠ Cut site (~${cutSite} bp) is only ${margin} bp from the trim boundary. Consider widening trim range to include more context around the cut site.`;
    }
  }

  return { trimStart, trimEnd, trimmedSeq, mode, warning };
}

function findTrimPoints(quality) {
  // Legacy wrapper — used internally by parseAB1 for initial parse
  if (!quality || quality.length === 0) return { trimStart: 0, trimEnd: quality?.length || 0 };
  const r = applyTrim('X'.repeat(quality.length), quality, 'standard', {}, null);
  return { trimStart: r.trimStart, trimEnd: r.trimEnd };
}

function assessQuality(avgQ, seqLen) {
  if (avgQ >= 40 && seqLen >= 400) return { grade:'Excellent',  color:'#00d4aa', score:5 };
  if (avgQ >= 30 && seqLen >= 300) return { grade:'Good',       color:'#4d9fff', score:4 };
  if (avgQ >= 20 && seqLen >= 200) return { grade:'Acceptable', color:'#f0c040', score:3 };
  if (avgQ >= 15 && seqLen >= 100) return { grade:'Poor',       color:'#ff8c42', score:2 };
  return                                  { grade:'Failed',     color:'#ff4d6d', score:1 };
}

// ── Trace processing ──────────────────────────────────────

function normalizeTraces(traces) {
  const out = {};
  for (const [base, data] of Object.entries(traces)) {
    if (!data || data.length === 0) continue;
    const max = Math.max(...data.map(v => Math.abs(v)));
    out[base] = max === 0 ? data.map(()=>0) : data.map(v => Math.max(0,v)/max);
  }
  return out;
}

function generateSyntheticTraces(sequence, quality) {
  if (!sequence) return {};
  const LEN   = sequence.length * 12;
  const bases = ['A','T','G','C'];
  const traces = {};
  bases.forEach(b => { traces[b] = Array.from({length:LEN}, ()=>Math.random()*0.03); });

  for (let i = 0; i < sequence.length; i++) {
    const base   = sequence[i];
    const qual   = quality[i] || 30;
    const center = Math.round((i + 0.5) * 12);
    const height = Math.min(0.95, 0.3 + (qual/60)*0.65);
    const W = 5;
    if (!traces[base]) continue;
    for (let j = -W*2; j <= W*2; j++) {
      const pos = center + j;
      if (pos < 0 || pos >= LEN) continue;
      const peak = height * Math.exp(-(j*j)/(2*W*W));
      traces[base][pos] = Math.max(traces[base][pos], peak);
      bases.forEach(b => {
        if (b !== base) {
          traces[b][pos] = Math.max(traces[b][pos], peak*(0.04+Math.random()*0.04));
        }
      });
    }
  }
  return traces;
}

// ── Public API ────────────────────────────────────────────

export function readAB1File(file) {
  return new Promise((resolve, reject) => {
    if (!file) return reject(new Error('No file provided'));
    if (!file.name.toLowerCase().endsWith('.ab1'))
      return reject(new Error(`Expected .ab1 file, got: ${file.name}`));
    const reader = new FileReader();
    reader.onload  = e => {
      const r = parseAB1(e.target.result, file.name);
      r.success ? resolve(r) : reject(new Error(r.error));
    };
    reader.onerror = () => reject(new Error('Failed to read file'));
    reader.readAsArrayBuffer(file);
  });
}

export function readMultipleAB1Files(files) {
  return Promise.all(
    Array.from(files).map(f =>
      readAB1File(f).catch(err => ({ success:false, filename:f.name, error:err.message }))
    )
  );
}

// Convert parsed AB1 → PlantaCRISPAI sample format
// trimMode / trimOptions allow the UI to re-trim without re-parsing the binary
export function ab1ToSample(parsed, generation = 'T0', ploidy = 'diploid',
                             trimMode = 'standard', trimOptions = {}) {
  if (!parsed.success) return null;

  // Apply the requested trim mode (may differ from the default standard trim)
  const trimResult = applyTrim(
    parsed.sequence,
    parsed.quality,
    trimMode,
    trimOptions,
    null  // cut site unknown at parse time; warning added after analysis
  );

  return {
    id:          parsed.sampleName || parsed.filename.replace(/\.ab1$/i, ''),
    // sequence = trimmed (what the engine uses for analysis)
    sequence:    trimResult.trimmedSeq,
    // fullSeq  = raw complete read (always preserved for re-trimming)
    fullSeq:     parsed.sequence,
    quality:     parsed.quality,
    traces:      parsed.traces,
    peakLocs:    parsed.peakLocations,
    // Trim metadata
    trimStart:   trimResult.trimStart,
    trimEnd:     trimResult.trimEnd,
    trimMode:    trimMode,
    trimOptions: trimOptions,
    trimWarning: trimResult.warning,
    // Quality stats
    avgQuality:  parsed.avgQuality,
    qualGrade:   parsed.qualityGrade,
    seqLength:   trimResult.trimmedSeq.length,
    rawSeqLength: parsed.seqLength,
    // File metadata
    generation,
    ploidy,
    source:      'ab1',
    runDate:     parsed.runDate,
    filename:    parsed.filename,
  };
}

/**
 * Re-apply trimming to an existing sample without re-parsing the .ab1 file.
 * Called when the user changes trim mode or slider positions in the UI.
 *
 * @param {Object} sample      - Existing sample with fullSeq + quality
 * @param {string} trimMode    - New mode: 'none'|'standard'|'quality'|'manual'
 * @param {Object} trimOptions - { phredThreshold, window, manualStart, manualEnd }
 * @param {number} [cutSite]   - Predicted cut site for proximity warning
 * @returns {Object} Updated sample with new trim applied
 */
export function retrimSample(sample, trimMode, trimOptions = {}, cutSite = null) {
  const rawSeq = sample.fullSeq || sample.sequence || '';
  const trimResult = applyTrim(rawSeq, sample.quality || [], trimMode, trimOptions, cutSite);
  return {
    ...sample,
    sequence:    trimResult.trimmedSeq,
    trimStart:   trimResult.trimStart,
    trimEnd:     trimResult.trimEnd,
    trimMode,
    trimOptions,
    trimWarning: trimResult.warning,
    seqLength:   trimResult.trimmedSeq.length,
  };
}
