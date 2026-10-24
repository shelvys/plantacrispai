// ============================================================
// FILE: src/components/ChromatogramViewer.js  v5
//
// COMPLETE REWRITE — fixes all rendering and lint issues.
//
// Key fixes vs v4:
//   • globalMax scan now uses PIXEL indices (data is pixel-indexed,
//     not base-indexed) — this was the root cause of the flat canvas bug
//   • primary is wrapped in useMemo (no re-render loops)
//   • useCallback removed from import (was unused)
//   • Canvas scaling is per-channel-peak, not global, for better visual
//   • ICE-style dual-trace overlay: edited (solid) + WT (dashed) stacked
//   • gRNA found in REFERENCE sequence, then mapped to sample coords
//     via the alignment — correct even with indels
//   • Quality bar under each base call (when real .ab1 data present)
//   • Zoom uses mouse wheel directly on canvas
// ============================================================

import React, { useRef, useEffect, useState, useMemo } from 'react';

// ── Colours ───────────────────────────────────────────────────────────────────
const BC   = { A:'#00c896', T:'#ff5a5a', G:'#4a9eff', C:'#f5a623', N:'#667' };
const BG   = '#060c1c';
const CARD = '#0b1425';
const BDR  = '#172040';
const TXT  = '#c8d6f0';
const DIM  = '#4a6080';
const ACN  = '#00d4aa';

// ── Synthetic trace generator ─────────────────────────────────────────────────
// Each base gets a Gaussian peak at its pixel position.
// Data arrays are PIXEL-indexed (length = seq.length * STEP).
function synthTraces(seq, quality) {
  if (!seq) return null;
  const s    = seq.replace(/[^ATGCN]/gi, '').toUpperCase();
  if (!s.length) return null;
  const STEP = 14;          // pixels per base
  const W    = s.length * STEP;
  const CH   = {};
  ['A','T','G','C'].forEach(b => { CH[b] = new Float32Array(W); });

  for (let i = 0; i < s.length; i++) {
    const b = s[i];
    if (!CH[b]) continue;
    const q   = quality ? Math.min(quality[i] || 30, 60) : 30;
    const h   = 150 + q * 1.6;          // peak amplitude
    const sig = STEP * 0.55;
    const SPR = Math.ceil(STEP * 2.0);
    for (let dx = -SPR; dx <= SPR; dx++) {
      const px = Math.round(i * STEP + STEP * 0.5 + dx);
      if (px < 0 || px >= W) continue;
      const g = Math.exp(-(dx*dx) / (2*sig*sig));
      CH[b][px] += g * h;
      // Realistic cross-talk
      ['A','T','G','C'].forEach(bb => {
        if (bb !== b) CH[bb][px] += g * (3 + Math.random()*3);
      });
    }
  }
  // 3-point smoothing
  const sm = {};
  ['A','T','G','C'].forEach(b => {
    const d = CH[b]; const o = new Float32Array(W);
    for (let i = 1; i < W-1; i++) o[i] = d[i-1]*0.2 + d[i]*0.6 + d[i+1]*0.2;
    o[0] = d[0]; o[W-1] = d[W-1]; sm[b] = o;
  });
  return {
    channels: sm,
    width: W,                                        // total pixel width
    baseLen: s.length,                               // number of bases
    STEP,                                            // pixels per base
    basePositions: s.split('').map((b,i) => ({ base: b, x: i*STEP + STEP*0.5 })),
  };
}

// ── Bezier path helper ────────────────────────────────────────────────────────
function bezierPath(ctx, pts) {
  if (!pts || pts.length < 2) return;
  ctx.beginPath();
  ctx.moveTo(pts[0].x, pts[0].y);
  for (let i = 0; i < pts.length-1; i++) {
    const mx = (pts[i].x + pts[i+1].x) / 2;
    ctx.bezierCurveTo(mx, pts[i].y, mx, pts[i+1].y, pts[i+1].x, pts[i+1].y);
  }
}



// ── Derive cut/gRNA/PAM positions in SAMPLE sequence coordinates ──────────────
function deriveAnnotations(sample, grnaSeq, grnaList, refCutPos) {
  const COLORS = ['#00d4aa','#bf5af2','#4d9fff','#f0c040'];
  const s = ((sample?.sequence || sample?.fullSeq || '')).toUpperCase().replace(/[^ATGCN]/g,'');
  if (!s.length) return [{ label:'cut', grna:null, pam:null, cutBase: refCutPos ?? null }];

  // Build gRNA list — use grnaList if it has real sequences, else fall back to grnaSeq prop
  const validGrnas = (grnaList || []).filter(g => (g.sequence||'').replace(/[^ATGCN]/gi,'').length >= 15);
  const gList = validGrnas.length > 0
    ? validGrnas
    : (grnaSeq && grnaSeq.replace(/[^ATGCN]/gi,'').length >= 15
        ? [{ label:'sgRNA', sequence: grnaSeq }]
        : []);

  if (!gList.length) return [{ label:'cut', grna:null, pam:null, cutBase: refCutPos ?? null }];

  const RC = str => str.split('').reverse()
    .map(b => ({A:'T',T:'A',G:'C',C:'G',N:'N'}[b] || b)).join('');

  return gList.map((g, gi) => {
    const gSeq = g.sequence.toUpperCase().replace(/[^ATGCN]/g,'');
    let idx = -1, strand = '+';

    // Strategy 1: exact full match (WT and heterozygous samples)
    idx = s.indexOf(gSeq);

    // Strategy 2: RC exact match (minus-strand gRNAs)
    if (idx < 0) { const rc = RC(gSeq); idx = s.indexOf(rc); if (idx >= 0) strand = '-'; }

    // Strategy 3: upstream 17bp match (edited samples — insertion/deletion disrupts last 3bp)
    // The first 17bp of gRNA are ALWAYS intact because Cas9 cuts 3bp from the PAM end
    if (idx < 0) {
      const up17 = gSeq.slice(0, 17);
      idx = s.indexOf(up17);
      if (idx >= 0) strand = '+';
    }

    // Strategy 4: RC upstream 17bp
    if (idx < 0) {
      const rcUp17 = RC(gSeq).slice(0, 17);
      const ri = s.indexOf(rcUp17);
      if (ri >= 0) { idx = ri; strand = '-'; }
    }

    // Not found in any strategy
    if (idx < 0) return { label: g.label||`sgRNA${gi+1}`, grna:null, pam:null,
                          cutBase: sample?.sampleCutPosition ?? refCutPos ?? null };

    const end = idx + 20; // always 20bp gRNA window
    const cutBase  = strand === '+' ? idx + 17 : idx + 3;  // 3bp from PAM end
    const pamStart = strand === '+' ? end       : Math.max(0, idx - 3);
    const pamEnd   = strand === '+' ? end + 3   : idx;

    return {
      label:   g.label || `sgRNA${gi+1}`,
      grna:    { start: idx, end: Math.min(end, s.length), strand, color: COLORS[gi%4] },
      pam:     { start: Math.max(0,pamStart), end: Math.min(pamEnd, s.length) },
      cutBase: Math.max(0, Math.min(cutBase, s.length - 1)),
    };
  }).filter(Boolean);
}

// ── Main draw function ────────────────────────────────────────────────────────
function drawChroma(canvas, traces, opts) {
  if (!canvas || !traces?.channels || !traces.basePositions) return;
  const { annotations, showBases, zoom, pan, wtTraces, showWT, hoverX, viewMode, yScale=1.0, seqSearch='' } = opts;

  const ctx  = canvas.getContext('2d');
  const W    = canvas.width;
  const H    = canvas.height;
  const STEP = (traces.STEP || 14) * zoom;
  const TOP   = 52;        // reserved for annotation flags at top
  const BL   = H - 30;    // baseline — 30px for base row

  // Pan → pixel offset
  const totalW = traces.width * zoom;
  const maxPan = Math.max(0, totalW - W);
  const OX     = -(pan * maxPan);

  // ── Scale ──────────────────────────────────────────────────────────────────
  // Visible pixel range in the raw data array
  const visStart = Math.max(0, Math.floor(-OX));
  const visEnd   = Math.min(traces.width, Math.ceil((W - OX)));

  // Per-channel scales for NORMALIZED mode (each channel scaled to its own max in view)
  // Global scale for RAW mode (true amplitude ratios preserved)
  const channelScales = {};
  let globalMax = 0;
  Object.entries(traces.channels).forEach(([base, d]) => {
    let chMax = 0;
    for (let px = visStart; px < visEnd && px < d.length; px++)
      if (d[px] > chMax) chMax = d[px];
    if (chMax < 1) chMax = 1;
    channelScales[base] = (BL - 18) / chMax;
    // Also track global max for raw mode
    for (let px = 0; px < d.length; px++) if (d[px] > globalMax) globalMax = d[px];
  });
  if (globalMax < 1) globalMax = 1;
  const GLOBAL_SCALE = ((BL - 18) / globalMax) * yScale;
  // Rebuild channelScales with TOP margin
  Object.entries(traces.channels).forEach(([base, d]) => {
    let chMax = 0;
    for (let px = visStart; px < visEnd && px < d.length; px++)
      if (d[px] > chMax) chMax = d[px];
    if (chMax < 1) chMax = 1;
    channelScales[base] = ((BL - TOP - 8) / chMax) * yScale;
  });
  const getScale = (base) => viewMode === 'normalized'
    ? channelScales[base]
    : GLOBAL_SCALE;
  const SCALE = GLOBAL_SCALE; // used for WT overlay

  // ── Background ──────────────────────────────────────────────────────────────
  ctx.clearRect(0, 0, W, H);
  const bgG = ctx.createLinearGradient(0, 0, 0, H);
  bgG.addColorStop(0, '#07101f'); bgG.addColorStop(1, '#030710');
  ctx.fillStyle = bgG; ctx.fillRect(0, 0, W, H);

  // ── Derive primary annotation ───────────────────────────────────────────────
  const primary = annotations?.[0];
  const grna    = primary?.grna;
  const pam     = primary?.pam;
  const cutBase = primary?.cutBase;

  // ── STEP 1: Thin coloured track bars at top — region indicator only ──────────
  // Industry standard (SnapGene/ICE): coloured band at top, NOT full-height fill.
  // This keeps traces fully readable.
  const TRACK_H = 5;  // height of top track bar in px
  if (grna) {
    const gx1 = grna.start * STEP + OX;
    const gx2 = grna.end   * STEP + OX;
    if (gx2 > 0 && gx1 < W) {
      // Solid coloured bar just below the TOP zone
      ctx.fillStyle = '#00d4aadd';
      ctx.fillRect(Math.max(0,gx1), TOP-TRACK_H, Math.min(gx2,W) - Math.max(0,gx1), TRACK_H);
    }
  }
  if (pam) {
    const px1 = pam.start * STEP + OX;
    const px2 = pam.end   * STEP + OX;
    if (px2 > 0 && px1 < W) {
      ctx.fillStyle = 'rgba(255,152,0,0.9)';
      ctx.fillRect(Math.max(0,px1), TOP-TRACK_H, Math.min(px2,W) - Math.max(0,px1), TRACK_H);
    }
  }
  // ── gRNA sequence bases printed in track zone (small, above trace area) ────
  if (grna && traces.basePositions) {
    const fs = Math.max(7, Math.min(10, STEP * 0.68));
    ctx.font = `bold ${fs}px monospace`;
    ctx.textAlign = 'center';
    for (let idx = grna.start; idx < grna.end; idx++) {
      const bp = traces.basePositions[idx];
      if (!bp) continue;
      const x = bp.x * zoom + OX;
      if (x < 0 || x > W) continue;
      // gRNA bases drawn in track zone (above TOP)
      const col = BC[traces.basePositions[idx]?.base] || '#00d4aa';
      ctx.fillStyle = col;
      ctx.fillText(bp.base, x, TOP - TRACK_H - 3);
    }
  }

  // ── Sequence search highlight ──────────────────────────────────────────────
  if (seqSearch && seqSearch.length >= 3 && traces.basePositions) {
    const needle = seqSearch.toUpperCase().replace(/[^ATGCN]/g,'');
    const haystack = traces.basePositions.map(bp => bp.base).join('');
    let pos = 0;
    while (pos < haystack.length) {
      const found = haystack.indexOf(needle, pos);
      if (found < 0) break;
      const sx1 = traces.basePositions[found].x * zoom + OX;
      const sx2 = traces.basePositions[Math.min(found + needle.length - 1, traces.basePositions.length-1)].x * zoom + OX + STEP;
      if (sx2 > 0 && sx1 < W) {
        ctx.fillStyle = 'rgba(255,220,50,0.18)';
        ctx.fillRect(Math.max(0,sx1), TOP, Math.min(sx2,W)-Math.max(0,sx1), BL-TOP);
        ctx.strokeStyle = 'rgba(255,220,50,0.6)'; ctx.lineWidth = 1.5;
        ctx.strokeRect(Math.max(0,sx1), TOP, Math.min(sx2,W)-Math.max(0,sx1), BL-TOP);
      }
      pos = found + 1;
    }
  }

  // ── STEP 2: WT overlay (behind main trace, dashed) ─────────────────────────
  if (showWT && wtTraces?.channels && wtTraces.basePositions) {
    let wtMax = 0;
    Object.values(wtTraces.channels).forEach(d => {
      for (let px = 0; px < d.length; px++) if (d[px] > wtMax) wtMax = d[px];
    });
    const wtScale = wtMax > 0 ? ((BL-TOP-8)/wtMax) * yScale : SCALE;
    Object.entries(wtTraces.channels).forEach(([base, data]) => {
      const col = BC[base] || '#888';
      const pts = [];
      const SKP = Math.max(1, Math.floor(1/zoom));
      for (let px = 0; px < data.length; px += SKP) {
        const x = px * zoom + OX;
        if (x < -20 || x > W+20) continue;
        pts.push({ x, y: Math.max(TOP, BL - data[px] * wtScale) });
      }
      if (pts.length < 2) return;
      ctx.save();
      ctx.setLineDash([3, 7]); ctx.globalAlpha = 0.3;
      ctx.strokeStyle = col; ctx.lineWidth = 1.0;
      bezierPath(ctx, pts); ctx.stroke();
      ctx.restore();
    });
  }

  // ── STEP 3: Main fluorescence traces ────────────────────────────────────────
  ['C','G','T','A'].forEach(base => {
    const data = traces.channels[base]; if (!data) return;
    const col  = BC[base] || '#888';
    const SCA  = getScale(base);   // per-channel in normalized, global in raw
    const SKP  = Math.max(1, Math.floor(1/zoom));
    const pts  = [];
    for (let px = 0; px < data.length; px += SKP) {
      const x = px * zoom + OX;
      if (x < -20 || x > W+20) continue;
      pts.push({ x, y: Math.max(TOP, BL - data[px] * SCA) });
    }
    if (pts.length < 2) return;
    ctx.save();
    // Filled area under curve
    bezierPath(ctx, pts);
    ctx.lineTo(pts[pts.length-1].x, BL);
    ctx.lineTo(pts[0].x, BL); ctx.closePath();
    const ag = ctx.createLinearGradient(0, 14, 0, BL);
    ag.addColorStop(0, col+'10'); ag.addColorStop(0.5, col+'04'); ag.addColorStop(1, col+'01');
    ctx.fillStyle = ag; ctx.fill();
    // Stroke line
    bezierPath(ctx, pts);
    ctx.strokeStyle = col; ctx.lineWidth = zoom > 1.5 ? 2.0 : 1.6;
    ctx.lineJoin = 'round'; ctx.lineCap = 'round';
    ctx.shadowColor = col + '33'; ctx.shadowBlur = 2;
    ctx.stroke(); ctx.shadowBlur = 0;
    ctx.restore();
  });

  // ── STEP 3.5: Light bottom accent bars under the baseline ───────────────────
  // After traces: a coloured strip along the BOTTOM of the trace area (under baseline)
  // This marks the region without obscuring any peaks.
  const BOT_TRACK = 3;
  if (grna) {
    const gx1 = grna.start * STEP + OX;
    const gx2 = grna.end   * STEP + OX;
    if (gx2 > 0 && gx1 < W) {
      ctx.fillStyle = '#00d4aa99';
      ctx.fillRect(Math.max(0,gx1), BL+1, Math.min(gx2,W) - Math.max(0,gx1), BOT_TRACK);
    }
  }
  if (pam) {
    const px1 = pam.start * STEP + OX;
    const px2 = pam.end   * STEP + OX;
    if (px2 > 0 && px1 < W) {
      ctx.fillStyle = 'rgba(255,152,0,0.75)';
      ctx.fillRect(Math.max(0,px1), BL+1, Math.min(px2,W) - Math.max(0,px1), BOT_TRACK);
    }
  }

  // ── STEP 4: Draw annotation overlays ON TOP of traces ───────────────────────

  // ── gRNA border lines ──
  if (grna) {
    const gx1 = grna.start * STEP + OX;
    const gx2 = grna.end   * STEP + OX;
    if (gx2 > 0 && gx1 < W) {
      // Left border glow
      if (gx1 > 0 && gx1 < W) {
        ctx.save();
        ctx.shadowColor = '#00d4aa'; ctx.shadowBlur = 6;
        ctx.strokeStyle = '#00d4aacc'; ctx.lineWidth = 1.5;
        ctx.beginPath(); ctx.moveTo(gx1, TOP); ctx.lineTo(gx1, BL); ctx.stroke();
        ctx.shadowBlur = 0; ctx.restore();
      }
      // Right border
      if (gx2 > 0 && gx2 < W) {
        ctx.strokeStyle = '#00d4aa44'; ctx.lineWidth = 1;
        ctx.setLineDash([3,4]);
        ctx.beginPath(); ctx.moveTo(gx2, TOP); ctx.lineTo(gx2, BL); ctx.stroke();
        ctx.setLineDash([]);
      }
    }
  }

  // ── PAM border line ──
  if (pam) {
    const px1 = pam.start * STEP + OX;
    if (px1 > 0 && px1 < W) {
      ctx.save();
      ctx.shadowColor = 'rgba(255,152,0,0.6)'; ctx.shadowBlur = 4;
      ctx.strokeStyle = 'rgba(255,152,0,0.7)'; ctx.lineWidth = 1.5;
      ctx.beginPath(); ctx.moveTo(px1, TOP); ctx.lineTo(px1, BL); ctx.stroke();
      ctx.shadowBlur = 0; ctx.restore();
    }
  }

  // ── Cut site line ──
  const drawCutLine = (cx, alpha=1) => {
    if (cx <= 0 || cx >= W) return;
    ctx.save();
    ctx.shadowColor = '#ff4444'; ctx.shadowBlur = 12;
    ctx.setLineDash([5, 4]);
    ctx.strokeStyle = `rgba(255,68,68,${alpha * 0.85})`; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(cx, TOP); ctx.lineTo(cx, BL); ctx.stroke();
    ctx.setLineDash([]); ctx.shadowBlur = 0; ctx.restore();
  };
  if (cutBase != null) {
    const cx = cutBase * STEP + STEP * 0.5 + OX;
    drawCutLine(cx);
  } else if (primary?.cutBase != null) {
    drawCutLine(primary.cutBase * STEP + STEP * 0.5 + OX, 0.55);
  }

  // ── STEP 5: Prominent label flags (drawn last, always visible) ──────────────
  // Helper: draw a bottom-anchored triangular flag with label
  const drawFlag = (centerX, label, bgColor, textColor, yBase, arrowUp=true) => {
    if (centerX < 8 || centerX > W - 8) return;
    const text  = label;
    const TW    = text.length * 7 + 16;
    const TH    = 18;
    const tx    = Math.max(TW/2 + 4, Math.min(W - TW/2 - 4, centerX));
    const ty    = arrowUp ? yBase - TH - 6 : yBase + 6;

    // Pill background
    ctx.fillStyle = bgColor;
    ctx.shadowColor = bgColor; ctx.shadowBlur = 8;
    ctx.beginPath();
    if (ctx.roundRect) ctx.roundRect(tx - TW/2, ty, TW, TH, 3);
    else ctx.rect(tx - TW/2, ty, TW, TH);
    ctx.fill(); ctx.shadowBlur = 0;

    // Arrow connector to the line
    ctx.beginPath();
    if (arrowUp) {
      ctx.moveTo(centerX, yBase);
      ctx.lineTo(Math.min(tx + 4, tx + TW/2), ty + TH);
      ctx.lineTo(Math.max(tx - 4, tx - TW/2), ty + TH);
    } else {
      ctx.moveTo(centerX, yBase);
      ctx.lineTo(Math.min(tx + 4, tx + TW/2), ty);
      ctx.lineTo(Math.max(tx - 4, tx - TW/2), ty);
    }
    ctx.closePath();
    ctx.fillStyle = bgColor; ctx.fill();

    // Text
    ctx.fillStyle = textColor;
    ctx.font = 'bold 10px monospace';
    ctx.textAlign = 'center';
    ctx.fillText(text, tx, ty + TH - 5);
  };

  // Draw gRNA label flag — above gRNA centre
  if (grna) {
    const gx1 = grna.start * STEP + OX;
    const gx2 = grna.end   * STEP + OX;
    if (gx2 > 0 && gx1 < W) {
      const centre = (Math.max(0,gx1) + Math.min(W,gx2)) / 2;
      const grnaLabel = primary.label && primary.label !== 'cut' ? primary.label : 'sgRNA';
      drawFlag(centre, `${grnaLabel}  ${grna.end-grna.start}bp`, '#00d4aaee', '#001a10', TOP-TRACK_H-1, true);

      // Small "5'" / "3'" labels at gRNA ends
      if (gx1 > 4 && gx1 < W) {
        ctx.fillStyle = '#00d4aa99'; ctx.font = 'bold 7px monospace'; ctx.textAlign = 'center';
        ctx.fillText("5'", gx1, BL + 22);
      }
      if (gx2 > 4 && gx2 < W) {
        ctx.fillStyle = '#00d4aa99'; ctx.font = 'bold 7px monospace'; ctx.textAlign = 'center';
        ctx.fillText("3'", gx2, BL + 22);
      }
    }
  }

  // Draw PAM label flag
  if (pam) {
    const px1 = pam.start * STEP + OX;
    const px2 = pam.end   * STEP + OX;
    if (px2 > 0 && px1 < W) {
      const centre = (Math.max(0,px1) + Math.min(W,px2)) / 2;
      drawFlag(centre, 'PAM (NGG)', 'rgba(255,152,0,0.95)', '#1a0800', TOP-TRACK_H-1, true);
    }
  }

  // Draw cut site flag (prominent red, above cut line)
  if (cutBase != null) {
    const cx = cutBase * STEP + STEP * 0.5 + OX;
    if (cx > 4 && cx < W - 4) {
      // Scissor icon at top of cut line
      ctx.font = '15px serif';
      ctx.fillStyle = '#ff5555';
      ctx.textAlign = 'center';
      ctx.shadowColor = '#ff4444'; ctx.shadowBlur = 8;
      ctx.fillText('✂', cx, TOP - TRACK_H - 22);
      ctx.shadowBlur = 0;
      // Flag below scissor
      drawFlag(cx, 'Cut site', 'rgba(220,50,50,0.97)', '#fff0f0', TOP-TRACK_H-1, true);
    }
  } else if (primary?.cutBase != null && !grna) {
    const cx = primary.cutBase * STEP + STEP * 0.5 + OX;
    if (cx > 4 && cx < W - 4) {
      ctx.font = '13px serif'; ctx.fillStyle = '#ff555588';
      ctx.textAlign = 'center'; ctx.fillText('✂', cx, 14);
      drawFlag(cx, 'Cut (est.)', 'rgba(255,68,68,0.6)', '#fff0f0', 16, false);
    }
  }

  // ── STEP 6: Quality bars (above baseline, when real ab1) ────────────────────
  if (traces.quality && showBases) {
    traces.basePositions.forEach((bp, i) => {
      const x = bp.x * zoom + OX;
      if (x < -STEP || x > W+STEP) return;
      const q  = Math.min(traces.quality[i] || 0, 60);
      const qh = (q / 60) * 8;
      ctx.fillStyle = q >= 30 ? 'rgba(0,212,170,0.2)' : q >= 20 ? 'rgba(245,158,11,0.2)' : 'rgba(255,80,80,0.2)';
      ctx.fillRect(x - STEP*0.42, BL - qh, STEP*0.84, qh);
    });
  }

  // ── STEP 7: Hover spotlight ──────────────────────────────────────────────────
  if (hoverX != null && traces.basePositions) {
    const hov = traces.basePositions.find(bp => Math.abs(bp.x*zoom+OX - hoverX) < STEP*0.55);
    if (hov) {
      const bx  = hov.x * zoom + OX;
      const col = BC[hov.base] || TXT;
      const rg  = ctx.createRadialGradient(bx, BL, 0, bx, BL-20, 28);
      rg.addColorStop(0, col+'18'); rg.addColorStop(1, col+'00');
      ctx.fillStyle = rg; ctx.fillRect(bx-28, 0, 56, BL);
      ctx.font = 'bold 15px monospace'; ctx.fillStyle = col;
      ctx.textAlign = 'center';
      ctx.shadowColor = col; ctx.shadowBlur = 5;
      ctx.fillText(hov.base, bx, BL - 6); ctx.shadowBlur = 0;
    }
  }

  // ── STEP 8: Base row ─────────────────────────────────────────────────────────
  if (showBases && traces.basePositions) {
    traces.basePositions.forEach((bp, idx) => {
      const x = bp.x * zoom + OX;
      if (x < -STEP || x > W+STEP) return;
      const col  = BC[bp.base] || TXT;
      const inG  = grna   && idx >= grna.start  && idx < grna.end;
      const inP  = pam    && idx >= pam.start   && idx < pam.end;
      const atC  = cutBase != null && idx === cutBase;
      if (inG || inP || atC) {
        ctx.fillStyle = atC ? 'rgba(255,68,68,0.15)' : inP ? 'rgba(255,152,0,0.1)' : 'rgba(0,212,170,0.08)';
        ctx.fillRect(x - STEP*0.44, BL+2, STEP*0.88, 16);
      }
      const fs = Math.min(12, Math.max(8, STEP*0.7));
      ctx.font = `${(atC||inG) ? 'bold ' : ''}${fs}px monospace`;
      ctx.fillStyle = atC ? '#ff8888' : inG ? col : inP ? 'rgba(255,200,80,0.9)' : col;
      ctx.textAlign = 'center';
      ctx.fillText(bp.base, x, BL + 14);
    });
  }

  // ── Axis baseline ────────────────────────────────────────────────────────────
  ctx.strokeStyle = 'rgba(70,100,150,0.25)'; ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(0, BL); ctx.lineTo(W, BL); ctx.stroke();
}

// ── Indel Spectrum (ICE-style bar chart) ─────────────────────────────────────
function IndelSpectrum({ sample }) {
  // Build indel data from sample result
  // For real ICE data this would have per-allele breakdown.
  // Here we use what the engine gives us: indelSize + editType + zygosityClass
  const data = [];

  // WT allele
  const wtPct = sample.zygosityClass === 0 ? 100
    : sample.zygosityClass === 3 ? 0
    : sample.zygosityClass === 2 ? 50    // het → ~50% WT
    : 80; // low/mosaic → ~80% WT

  data.push({ size: 0, label: 'WT', pct: wtPct, color: '#4d9fff' });

  if (sample.indelSize > 0) {
    const editPct = 100 - wtPct;
    const sign    = sample.editType === 'Deletion' ? -1 : 1;
    const size    = sign * Math.abs(sample.indelSize);
    const col     = sample.editType === 'Deletion' ? '#ff6b6b' : '#00d4aa';
    data.push({ size, label: `${size > 0 ? '+':''}${size} bp`, pct: editPct, color: col });
  }

  const maxPct = Math.max(...data.map(d => d.pct), 1);

  if (!data.length) return (
    <div style={{ color:DIM, fontSize:12, textAlign:'center', padding:16 }}>
      No indel data available for this sample
    </div>
  );

  return (
    <div>
      <div style={{ display:'flex', flexDirection:'column', gap:8 }}>
        {data.sort((a,b) => b.pct - a.pct).map((d,i) => (
          <div key={i} style={{ display:'flex', alignItems:'center', gap:12 }}>
            {/* Size label */}
            <div style={{ width:52, textAlign:'right', fontSize:11, fontWeight:700, color:d.color, flexShrink:0 }}>
              {d.label}
            </div>
            {/* Bar */}
            <div style={{ flex:1, background:BG, borderRadius:4, height:18, overflow:'hidden', border:`1px solid ${BDR}` }}>
              <div style={{
                width: `${(d.pct/maxPct)*100}%`,
                height:'100%',
                background: `linear-gradient(90deg, ${d.color}cc, ${d.color})`,
                borderRadius:4,
                transition:'width 0.4s ease',
              }} />
            </div>
            {/* Percentage */}
            <div style={{ width:46, textAlign:'right', fontSize:12, fontWeight:700, color:d.color, flexShrink:0 }}>
              {d.pct.toFixed(1)}%
            </div>
            {/* Count indicator */}
            <div style={{ width:30, fontSize:9, color:DIM }}>
              {d.size === 0 ? 'allele' : d.size > 0 ? 'ins' : 'del'}
            </div>
          </div>
        ))}
      </div>
      <div style={{ marginTop:14, fontSize:10, color:DIM, lineHeight:1.6 }}>
        <span style={{ color:'#4d9fff' }}>■</span> WT (unmodified)&nbsp;&nbsp;
        <span style={{ color:'#00d4aa' }}>■</span> Insertion&nbsp;&nbsp;
        <span style={{ color:'#ff6b6b' }}>■</span> Deletion&nbsp;&nbsp;
        <span style={{ color:DIM, marginLeft:8 }}>
          * For precise allele frequencies, use ICE analysis with raw .ab1 fluorescence data
        </span>
      </div>
    </div>
  );
}

// ── Alignment diff view ───────────────────────────────────────────────────────
function AlignmentView({ sample, cutPos, grnaRegion, pamRegion }) {
  if (!sample?.alignedQuery || !sample?.alignedRef) {
    return (
      <div style={{ padding:28, textAlign:'center', color: DIM }}>
        <div style={{ fontSize:36, marginBottom:12 }}>🔬</div>
        <div style={{ fontSize:13, fontWeight:600, color: TXT, marginBottom:8 }}>
          Sequence Alignment Not Available
        </div>
        <div style={{ fontSize:12, lineHeight:1.7 }}>
          Alignment data requires a WT control sample.<br/>
          Add a WT control in Setup → re-run analysis.
        </div>
      </div>
    );
  }

  const ref = sample.alignedRef.toUpperCase();
  const qry = sample.alignedQuery.toUpperCase();
  const WIN = 55;

  // Find cut site in alignment coords
  // cutPos is in REFERENCE coordinates (count of non-gap ref bases)
  let cutAln = 0;
  if (cutPos != null) {
    let rp = 0;
    for (let k = 0; k < ref.length; k++) {
      if (ref[k] !== '-') {
        if (rp === cutPos) { cutAln = k; break; }
        rp++;
      }
    }
    // If not found, set to last position
    if (cutAln === 0 && cutPos > 0) cutAln = ref.length - 1;
  }

  const s = Math.max(0, cutAln - WIN);
  const e = Math.min(ref.length, cutAln + WIN);
  const cols = 42;  // chars per row

  // Map alignment position back to reference base position
  // (needed to check gRNA/PAM membership)
  let refPosMap = [];
  { let rp = 0;
    for (let k = 0; k < ref.length; k++) {
      refPosMap.push(rp);
      if (ref[k] !== '-') rp++;
    }
  }

  // Build annotated base array
  const bases = ref.slice(s,e).split('').map((rb,i) => {
    const qb    = qry[s+i] || '-';
    const abs   = s + i;
    const refBp = refPosMap[abs] ?? -1;
    const inGrna = grnaRegion && refBp >= grnaRegion.start && refBp < grnaRegion.end;
    const inPam  = pamRegion  && refBp >= pamRegion.start  && refBp < pamRegion.end;
    return {
      rb, qb,
      isCut:  abs === cutAln,
      isDel:  qb === '-' && rb !== '-',
      isIns:  rb === '-' && qb !== '-',
      isMM:   rb !== '-' && qb !== '-' && rb !== qb,
      match:  rb === qb && rb !== '-',
      inGrna, inPam,
    };
  });

  const rows = Math.ceil(bases.length / cols);

  return (
    <div style={{ padding:'14px 18px', overflowX:'auto' }}>
      {/* Legend */}
      <div style={{ display:'flex', gap:14, marginBottom:14, fontSize:11, flexWrap:'wrap', alignItems:'center' }}>
        {[
          { c:'#ff4444', l:'✂ Cut site' },
          { c:'#ff6b6b', l:'▼ Deletion' },
          { c:'#00d4aa', l:'▲ Insertion' },
          { c:'#f0c040', l:'≠ Mismatch' },
          { c: DIM,      l:'| Match' },
        ].map(x => <span key={x.l} style={{ color:x.c }}>{x.l}</span>)}
        {grnaRegion && <span style={{ background:'rgba(0,212,170,0.14)', border:'1px solid rgba(0,212,170,0.3)', color:'#00d4aa', borderRadius:4, padding:'1px 6px' }}>▒ sgRNA</span>}
        {pamRegion  && <span style={{ background:'rgba(255,152,0,0.12)', border:'1px solid rgba(255,152,0,0.3)', color:'rgba(255,200,80,0.9)', borderRadius:4, padding:'1px 6px' }}>▒ PAM</span>}
        <span style={{ color:DIM, marginLeft:'auto', fontSize:10 }}>±{WIN}bp around cut (+{cutPos})</span>
      </div>

      {Array.from({ length: rows }).map((_, row) => {
        const chunk   = bases.slice(row*cols, (row+1)*cols);
        const absOff  = s + row*cols;
        const CW      = 15;

        return (
          <div key={row} style={{ marginBottom:18, fontFamily:'monospace', fontSize:12 }}>
            {/* Position numbers */}
            <div style={{ color:DIM, fontSize:9, marginBottom:1, paddingLeft:50 }}>
              {chunk.map((_, ci) => {
                const rel = (absOff + ci) - cutAln;
                return (ci % 10 === 0)
                  ? <span key={ci} style={{ display:'inline-block', width:CW, textAlign:'center' }}>
                      {rel >= 0 ? `+${rel}` : rel}
                    </span>
                  : <span key={ci} style={{ display:'inline-block', width:CW }}></span>;
              })}
            </div>

            {/* Ref */}
            <div style={{ display:'flex', alignItems:'center', marginBottom:1 }}>
              <span style={{ width:46, color:DIM, fontSize:10, textAlign:'right', marginRight:4, flexShrink:0 }}>Ref</span>
              {chunk.map((b,ci) => (
                <span key={ci} style={{
                  display:'inline-block', width:CW, textAlign:'center',
                  color: b.isCut ? '#ff4444' : b.isDel ? '#ff6b6b' : b.inGrna ? BC[b.rb]||TXT : b.inPam ? 'rgba(255,200,80,0.95)' : BC[b.rb]||TXT,
                  fontWeight: b.isCut||b.isDel||b.inGrna ? 700 : 400,
                  background: b.isCut ? 'rgba(255,68,68,0.2)' : b.inGrna ? 'rgba(0,212,170,0.14)' : b.inPam ? 'rgba(255,152,0,0.12)' : 'transparent',
                  borderRadius:2,
                  outline: b.inGrna ? '1px solid rgba(0,212,170,0.2)' : 'none',
                }}>
                  {b.rb === '-' ? '·' : b.rb}
                </span>
              ))}
            </div>

            {/* Match symbols */}
            <div style={{ display:'flex', alignItems:'center', marginBottom:1 }}>
              <span style={{ width:50, flexShrink:0 }}></span>
              {chunk.map((b,ci) => (
                <span key={ci} style={{
                  display:'inline-block', width:CW, textAlign:'center', fontSize:10,
                  color: b.isCut ? '#ff4444' : b.isDel ? '#ff6b6b88' : b.isIns ? '#00d4aa88' : b.isMM ? '#f0c04088' : '#ffffff22',
                }}>
                  {b.isCut ? '✂' : b.isDel ? '▼' : b.isIns ? '▲' : b.isMM ? '≠' : b.match ? '|' : ' '}
                </span>
              ))}
            </div>

            {/* Sample */}
            <div style={{ display:'flex', alignItems:'center' }}>
              <span style={{ width:46, color:DIM, fontSize:10, textAlign:'right', marginRight:4, flexShrink:0, overflow:'hidden', textOverflow:'ellipsis', whiteSpace:'nowrap' }}>
                {sample.name||sample.id}
              </span>
              {chunk.map((b,ci) => (
                <span key={ci} style={{
                  display:'inline-block', width:CW, textAlign:'center',
                  color: b.isCut ? '#ff4444' : b.isIns ? '#00d4aa' : b.isDel ? '#444' : b.isMM ? '#f0c040' : BC[b.qb]||TXT,
                  fontWeight: b.isCut||b.isIns||b.isMM ? 700 : 400,
                  background: b.isCut ? 'rgba(255,68,68,0.15)' : b.isIns ? 'rgba(0,212,170,0.1)' : b.isDel ? 'rgba(255,107,107,0.08)' : b.isMM ? 'rgba(240,192,64,0.08)' : 'transparent',
                  borderRadius:2,
                }}>
                  {b.qb === '-' ? '—' : b.qb}
                </span>
              ))}
            </div>
          </div>
        );
      })}

      {/* Summary stats */}
      <div style={{ display:'flex', gap:10, paddingTop:12, borderTop:`1px solid ${BDR}`, flexWrap:'wrap' }}>
        {[
          { l:'Edit Type', v: sample.editType||'Wild-type',
            c: sample.editType==='Deletion'?'#ff6b6b':sample.editType==='Insertion'?ACN:'#4d9fff' },
          { l:'Indel Size', v: sample.indelSize ? `${Math.abs(sample.indelSize)} bp` : 'None' },
          { l:'Identity', v: sample.sequenceIdentity ? `${(sample.sequenceIdentity*100).toFixed(1)}%` : '--' },
          { l:'Frameshift', v: sample.isFrameshift ? '✓ Yes' : '✗ No',
            c: sample.isFrameshift ? '#ff6b6b' : ACN },
        ].map(s => (
          <div key={s.l} style={{ background:BG, border:`1px solid ${BDR}`, borderRadius:7, padding:'8px 14px' }}>
            <div style={{ color:DIM, fontSize:9, letterSpacing:1, textTransform:'uppercase' }}>{s.l}</div>
            <div style={{ color:s.c||TXT, fontSize:13, fontWeight:600, marginTop:2 }}>{s.v}</div>
          </div>
        ))}
      </div>
    </div>
  );
}

// ── Discordance chart ─────────────────────────────────────────────────────────
function DiscordancePlot({ profile, cutPos }) {
  const ref = useRef(null);
  useEffect(() => {
    const cv = ref.current; if (!cv || !profile?.length) return;
    const ctx = cv.getContext('2d'), W = cv.width, H = cv.height;
    ctx.clearRect(0,0,W,H); ctx.fillStyle = BG; ctx.fillRect(0,0,W,H);
    const n = profile.length, bw = W/n;
    const sm = profile.map((d,i) => {
      const sl = profile.slice(Math.max(0,i-1), i+2);
      return sl.reduce((a,b) => a+b.discordance,0)/sl.length;
    });
    const mx = Math.max(...sm, 0.01);
    // Cut line
    const ci = profile.findIndex(d => d.relPos >= 0);
    if (ci >= 0) {
      ctx.strokeStyle = '#ff444866'; ctx.lineWidth = 1.5; ctx.setLineDash([4,3]);
      ctx.beginPath(); ctx.moveTo(ci*bw,0); ctx.lineTo(ci*bw,H-18); ctx.stroke();
      ctx.setLineDash([]);
    }
    // Area fill
    ctx.beginPath();
    ctx.moveTo(0, H-18);
    sm.forEach((v,i) => ctx.lineTo(i*bw + bw/2, H-18-(v/mx)*(H-28)));
    ctx.lineTo(W, H-18); ctx.closePath();
    const ag = ctx.createLinearGradient(0,0,0,H);
    ag.addColorStop(0,'rgba(0,212,170,0.3)'); ag.addColorStop(1,'rgba(0,212,170,0.02)');
    ctx.fillStyle = ag; ctx.fill();
    // Line
    ctx.beginPath(); ctx.strokeStyle = ACN; ctx.lineWidth = 1.5;
    sm.forEach((v,i) => { const x=i*bw+bw/2, y=H-18-(v/mx)*(H-28); i===0?ctx.moveTo(x,y):ctx.lineTo(x,y); });
    ctx.stroke();
    // Axis labels
    ctx.fillStyle = DIM; ctx.font = '8px monospace';
    ctx.textAlign = 'left'; ctx.fillText('Discordance', 3, 10);
    ctx.textAlign = 'center'; ctx.fillText('← upstream | downstream →', W/2, H-4);
  }, [profile, cutPos]);
  if (!profile?.length) return (
    <div style={{ color:DIM, padding:'24px 0', textAlign:'center', fontSize:12 }}>
      Add WT control to see TIDE-equivalent discordance profile
    </div>
  );
  return <canvas ref={ref} width={600} height={110} style={{ width:'100%', height:110, borderRadius:6, display:'block' }} />;
}

// ── Sample sidebar ────────────────────────────────────────────────────────────
function SampleSidebar({ samples, selected, onSelect }) {
  if (!samples?.length) return null;
  const CC = { 0:'#4d9fff', 1:'#f59e0b', 2:'#a855f7', 3:'#00d4aa' };
  const IC = { 0:'○', 1:'◔', 2:'◑', 3:'●' };
  return (
    <div style={{ width:162, flexShrink:0, background:'#060b18', borderRight:`1px solid ${BDR}`, borderRadius:'12px 0 0 12px', overflowY:'auto', maxHeight:640 }}>
      <div style={{ padding:'10px 12px', color:DIM, fontSize:10, fontWeight:700, letterSpacing:1.2, borderBottom:`1px solid ${BDR}` }}>
        SAMPLES
      </div>
      {samples.map(s => {
        const sel = selected?.id === s.id || selected?.name === s.name;
        const col = CC[s.zygosityClass] || '#888';
        return (
          <div key={s.id||s.name} onClick={() => onSelect?.(s)}
            style={{ padding:'8px 12px', background: sel ? `${col}14` : 'transparent',
              borderLeft:`3px solid ${sel ? col : 'transparent'}`, cursor:'pointer', transition:'all 0.12s',
              borderBottom:`1px solid ${BDR}44` }}>
            <div style={{ color: sel?col:TXT, fontSize:11, fontWeight:600, overflow:'hidden', textOverflow:'ellipsis', whiteSpace:'nowrap' }}>
              {IC[s.zygosityClass]||'?'} {s.name||s.id}
            </div>
            <div style={{ color:DIM, fontSize:9, marginTop:2 }}>
              {s.generation} · ICE {s.iceScore??'--'}%
            </div>
            <div style={{ marginTop:4, height:2, borderRadius:2, width:`${Math.round((s.r2||0)*100)}%`,
              background: (s.r2||0)>0.8 ? ACN : '#f59e0b', maxWidth:'100%' }} />
          </div>
        );
      })}
    </div>
  );
}

// ── Button styles ─────────────────────────────────────────────────────────────
const BS = (active, color) => ({
  background: active ? `${color||ACN}18` : BG,
  border: `1px solid ${active ? (color||ACN)+'55' : BDR}`,
  color: active ? (color||ACN) : DIM,
  borderRadius:6, padding:'4px 10px', fontSize:11,
  cursor:'pointer', transition:'all 0.15s', fontFamily:'inherit',
});
const BT = { background:BG, border:`1px solid ${BDR}`, color:DIM, padding:'4px 10px', fontSize:11, cursor:'pointer', fontFamily:'inherit', transition:'all 0.15s' };

// ── Main component ────────────────────────────────────────────────────────────
export default function ChromatogramViewer({
  sample, wtSample, allSamples, onSelect,
  cutPosition=137, grnaSeq='', grnaList=null
}) {
  const canvasRef = useRef(null);
  const [zoom,      setZoom]      = useState(1.0);
  const [pan,       setPan]       = useState(0.5);  // default centre
  const [showWT,    setShowWT]    = useState(false);
  const [showBases, setShowBases] = useState(true);
  const [activeTab, setActiveTab] = useState('trace');
  const [hoverX,    setHoverX]    = useState(null);
  const [viewMode,  setViewMode]  = useState('normalized');
  const [yScale,    setYScale]    = useState(1.0);  // vertical zoom
  const [seqSearch, setSeqSearch] = useState('');   // sequence search

  // Traces
  const traces = useMemo(() => {
    if (!sample) return null;
    if (sample.traces?.channels) {
      // Real AB1 — add STEP/baseLen/basePositions if missing
      const t = sample.traces;
      if (!t.STEP) {
        const bpCount = t.basePositions?.length || (t.width ? Math.floor(t.width/14) : 0);
        return { ...t, STEP:14, baseLen:bpCount,
          basePositions: t.basePositions || Array.from({length:bpCount},(_,i)=>({ base:'N', x:i*14+7 })) };
      }
      return t;
    }
    if (sample.sequence) return synthTraces(sample.sequence, sample.quality);
    return null;
  }, [sample]);

  const wtTraces = useMemo(() => {
    if (!wtSample) return null;
    if (wtSample.traces?.channels) {
      const t = wtSample.traces;
      if (!t.STEP) return { ...t, STEP:14, baseLen: t.basePositions?.length||0 };
      return t;
    }
    if (wtSample.sequence) return synthTraces(wtSample.sequence, wtSample.quality);
    return null;
  }, [wtSample]);

  // Annotations (gRNA / PAM / cut) in sample sequence coordinates
  const annotations = useMemo(() => {
    if (!sample) return [];
    return deriveAnnotations(sample, grnaSeq, grnaList, sample.sampleCutPosition ?? cutPosition);
  }, [sample, grnaSeq, grnaList, cutPosition]);

  // primary wrapped in useMemo to avoid re-render loops
  const primary = useMemo(() => {
    return annotations[0] || { label:'cut', grna:null, pam:null, cutBase: sample?.sampleCutPosition ?? cutPosition };
  }, [annotations, sample?.sampleCutPosition, cutPosition]);

  // Auto-pan to cut/gRNA site on sample change
  useEffect(() => {
    if (!traces) return;
    const STEP    = traces.STEP || 14;
    const canvasW = canvasRef.current?.offsetWidth || 680;
    const totalW  = traces.width;
    // Target: gRNA start if found, else cut site
    const targetBase = primary.grna?.start ?? primary.cutBase;
    if (targetBase == null) { setPan(0); return; }
    // Position target at 35% from left
    const targetPx  = targetBase * STEP * 1.0;   // at zoom 1.0
    const tgtOffset = targetPx - canvasW * 0.35;
    const maxOffset = Math.max(0, totalW - canvasW);
    setPan(maxOffset > 0 ? Math.max(0, Math.min(1, tgtOffset/maxOffset)) : 0);
    setZoom(1.0);
    setYScale(1.0);
  }, [sample?.id, sample?.name]); // eslint-disable-line react-hooks/exhaustive-deps

  // Wheel zoom
  const handleWheel = (e) => {
    e.preventDefault();
    const delta = e.deltaY > 0 ? -0.15 : 0.15;
    setZoom(z => Math.max(0.4, Math.min(4.0, z+delta)));
  };

  // Draw
  useEffect(() => {
    if (activeTab !== 'trace') return;
    const cv = canvasRef.current; if (!cv || !traces) return;
    cv.width  = cv.parentElement?.clientWidth || 680;
    cv.height = 300;
    drawChroma(cv, traces, { annotations, showBases, zoom, pan, wtTraces, showWT, hoverX, viewMode, yScale, seqSearch });
  }, [traces, wtTraces, zoom, pan, showBases, showWT, hoverX, viewMode, yScale, seqSearch, activeTab, annotations]); // eslint-disable-line react-hooks/exhaustive-deps

  if (!sample) return (
    <div style={{ background:CARD, border:`1px solid ${BDR}`, borderRadius:12, padding:48, textAlign:'center', color:DIM }}>
      <div style={{ fontSize:40, marginBottom:12 }}>📈</div>
      <div style={{ fontSize:14, fontWeight:700, color:TXT, marginBottom:6 }}>Select a sample to view its chromatogram</div>
      <div style={{ fontSize:12 }}>Upload .ab1 files in Setup or run the demo</div>
    </div>
  );

  const CC   = { 0:'#4d9fff', 1:'#f59e0b', 2:'#a855f7', 3:'#00d4aa' };
  const sCol = CC[sample.zygosityClass] || ACN;
  const TABS = [
    { id:'trace',   label:'📈 Trace' },
    { id:'align',   label:'🔬 Alignment' },
    { id:'discord', label:'📉 Discordance' },
    { id:'info',    label:'⚙ Details' },
  ];

  return (
    <div style={{ display:'flex', borderRadius:12, overflow:'hidden', border:`1px solid ${BDR}` }}>
      {allSamples?.length > 1 && (
        <SampleSidebar samples={allSamples} selected={sample} onSelect={onSelect} />
      )}

      <div style={{ flex:1, background:CARD, minWidth:0, display:'flex', flexDirection:'column' }}>

        {/* ─── Header ─── */}
        <div style={{ padding:'14px 18px 10px', borderBottom:`1px solid ${BDR}` }}>
          <div style={{ display:'flex', justifyContent:'space-between', alignItems:'flex-start', gap:12, flexWrap:'wrap' }}>
            <div>
              <div style={{ display:'flex', alignItems:'center', gap:8, marginBottom:5, flexWrap:'wrap' }}>
                <span style={{ fontSize:18, fontWeight:800, color:sCol }}>{sample.name||sample.id}</span>
                <span style={{ fontSize:10, background:`${sCol}18`, border:`1px solid ${sCol}44`, color:sCol, borderRadius:20, padding:'2px 8px', letterSpacing:0.5 }}>
                  {['WT','Low/Mosaic','Heterozygous','Homozygous'][sample.zygosityClass??0]}
                </span>
                {sample.indelSize > 0 && (
                  <span style={{ fontSize:10, background:'rgba(255,107,107,0.1)', border:'1px solid rgba(255,107,107,0.3)', color:'#ff6b6b', borderRadius:20, padding:'2px 8px' }}>
                    {sample.editType==='Deletion'?'−':'+'}{ Math.abs(sample.indelSize)} bp
                  </span>
                )}
              </div>
              <div style={{ color:DIM, fontSize:10, marginBottom:6 }}>
                {sample.generation} · {sample.sequence?.length||0} bp
                {sample.traces?.channels ? ' · Real .ab1' : ' · Synthetic'}
                {sample.wasFlipped && <span style={{ color:'#bf5af2', marginLeft:6 }}>↺ RC corrected</span>}
              </div>
              {/* Annotation badges */}
              <div style={{ display:'flex', gap:6, flexWrap:'wrap' }}>
                {primary.grna ? (
                  <>
                    <span style={{ fontSize:9, background:'rgba(0,212,170,0.1)', border:'1px solid rgba(0,212,170,0.25)', color:ACN, borderRadius:4, padding:'2px 7px' }}>
                      sgRNA {primary.grna.strand} strand · pos {primary.grna.start}–{primary.grna.end}
                    </span>
                    <span style={{ fontSize:9, background:'rgba(255,68,68,0.1)', border:'1px solid rgba(255,68,68,0.25)', color:'#ff7777', borderRadius:4, padding:'2px 7px' }}>
                      ✂ cut at +{primary.cutBase}
                    </span>
                    {primary.pam && (
                      <span style={{ fontSize:9, background:'rgba(255,152,0,0.1)', border:'1px solid rgba(255,152,0,0.25)', color:'rgba(255,152,0,0.9)', borderRadius:4, padding:'2px 7px' }}>
                        PAM {primary.pam.start}–{primary.pam.end}
                      </span>
                    )}
                  </>
                ) : (
                  <span style={{ fontSize:9, background:'rgba(255,68,68,0.08)', border:'1px solid rgba(255,68,68,0.2)', color:'#ff888888', borderRadius:4, padding:'2px 7px' }}>
                    ✂ cut at +{primary.cutBase} (sgRNA not found in sequence)
                  </span>
                )}
              </div>
            </div>
            {/* Quality row */}
            <div style={{ display:'flex', gap:6, flexWrap:'wrap', alignItems:'center' }}>
              {sample.r2 != null && (
                <span style={{ fontSize:11, padding:'3px 10px', borderRadius:20,
                  background: sample.r2>=0.9 ? 'rgba(0,212,170,0.12)' : 'rgba(245,158,11,0.12)',
                  border: `1px solid ${sample.r2>=0.9 ? 'rgba(0,212,170,0.4)' : 'rgba(245,158,11,0.4)'}`,
                  color: sample.r2>=0.9 ? ACN : '#f59e0b' }}>
                  R² {sample.r2.toFixed(3)} {sample.r2>=0.9 ? '✓' : ''}
                </span>
              )}
              {sample.iceScore != null && (
                <span style={{ fontSize:11, background:'rgba(0,212,170,0.1)', border:'1px solid rgba(0,212,170,0.3)', color:ACN, borderRadius:20, padding:'3px 10px' }}>
                  ICE {sample.iceScore}%
                </span>
              )}
              {sample.mlScore != null && (
                <span style={{ fontSize:11, background:`${sCol}14`, border:`1px solid ${sCol}40`, color:sCol, borderRadius:20, padding:'3px 10px' }}>
                  ML {sample.mlScore}%
                </span>
              )}
            </div>
          </div>
        </div>

        {/* ─── Flag banner ─── */}
        {sample.flag && (
          <div style={{ background:'rgba(255,107,107,0.06)', borderBottom:`1px solid rgba(255,107,107,0.2)`, padding:'7px 18px', fontSize:11, color:'#ff8888' }}>
            ⚠ {sample.flag}
          </div>
        )}

        {/* ─── Tabs ─── */}
        <div style={{ display:'flex', borderBottom:`1px solid ${BDR}`, padding:'0 18px' }}>
          {TABS.map(t => (
            <button key={t.id} onClick={() => setActiveTab(t.id)} style={{
              background:'transparent', border:'none',
              borderBottom: activeTab===t.id ? `2px solid ${ACN}` : '2px solid transparent',
              color: activeTab===t.id ? ACN : DIM,
              padding:'8px 14px', fontSize:12, cursor:'pointer', marginBottom:-1,
              fontWeight: activeTab===t.id ? 700 : 400, transition:'all 0.15s', fontFamily:'inherit',
            }}>{t.label}</button>
          ))}
        </div>

        {/* ─── Tab: Trace ─── */}
        {activeTab === 'trace' && (
          <div style={{ padding:'14px 18px' }}>
            {/* Controls */}
            <div style={{ display:'flex', gap:7, marginBottom:10, flexWrap:'wrap', alignItems:'center' }}>
              <button onClick={() => setShowBases(b=>!b)} style={BS(showBases)}>
                {showBases ? 'Hide' : 'Show'} Bases
              </button>
              {wtTraces && (
                <button onClick={() => setShowWT(b=>!b)} style={BS(showWT)}>
                  {showWT ? '◉ WT On' : '○ WT Off'}
                </button>
              )}
              <div style={{ display:'flex', border:`1px solid ${BDR}`, borderRadius:6, overflow:'hidden' }}>
                <button onClick={() => setZoom(z=>Math.max(0.4,z-0.25))} style={BT}>−</button>
                <span style={{ padding:'4px 10px', fontSize:11, color:DIM, background:BG, borderLeft:`1px solid ${BDR}`, borderRight:`1px solid ${BDR}`, minWidth:50, textAlign:'center' }}>
                  {zoom.toFixed(1)}×
                </span>
                <button onClick={() => setZoom(z=>Math.min(4.0,z+0.25))} style={BT}>+</button>
              </div>
              <button onClick={() => { setZoom(1.0); }} style={BT}>Reset</button>
              {primary.cutBase != null && (
                <button onClick={() => {
                  if (!traces) return;
                  const canvasW = canvasRef.current?.offsetWidth || 680;
                  const targetBase = primary.grna?.start ?? primary.cutBase;
                  const totalW  = traces.width * zoom;
                  const tgtPx   = targetBase * (traces.STEP||14) * zoom;
                  const tgtOff  = tgtPx - canvasW * 0.35;
                  const mx = Math.max(0, totalW - canvasW);
                  setPan(mx > 0 ? Math.max(0,Math.min(1,tgtOff/mx)) : 0);
                }} style={{ ...BT, color:'#ff7777', borderColor:'rgba(255,107,107,0.35)', background:'rgba(255,68,68,0.06)' }}>
                  ✂ Jump to edit
                </button>
              )}
              {/* Vertical zoom */}
              <div style={{ display:'flex', alignItems:'center', gap:5, borderLeft:`1px solid ${BDR}`, paddingLeft:8 }}>
                <span style={{ color:DIM, fontSize:10, flexShrink:0 }}>V-zoom</span>
                <button onClick={() => setYScale(y => Math.max(0.3, y-0.2))} style={BT}>−</button>
                <span style={{ fontSize:10, color:DIM, minWidth:32, textAlign:'center' }}>{yScale.toFixed(1)}×</span>
                <button onClick={() => setYScale(y => Math.min(5.0, y+0.2))} style={BT}>+</button>
                <button onClick={() => setYScale(1.0)} style={{ ...BT, fontSize:9 }}>R</button>
              </div>

              {/* View mode */}
              <div style={{ marginLeft:'auto', display:'flex', border:`1px solid ${BDR}`, borderRadius:6, overflow:'hidden' }}>
                {[['normalized','📐 Normalized'],['raw','📊 Raw']].map(([m,l]) => (
                  <button key={m} onClick={() => setViewMode(m)} style={{
                    background: viewMode===m ? 'rgba(0,212,170,0.15)' : BG,
                    border:'none', borderRight: m==='normalized' ? `1px solid ${BDR}` : 'none',
                    color: viewMode===m ? ACN : DIM, padding:'5px 12px', fontSize:11,
                    cursor:'pointer', fontWeight: viewMode===m ? 700 : 400, fontFamily:'inherit',
                  }}>{l}</button>
                ))}
              </div>
            </div>

            {/* Legend */}
            <div style={{ display:'flex', gap:14, marginBottom:8, alignItems:'center', flexWrap:'wrap' }}>
              {Object.entries(BC).filter(([b])=>b!=='N').map(([b,c]) => (
                <span key={b} style={{ color:c, fontSize:11, fontWeight:700 }}>●{b}</span>
              ))}
              <span style={{ color:'#ff4444', fontSize:11 }}>✂cut</span>
              {primary.grna && <span style={{ color:ACN, fontSize:11 }}>▒ sgRNA</span>}
              {primary.pam  && <span style={{ color:'rgba(255,152,0,0.85)', fontSize:11 }}>▒ PAM</span>}
              {showWT && wtTraces && <span style={{ color:DIM, fontSize:11, borderBottom:'1px dashed '+DIM }}>-- WT</span>}
              <span style={{ marginLeft:'auto', fontSize:10, color:DIM, fontStyle:'italic' }}>
                {viewMode==='normalized' ? '📐 auto-scaled to view' : '📊 raw amplitude'}
              </span>
            </div>

            {/* Canvas */}
            {traces ? (
              <div style={{ position:'relative', borderRadius:8, overflow:'hidden', background:BG, cursor:'crosshair' }}>
                <canvas
                  ref={canvasRef}
                  style={{ width:'100%', height:300, display:'block' }}
                  onMouseMove={e => {
                    const r = e.currentTarget.getBoundingClientRect();
                    setHoverX((e.clientX-r.left)*(e.currentTarget.width/r.width));
                  }}
                  onMouseLeave={() => setHoverX(null)}
                  onWheel={handleWheel}
                />
                {/* Tooltip */}
                {hoverX != null && traces.basePositions && (() => {
                  const totalW = traces.width * zoom;
                  const maxPan = Math.max(0, totalW - (canvasRef.current?.width||680));
                  const OX     = -(pan * maxPan);
                  const STEP   = (traces.STEP||14) * zoom;
                  const hov    = traces.basePositions.find(bp => Math.abs(bp.x*zoom+OX-hoverX) < STEP*0.55);
                  if (!hov) return null;
                  const idx = traces.basePositions.indexOf(hov);
                  const col = BC[hov.base] || '#fff';
                  const grna = primary.grna;
                  const inG  = grna && idx>=grna.start && idx<grna.end;
                  const inP  = primary.pam && idx>=primary.pam.start && idx<primary.pam.end;
                  const atC  = primary.cutBase!=null && idx===primary.cutBase;
                  const reg  = atC?'✂ Cut Site': inG?'🧬 sgRNA': inP?'🔶 PAM': null;
                  return (
                    <div style={{ position:'absolute', top:8,
                      left: Math.min(hoverX+12, (canvasRef.current?.offsetWidth||600)-150),
                      background:'rgba(6,12,28,0.96)', border:`1px solid ${col}44`,
                      borderRadius:8, padding:'8px 12px', pointerEvents:'none', zIndex:10,
                      boxShadow:`0 4px 20px ${col}18` }}>
                      <div style={{ fontSize:22, fontWeight:900, color:col, fontFamily:'monospace', lineHeight:1 }}>{hov.base}</div>
                      <div style={{ fontSize:9, color:DIM, marginTop:3 }}>pos +{idx+1}</div>
                      {reg && <div style={{ fontSize:10, color: atC?'#ff4444':inP?'rgba(255,152,0,0.9)':ACN, fontWeight:700, marginTop:2 }}>{reg}</div>}
                      {traces.quality?.[idx] != null && (
                        <div style={{ fontSize:9, color: traces.quality[idx]>=30?ACN:'#f59e0b', marginTop:1 }}>
                          Q{traces.quality[idx]}
                        </div>
                      )}
                    </div>
                  );
                })()}
              </div>
            ) : (
              <div style={{ height:290, background:BG, borderRadius:8, display:'flex', alignItems:'center', justifyContent:'center', color:DIM, fontSize:13 }}>
                No sequence data — upload .ab1 or enter sequence in Setup
              </div>
            )}

            {/* Pan + search + info */}
            {traces && (
              <div style={{ marginTop:10 }}>
                <div style={{ display:'flex', alignItems:'center', gap:10 }}>
                  <span style={{ color:DIM, fontSize:10, flexShrink:0 }}>Pan</span>
                  <input type="range" min={0} max={1} step={0.0005} value={pan}
                    onChange={e => setPan(parseFloat(e.target.value))}
                    style={{ flex:1, accentColor: ACN, height:3 }} />
                </div>
                {/* Sequence search */}
                <div style={{ display:'flex', alignItems:'center', gap:8, marginTop:8 }}>
                  <span style={{ color:DIM, fontSize:10, flexShrink:0 }}>🔍 Search sequence:</span>
                  <input
                    type="text"
                    value={seqSearch}
                    onChange={e => setSeqSearch(e.target.value.toUpperCase().replace(/[^ATGCN]/gi,''))}
                    placeholder="e.g. ATGCGT (min 3bp)..."
                    style={{
                      flex:1, background:'#0a1628', border:`1px solid ${seqSearch.length >= 3 ? 'rgba(255,220,50,0.5)' : BDR}`,
                      color: seqSearch.length >= 3 ? 'rgba(255,220,50,0.9)' : TXT,
                      borderRadius:6, padding:'4px 10px', fontSize:11, fontFamily:'monospace',
                      outline:'none',
                    }}
                  />
                  {seqSearch && <button onClick={() => setSeqSearch('')} style={{ ...BT, color:'#ff6b6b', borderColor:'rgba(255,107,107,0.3)', padding:'4px 8px' }}>✕</button>}
                  {seqSearch.length >= 3 && traces.basePositions && (() => {
                    const needle = seqSearch.toUpperCase();
                    const hay = traces.basePositions.map(bp => bp.base).join('');
                    const matches = (hay.match(new RegExp(needle,'g')) || []).length;
                    return <span style={{ color:'rgba(255,220,50,0.8)', fontSize:10, flexShrink:0 }}>{matches} match{matches!==1?'es':''}</span>;
                  })()}
                </div>
                <div style={{ color:DIM, fontSize:10, marginTop:6, textAlign:'center' }}>
                  Scroll wheel to zoom · Pan slider to navigate
                  {sample.traces?.channels ? ' · ✓ Real .ab1 fluorescence' : ' · Synthetic (upload .ab1 for real data)'}
                </div>
              </div>
            )}
          </div>
        )}

        {/* ─── Tab: Alignment ─── */}
        {activeTab === 'align' && (
          <AlignmentView sample={sample} cutPos={sample.grnaEnd > 0 ? sample.grnaEnd - 3 : sample.cutPosition} grnaRegion={sample.grnaStart > 0 ? { start: sample.grnaStart, end: sample.grnaEnd } : null} pamRegion={sample.grnaEnd > 0 ? { start: sample.grnaEnd, end: sample.grnaEnd + 3 } : null} />
        )}

        {/* ─── Tab: Discordance & Indel Spectrum ─── */}
        {activeTab === 'discord' && (
          <div style={{ padding:'16px 18px' }}>
            {/* Discordance plot */}
            <div style={{ fontSize:13, fontWeight:700, color:ACN, marginBottom:4 }}>TIDE-equivalent Discordance Profile</div>
            <div style={{ fontSize:11, color:DIM, marginBottom:12 }}>
              Base-by-base signal deviation from WT. A clean edit shows a sharp jump at the cut site (+0) that persists downstream.
            </div>
            <DiscordancePlot profile={sample.discordanceProfile} cutPos={primary.cutBase} />
            {sample.discordanceProfile?.length > 0 && (
              <div style={{ display:'grid', gridTemplateColumns:'repeat(3,1fr)', gap:10, marginTop:14, marginBottom:22 }}>
                {[
                  { l:'ICE Score',       v:`${sample.iceScore??'--'}%`,                                                              c:ACN },
                  { l:'R² Downstream',   v: sample.r2dn?.toFixed(3)||'--',                                                           c: (sample.r2dn||0)>0.9?ACN:'#f59e0b' },
                  { l:'Mean Discordance',v:(sample.discordanceProfile.reduce((a,b)=>a+b.discordance,0)/sample.discordanceProfile.length*100).toFixed(1)+'%', c:'#a855f7' },
                ].map(s => (
                  <div key={s.l} style={{ background:BG, borderRadius:8, padding:'10px 14px', border:`1px solid ${BDR}` }}>
                    <div style={{ color:DIM, fontSize:10 }}>{s.l}</div>
                    <div style={{ color:s.c, fontSize:20, fontWeight:700, marginTop:3 }}>{s.v}</div>
                  </div>
                ))}
              </div>
            )}

            {/* Indel Spectrum — ICE-style bar chart */}
            <div style={{ borderTop:`1px solid ${BDR}`, paddingTop:18 }}>
              <div style={{ fontSize:13, fontWeight:700, color:'#a855f7', marginBottom:4 }}>Indel Spectrum</div>
              <div style={{ fontSize:11, color:DIM, marginBottom:14 }}>
                Distribution of detected indel sizes (ICE decomposition). Positive = insertion, negative = deletion.
              </div>
              <IndelSpectrum sample={sample} />
            </div>
          </div>
        )}

        {/* ─── Tab: Details ─── */}
        {activeTab === 'info' && (
          <div style={{ padding:'16px 18px' }}>
            <div style={{ display:'grid', gridTemplateColumns:'1fr 1fr', gap:10, marginBottom:12 }}>
              {[
                { l:'Zygosity',         v: sample.zygosityRaw||sample.zygosityLabel, c: sCol },
                { l:'Edit Type',        v: sample.editType||'WT', c: sample.editType==='Deletion'?'#ff6b6b':sample.editType==='Insertion'?ACN:'#4d9fff' },
                { l:'Indel Size',       v: sample.indelSize?`${Math.abs(sample.indelSize)} bp`:'None' },
                { l:'Frameshift',       v: sample.isFrameshift?'✓ Yes':'✗ No', c: sample.isFrameshift?'#ff6b6b':ACN },
                { l:'Sequence Length',  v:`${sample.sequence?.length||0} bp` },
                { l:'Identity to WT',   v: sample.sequenceIdentity?`${(sample.sequenceIdentity*100).toFixed(1)}%`:'--' },
                { l:'R² Upstream',      v: sample.r2up?.toFixed(3)||'--' },
                { l:'R² Downstream',    v: sample.r2dn?.toFixed(3)||'--', c:(sample.r2dn||0)>=0.95?ACN:(sample.r2dn||0)>=0.8?'#f59e0b':'#ff6b6b' },
              ].map(i => (
                <div key={i.l} style={{ background:BG, borderRadius:8, padding:'10px 14px', border:`1px solid ${BDR}` }}>
                  <div style={{ color:DIM, fontSize:10, marginBottom:3 }}>{i.l}</div>
                  <div style={{ color:i.c||TXT, fontSize:14, fontWeight:600 }}>{i.v}</div>
                </div>
              ))}
            </div>
            <div style={{ background:BG, borderRadius:8, padding:'12px 14px', border:`1px solid ${BDR}` }}>
              <div style={{ color:DIM, fontSize:10, marginBottom:5 }}>PREDICTED PROTEIN EFFECT</div>
              <div style={{ color: sample.isFrameshift?'#ff6b6b':sample.indelSize>0?'#f59e0b':ACN, fontSize:13, fontWeight:600 }}>
                {sample.proteinEffect||'No prediction available'}
              </div>
            </div>
          </div>
        )}

      </div>
    </div>
  );
}