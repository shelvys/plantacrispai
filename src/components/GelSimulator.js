// ============================================================
// FILE: src/components/GelSimulator.js
// PURPOSE: Visual agarose gel electrophoresis simulator
//
// Renders a realistic-looking gel image showing:
//   - DNA ladder lane (1kb, 750, 500, 400, 300, 250, 200, 150, 100)
//   - Wild-type control lane (2 bands — site cut)
//   - Each sample lane (1, 2 or 3 bands depending on zygosity)
//
// Band positions are calculated from real gel physics:
//   - Band migration ∝ log(size) at constant voltage
//   - Larger fragments migrate slower (higher on gel)
//   - Band thickness ∝ DNA amount (uniform here)
// ============================================================

import React, { useRef, useEffect, useCallback } from 'react';

// Reference ladder sizes for scale
const LADDER_SIZES = [1000, 750, 500, 400, 300, 250, 200, 150, 100];

// Lane color palette
const LANE_COLORS = {
  ladder:     '#8899BB',
  wt_control: '#4D9FFF',
  wt:         '#4D9FFF',
  edited:     '#FF4D6D',
  mixed:      '#F0C040',
};

// Gel theme
const GEL = {
  background:    '#0A0F1A',
  gelBackground: '#080C18',
  gelBorder:     '#1a3060',
  laneDiv:       'rgba(26,48,96,0.4)',
  bandGlow:      0.6,
  textColor:     '#8899BB',
};

/**
 * Convert band size (bp) to gel Y position
 * Gel migration follows logarithmic relationship with size
 * y = a - b * log(size)
 *
 * @param {number} sizeInBp - band size in base pairs
 * @param {number} gelHeight - total gel height in pixels
 * @param {number} topMargin - top margin for labels
 * @param {number} bottomMargin - bottom margin
 * @param {number} maxSize - largest fragment (reference for top)
 * @param {number} minSize - smallest fragment (reference for bottom)
 */
function bpToY(sizeInBp, gelHeight, topMargin, bottomMargin, maxSize = 1200, minSize = 80) {
  const logMax = Math.log10(maxSize);
  const logMin = Math.log10(minSize);
  const logSize = Math.log10(Math.max(sizeInBp, minSize));

  // Normalize: 0 = top (large fragments), 1 = bottom (small fragments)
  const norm = (logMax - logSize) / (logMax - logMin);

  const drawHeight = gelHeight - topMargin - bottomMargin;
  return topMargin + norm * drawHeight;
}

export default function GelSimulator({ lanes, enzyme, ampliconLength }) {
  const canvasRef = useRef(null);

  const draw = useCallback(() => {
    const canvas = canvasRef.current;
    if (!canvas || !lanes || lanes.length === 0) return;

    const ctx    = canvas.getContext('2d');
    const W      = canvas.width;
    const H      = canvas.height;

    // Layout constants
    const TOP_MARGIN    = 48;  // space for lane labels
    const BOTTOM_MARGIN = 36;  // space for base labels
    const LEFT_MARGIN   = 50;  // space for bp scale
    const RIGHT_MARGIN  = 16;
    const GEL_W         = W - LEFT_MARGIN - RIGHT_MARGIN;
    const GEL_H         = H - TOP_MARGIN - BOTTOM_MARGIN;
    const LANE_W        = GEL_W / lanes.length;
    const BAND_H        = 6;
    const MAX_SIZE      = 1200;
    const MIN_SIZE      = 80;

    // Clear
    ctx.clearRect(0, 0, W, H);
    ctx.fillStyle = GEL.background;
    ctx.fillRect(0, 0, W, H);

    // ── Gel background ────────────────────────────────────
    ctx.fillStyle = GEL.gelBackground;
    ctx.fillRect(LEFT_MARGIN, TOP_MARGIN, GEL_W, GEL_H);
    ctx.strokeStyle = GEL.gelBorder;
    ctx.lineWidth   = 1;
    ctx.strokeRect(LEFT_MARGIN, TOP_MARGIN, GEL_W, GEL_H);

    // ── Horizontal reference lines (dashed) ──────────────
    ctx.setLineDash([3, 6]);
    ctx.strokeStyle = 'rgba(26,48,96,0.6)';
    ctx.lineWidth   = 0.8;
    LADDER_SIZES.forEach(size => {
      const y = bpToY(size, H, TOP_MARGIN, BOTTOM_MARGIN, MAX_SIZE, MIN_SIZE);
      ctx.beginPath();
      ctx.moveTo(LEFT_MARGIN, y);
      ctx.lineTo(LEFT_MARGIN + GEL_W, y);
      ctx.stroke();
    });
    ctx.setLineDash([]);

    // ── BP scale (left axis) ─────────────────────────────
    ctx.fillStyle  = GEL.textColor;
    ctx.font       = '9px monospace';
    ctx.textAlign  = 'right';
    LADDER_SIZES.forEach(size => {
      const y = bpToY(size, H, TOP_MARGIN, BOTTOM_MARGIN, MAX_SIZE, MIN_SIZE);
      ctx.fillText(`${size}`, LEFT_MARGIN - 6, y + 4);
    });

    // ── Lane dividers ─────────────────────────────────────
    ctx.strokeStyle = GEL.laneDiv;
    ctx.lineWidth   = 0.5;
    for (let i = 1; i < lanes.length; i++) {
      const x = LEFT_MARGIN + i * LANE_W;
      ctx.beginPath();
      ctx.moveTo(x, TOP_MARGIN);
      ctx.lineTo(x, TOP_MARGIN + GEL_H);
      ctx.stroke();
    }

    // ── Draw each lane ────────────────────────────────────
    lanes.forEach((lane, laneIdx) => {
      const laneX    = LEFT_MARGIN + laneIdx * LANE_W;
      const laneCX   = laneX + LANE_W / 2;
      const color    = LANE_COLORS[lane.type] || '#8899BB';
      const bandW    = Math.max(8, LANE_W * 0.55);

      // Lane label (top)
      ctx.fillStyle = lane.type === 'ladder' ? GEL.textColor : color;
      ctx.font      = `bold 9px monospace`;
      ctx.textAlign = 'center';
      const labelText = lane.label.length > 10
        ? lane.label.slice(0, 9) + '…'
        : lane.label;
      ctx.fillText(labelText, laneCX, TOP_MARGIN - 20);

      // Lane type indicator
      if (lane.note) {
        ctx.font      = '8px monospace';
        ctx.fillStyle = GEL.textColor;
        ctx.fillText(lane.note, laneCX, TOP_MARGIN - 10);
      }

      // Draw bands
      (lane.bands || []).forEach(size => {
        if (!size || size <= 0) return;

        const y = bpToY(size, H, TOP_MARGIN, BOTTOM_MARGIN, MAX_SIZE, MIN_SIZE);

        // Band glow
        const gradient = ctx.createLinearGradient(
          laneCX - bandW/2, y - BAND_H,
          laneCX + bandW/2, y + BAND_H
        );
        gradient.addColorStop(0, `${color}00`);
        gradient.addColorStop(0.3, `${color}CC`);
        gradient.addColorStop(0.5, `${color}FF`);
        gradient.addColorStop(0.7, `${color}CC`);
        gradient.addColorStop(1, `${color}00`);

        // Outer glow
        ctx.shadowColor = color;
        ctx.shadowBlur  = 8;
        ctx.fillStyle   = gradient;
        ctx.fillRect(laneCX - bandW/2, y - BAND_H/2, bandW, BAND_H);
        ctx.shadowBlur  = 0;

        // Band size label (right of band for ladder, below for others)
        ctx.font       = '8px monospace';
        ctx.textAlign  = laneIdx === 0 ? 'left' : 'center';
        ctx.fillStyle  = lane.type === 'ladder' ? GEL.textColor : `${color}BB`;

        if (lane.type === 'ladder') {
          ctx.fillText(`${size}`, laneCX + bandW/2 + 3, y + 3);
        } else {
          ctx.fillText(`${size}bp`, laneCX, y + BAND_H + 10);
        }
      });

      // Well indicator at top
      ctx.fillStyle = 'rgba(26,48,96,0.8)';
      ctx.fillRect(laneCX - 8, TOP_MARGIN, 16, 6);
      ctx.fillStyle = 'rgba(0,212,170,0.3)';
      ctx.fillRect(laneCX - 7, TOP_MARGIN + 1, 14, 4);
    });

    // ── Title ─────────────────────────────────────────────
    ctx.fillStyle = '#8899BB';
    ctx.font      = '10px monospace';
    ctx.textAlign = 'center';
    ctx.fillText(
      `${enzyme || 'RE'} digestion | ${ampliconLength || '—'}bp amplicon | Predicted gel`,
      W/2, H - 10
    );

    // ── Scale bar ─────────────────────────────────────────
    ctx.fillStyle  = GEL.textColor;
    ctx.font       = '8px monospace';
    ctx.textAlign  = 'left';
    ctx.fillText('bp', 4, TOP_MARGIN + GEL_H/2);

  }, [lanes, enzyme, ampliconLength]);

  useEffect(() => { draw(); }, [draw]);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const parent = canvas.parentElement;
    const resize = () => {
      canvas.width  = parent.clientWidth - 4;
      canvas.height = 380;
      draw();
    };
    resize();
    const ro = new ResizeObserver(resize);
    ro.observe(parent);
    return () => ro.disconnect();
  }, [draw]);

  if (!lanes || lanes.length === 0) {
    return (
      <div style={{
        background: '#0D1B35', border: '1px solid #1a3060',
        borderRadius: 12, padding: 40, textAlign: 'center', color: '#8899BB',
      }}>
        <div style={{ fontSize: 32, marginBottom: 10 }}>🧫</div>
        <div style={{ fontSize: 14 }}>Select a CAPS candidate to see gel simulation</div>
        <div style={{ fontSize: 11, marginTop: 6, opacity: 0.6 }}>
          Band positions calculated from logarithmic gel migration model
        </div>
      </div>
    );
  }

  return (
    <div>
      {/* Gel Canvas */}
      <div style={{
        background: GEL.background,
        border: '1px solid #1a3060',
        borderRadius: 10, overflow: 'hidden',
        position: 'relative',
      }}>
        <canvas ref={canvasRef} style={{ display: 'block' }} />
      </div>

      {/* Legend */}
      <div style={{
        marginTop: 12, display: 'flex', gap: 16, flexWrap: 'wrap',
        background: '#0D1B35', border: '1px solid #1a3060',
        borderRadius: 8, padding: '10px 14px',
      }}>
        {[
          { color: '#8899BB', label: 'DNA Ladder' },
          { color: '#4D9FFF', label: 'Wild-type (2 bands — site cut)' },
          { color: '#F0C040', label: 'Heterozygous (3 bands)' },
          { color: '#FF4D6D', label: 'Homozygous edited (1 band — uncut)' },
        ].map(({ color, label }) => (
          <div key={label} style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
            <div style={{ width: 24, height: 4, background: color, borderRadius: 2,
                          boxShadow: `0 0 6px ${color}` }} />
            <span style={{ fontSize: 11, color: '#8899BB' }}>{label}</span>
          </div>
        ))}
      </div>

      {/* Interpretation guide */}
      <div style={{
        marginTop: 8, background: 'rgba(0,212,170,0.06)',
        border: '1px solid rgba(0,212,170,0.2)',
        borderRadius: 8, padding: '10px 14px',
      }}>
        <div style={{ fontSize: 11, color: '#00D4AA', fontWeight: 700, marginBottom: 6 }}>
          💡 How to interpret your CAPS gel
        </div>
        <div style={{ fontSize: 11, color: '#8899BB', lineHeight: 1.7 }}>
          <strong style={{ color: '#4D9FFF' }}>Wild-type:</strong> Two bands — restriction enzyme cuts the amplicon. Both alleles intact. &nbsp;|&nbsp;
          <strong style={{ color: '#FF4D6D' }}>Homozygous edit:</strong> One band — enzyme cannot cut. Both alleles disrupted. &nbsp;|&nbsp;
          <strong style={{ color: '#F0C040' }}>Heterozygous:</strong> Three bands — one allele cut, one uncut.
        </div>
      </div>
    </div>
  );
}
