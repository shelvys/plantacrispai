// ============================================================
// FILE: src/pages/AnalysisApp.js  (v3.0 — full rewrite)
// All issues from screenshots fixed:
// - MTU_FWP shown as WT Reference in all tabs
// - 1219_FWP highlighted as confirmed homozygous star
// - Zygosity colors/labels correct per class
// - Dashboard shows insertion/deletion explicitly
// - AI Insights data-driven, not static text
// - Report: PDF download + Excel CSV + detailed per-sample
// - Chromatogram: correct sample list from results
// - Explorer: correct metrics, R² up/down, ICE, identity
// - runAnalysis passes all samples including WT control
// ============================================================

import React, { useState, useRef, useCallback } from 'react';
import CAPSDesignerPanel from '../components/CAPSDesignerPanel';
import ChromatogramViewer from '../components/ChromatogramViewer';
import AB1UploadPanel from '../components/AB1UploadPanel';
import MLQualityPanel from '../components/MLQualityPanel';
import { analyseSample, deriveCutPosition } from '../engine/analysisEngine';
import { parseAB1, ab1ToSample } from '../engine/ab1Parser';
import { DEMO_CONFIG, DEMO_SAMPLES, ORGANISM_OPTIONS, PLOIDY_OPTIONS } from '../engine/demoData';

// ── Colour system ──────────────────────────────────────────
const C = {
  bg:     '#0a0f1a',
  card:   '#111827',
  card2:  '#1a2235',
  border: '#1e3a5f',
  cyan:   '#00d4aa',
  blue:   '#4d9fff',
  purple: '#b060ff',
  gold:   '#f0c040',
  red:    '#ff4d6d',
  green:  '#00d4aa',
  text:   '#e8f0fe',
  muted:  '#8899bb',
  dim:    '#3a4a5a',
};

// Zygosity colours and labels (single source of truth)
const ZYG = {
  0: { color: C.muted,  label: 'Wild-type Segregant',      badge: 'WT',     icon: '○' },
  1: { color: C.gold,   label: 'Low / Mosaic Edit',         badge: 'Low',    icon: '◑' },
  2: { color: C.blue,   label: 'Heterozygous / Chimeric',   badge: 'Het',    icon: '◑' },
  3: { color: C.cyan,   label: 'Homozygous Edited',         badge: 'Homo ★', icon: '●' },
};

// ── Reusable components ────────────────────────────────────

function Chip({ color = C.muted, children, small }) {
  return (
    <span style={{
      background: `${color}22`, border: `1px solid ${color}55`,
      color, padding: small ? '1px 7px' : '2px 10px',
      borderRadius: 20, fontSize: small ? 10 : 11,
      fontWeight: 600, display: 'inline-block', whiteSpace: 'nowrap'
    }}>{children}</span>
  );
}

function ScoreBar({ score, width = 60 }) {
  const color = score >= 80 ? C.cyan : score >= 60 ? C.gold : C.red;
  return (
    <div style={{ display: 'flex', alignItems: 'center', gap: 6 }}>
      <div style={{ width, height: 6, background: 'rgba(255,255,255,0.07)', borderRadius: 3, overflow: 'hidden' }}>
        <div style={{ width: `${score}%`, height: '100%', background: color, borderRadius: 3 }} />
      </div>
      <span style={{ color, fontSize: 11, fontWeight: 700, minWidth: 28 }}>{score}%</span>
    </div>
  );
}

function ProgressBar({ value, color = C.cyan }) {
  return (
    <div style={{ width: '100%', height: 6, background: 'rgba(255,255,255,0.05)', borderRadius: 3, overflow: 'hidden' }}>
      <div style={{ width: `${Math.min(100, value)}%`, height: '100%', background: color, borderRadius: 3, transition: 'width 0.6s ease' }} />
    </div>
  );
}

function StatCard({ label, value, color, sub, onClick }) {
  const [hov, setHov] = React.useState(false);
  return (
    <div
      onClick={onClick}
      onMouseEnter={() => setHov(true)}
      onMouseLeave={() => setHov(false)}
      style={{
        background: hov && onClick ? `${color || C.cyan}12` : C.card,
        border: `1px solid ${hov && onClick ? (color || C.cyan) + '60' : C.border}`,
        borderRadius: 14, padding: '20px 22px', textAlign: 'center',
        cursor: onClick ? 'pointer' : 'default',
        transition: 'all 0.18s',
        transform: hov && onClick ? 'translateY(-2px)' : 'none',
        boxShadow: hov && onClick ? `0 6px 20px ${(color || C.cyan)}22` : 'none',
      }}
    >
      <div style={{ fontSize: 32, fontWeight: 800, color: color || C.text, lineHeight: 1.1, marginBottom: 4 }}>{value}</div>
      <div style={{ fontSize: 10, color: C.muted, letterSpacing: 1, textTransform: 'uppercase' }}>{label}</div>
      {sub && <div style={{ fontSize: 10, color: color || C.muted, marginTop: 3 }}>{sub}</div>}
      {onClick && <div style={{ fontSize: 9, color: (color||C.cyan)+'88', marginTop: 4 }}>Click to explore →</div>}
    </div>
  );
}

function SectionTitle({ children }) {
  return <div style={{ fontSize: 13, fontWeight: 700, color: C.cyan, marginBottom: 16 }}>{children}</div>;
}

function ZygBadge({ cls }) {
  const z = ZYG[cls] || ZYG[0];
  return <Chip color={z.color}>{z.icon} {z.badge}</Chip>;
}

function EditChip({ editType, indelSize }) {
  if (!editType || editType === 'Wild-type') return <Chip color={C.dim}>Wild-type</Chip>;
  const color = editType === 'Deletion' ? C.red : editType === 'Insertion' ? C.cyan : C.gold;
  return <Chip color={color}>{editType}{indelSize > 0 ? ` (${indelSize} bp)` : ''}</Chip>;
}

// ── Chromatogram renderer → base64 PNG ─────────────────────
const BASE_COLORS_PDF = { A: '#00d4aa', T: '#ff4d6d', G: '#4d9fff', C: '#f0c040', N: '#8899bb' };

// ── Build allele frequency data from results ──────────────────────────────────
// Derives a synthetic allele table from the engine output (indel sizes + types)
// For real ICE data this would come from deconvolution; here we build it from
// the per-sample classifications which is accurate for Sanger genotyping.
function buildAlleleFreqData(results, cfg) {
  if (!results || results.length === 0) return { alleles: [], ref: '' };
  const ref    = (cfg?.ref || '').toUpperCase();   // full reference
  const cutPos = cfg?.cutPos || Math.floor(ref.length / 2);

  // Group results by edit type + size
  const groups = {};
  results.forEach(r => {
    const key = r.indelSize > 0
      ? `${r.editType === 'Deletion' ? 'del' : 'ins'}_${r.indelSize}`
      : 'wt';
    if (!groups[key]) groups[key] = { type: r.editType || 'WT', size: r.indelSize || 0, count: 0, samples: [] };
    groups[key].count++;
    groups[key].samples.push(r.id || r.name || 'unknown');
  });

  const total = results.length;
  const alleles = Object.entries(groups)
    .map(([key, g]) => ({
      ...g, freq: (g.count / total * 100).toFixed(1),
      isWT: key === 'wt',
    }))
    .sort((a, b) => b.count - a.count);

  return { alleles, ref, cutPos, total };
}

// ── Allele Plot component ─────────────────────────────────────────────────────
function AllelePlot({ results, cfg }) {
  const { alleles, ref, cutPos, total } = buildAlleleFreqData(results, cfg);
  if (!alleles || alleles.length === 0) return null;

  const BC = { A:'#00d4aa', T:'#ff6b6b', G:'#4d9fff', C:'#f0c040' };
  const WIN = 20; // bases each side of cut
  const refLeft  = (ref || '').slice(Math.max(0, cutPos - WIN), cutPos);
  const refRight = (ref || '').slice(cutPos, cutPos + WIN);
  const refSlice = refLeft + refRight;
  // const leftLen = refLeft.length; // unused

  // Build allele sequence display
  function buildSeq(a) {
    if (a.isWT) return { left: refLeft, variant: '', right: refRight, isDel: false, isIns: false };
    if (a.type === 'Deletion') {
      const delSize = Math.min(a.size, refRight.length);
      return { left: refLeft, variant: '-'.repeat(delSize), right: refRight.slice(delSize), isDel: true, isIns: false };
    }
    // Insertion — show inserted bases as ▲ blocks
    return { left: refLeft, variant: '▲'.repeat(Math.min(a.size, 6)), right: refRight, isDel: false, isIns: true };
  }

  const maxFreq = Math.max(...alleles.map(a => parseFloat(a.freq)));

  return (
    <div style={{ background:'#0d1b35', border:'1px solid #1a3060', borderRadius:14, padding:'22px 24px', marginTop:20 }}>

      {/* Header */}
      <div style={{ display:'flex', justifyContent:'space-between', alignItems:'center', marginBottom:20 }}>
        <div>
          <div style={{ fontSize:15, fontWeight:700, color:'#00D4AA' }}>🧬 Allele Distribution — Around Cut Site</div>
          <div style={{ fontSize:11, color:'#8899bb', marginTop:3 }}>
            Visualisation of identified alleles relative to the predicted Cas9 cleavage position
          </div>
        </div>
        <div style={{ display:'flex', gap:10, fontSize:11 }}>
          <span style={{ background:'rgba(0,212,170,0.1)', border:'1px solid rgba(0,212,170,0.3)', borderRadius:6, padding:'4px 10px', color:'#00d4aa' }}>{total} samples</span>
          <span style={{ background:'rgba(77,159,255,0.1)', border:'1px solid rgba(77,159,255,0.3)', borderRadius:6, padding:'4px 10px', color:'#4d9fff' }}>{alleles.length} allele classes</span>
        </div>
      </div>

      {/* Main table — full width, 3 columns: label | sequence | frequency bar */}
      <div style={{ overflowX:'auto' }}>
        <table style={{ width:'100%', borderCollapse:'collapse', fontFamily:'monospace' }}>

          {/* Column headers */}
          <thead>
            <tr>
              <th style={{ width:120, padding:'0 12px 8px 0', textAlign:'right', fontSize:10, color:'#6688aa', fontWeight:600, letterSpacing:1 }}>ALLELE</th>
              <th style={{ padding:'0 16px 8px', textAlign:'left', fontSize:10, color:'#6688aa', fontWeight:600, letterSpacing:1 }}>
                SEQUENCE AROUND CUT SITE&nbsp;&nbsp;
                <span style={{ color:'#4d9fff', fontWeight:400 }}>← {WIN}bp</span>
                <span style={{ color:'#ff4444', margin:'0 6px', borderBottom:'1px dashed #ff4444', padding:'0 2px' }}>✂ cut</span>
                <span style={{ color:'#4d9fff', fontWeight:400 }}>{WIN}bp →</span>
              </th>
              <th style={{ width:220, padding:'0 0 8px 16px', textAlign:'left', fontSize:10, color:'#6688aa', fontWeight:600, letterSpacing:1 }}>FREQUENCY</th>
              <th style={{ width:80, padding:'0 0 8px 8px', textAlign:'right', fontSize:10, color:'#6688aa', fontWeight:600, letterSpacing:1 }}>COUNT</th>
            </tr>
          </thead>

          {/* Reference row */}
          <tbody>
            <tr style={{ borderBottom:'2px solid rgba(77,159,255,0.2)' }}>
              <td style={{ padding:'6px 12px 10px 0', textAlign:'right', fontSize:11, color:'#4d9fff', fontWeight:700, whiteSpace:'nowrap' }}>
                Reference
              </td>
              <td style={{ padding:'6px 16px 10px' }}>
                <div style={{ display:'flex', gap:0, alignItems:'center' }}>
                  {/* Left side of cut */}
                  {refLeft.split('').map((b, i) => (
                    <span key={i} style={{ fontSize:13, fontWeight:700, color:BC[b]||'#fff', padding:'2px 3px', letterSpacing:'0.5px' }}>{b}</span>
                  ))}
                  {/* Cut marker */}
                  <span style={{ display:'inline-flex', flexDirection:'column', alignItems:'center', margin:'0 4px', flexShrink:0 }}>
                    <span style={{ fontSize:11, color:'#ff4444' }}>✂</span>
                    <span style={{ width:1, height:16, background:'#ff444488', display:'block', borderLeft:'2px dashed #ff4444' }} />
                  </span>
                  {/* Right side */}
                  {refRight.split('').map((b, i) => (
                    <span key={i} style={{ fontSize:13, fontWeight:700, color:BC[b]||'#fff', padding:'2px 3px', letterSpacing:'0.5px' }}>{b}</span>
                  ))}
                </div>
                {/* sgRNA underline */}
                {cfg?.grnaSeq && (
                  <div style={{ fontSize:9, color:'#4d9fff', marginTop:2, paddingLeft:4 }}>
                    {'·'.repeat(Math.min(20, refSlice.length))} ← sgRNA region
                  </div>
                )}
              </td>
              <td style={{ padding:'6px 0 10px 16px' }}>
                <div style={{ height:8, background:'rgba(77,159,255,0.15)', borderRadius:4, border:'1px solid rgba(77,159,255,0.2)' }}>
                  <div style={{ width:'100%', height:'100%', background:'rgba(77,159,255,0.4)', borderRadius:4 }} />
                </div>
              </td>
              <td style={{ padding:'6px 0 10px 8px', fontSize:11, color:'#4d9fff', textAlign:'right' }}>—</td>
            </tr>

            {/* Allele rows */}
            {alleles.slice(0, 14).map((a, ai) => {
              const seq = buildSeq(a);
              const barColor = a.isWT ? '#4d9fff' : a.type === 'Deletion' ? '#ff6b6b' : '#00d4aa';
              const barWidth = `${(parseFloat(a.freq) / maxFreq) * 100}%`;
              const rowBg = ai % 2 === 0 ? 'rgba(255,255,255,0.015)' : 'transparent';
              return (
                <tr key={ai} style={{ background: rowBg, borderBottom:'1px solid rgba(26,48,96,0.5)' }}>
                  <td style={{ padding:'5px 12px 5px 0', textAlign:'right', whiteSpace:'nowrap' }}>
                    <span style={{
                      display:'inline-block', fontSize:11, fontWeight:700,
                      color: barColor,
                      background: `${barColor}18`,
                      border: `1px solid ${barColor}44`,
                      borderRadius:5, padding:'1px 7px',
                    }}>
                      {a.isWT ? 'WT' : `${a.type === 'Deletion' ? '−' : '+'}${a.size} bp`}
                    </span>
                  </td>
                  <td style={{ padding:'5px 16px' }}>
                    <div style={{ display:'flex', gap:0, alignItems:'center', flexWrap:'nowrap' }}>
                      {/* Left bases */}
                      {seq.left.split('').map((b, i) => (
                        <span key={i} style={{ fontSize:12, fontWeight:400, color:BC[b]||'#aaa', padding:'1px 2px', opacity:0.85 }}>{b}</span>
                      ))}
                      {/* Variant region */}
                      <span style={{ margin:'0 2px' }}>
                        {seq.isDel ? (
                          seq.variant.split('').map((_, i) => (
                            <span key={i} style={{ fontSize:12, color:'#ff6b6b', textDecoration:'line-through', padding:'1px 2px', background:'rgba(255,107,107,0.1)' }}>·</span>
                          ))
                        ) : seq.isIns ? (
                          seq.variant.split('').map((b, i) => (
                            <span key={i} style={{ fontSize:12, color:'#00d4aa', fontWeight:900, padding:'1px 3px', background:'rgba(0,212,170,0.15)', border:'1px solid rgba(0,212,170,0.4)', borderRadius:3, margin:'0 1px' }}>{b}</span>
                          ))
                        ) : null}
                      </span>
                      {/* Cut marker */}
                      {!seq.isDel && !seq.isIns && (
                        <span style={{ width:1, height:14, borderLeft:'2px dashed #ff444466', display:'inline-block', margin:'0 3px', verticalAlign:'middle' }} />
                      )}
                      {/* Right bases */}
                      {seq.right.split('').map((b, i) => (
                        <span key={i} style={{ fontSize:12, fontWeight:400, color:BC[b]||'#aaa', padding:'1px 2px', opacity:0.85 }}>{b}</span>
                      ))}
                    </div>
                    {/* Sample names */}
                    <div style={{ fontSize:9, color:'#4d6a9a', marginTop:1 }}>
                      {a.samples.slice(0, 5).join(', ')}{a.samples.length > 5 ? ` +${a.samples.length - 5} more` : ''}
                    </div>
                  </td>
                  <td style={{ padding:'5px 0 5px 16px' }}>
                    <div style={{ display:'flex', alignItems:'center', gap:8 }}>
                      <div style={{ flex:1, height:10, background:'#0d1b35', borderRadius:5, border:'1px solid #1a3060', overflow:'hidden', minWidth:100 }}>
                        <div style={{
                          width: barWidth, height:'100%',
                          background: `linear-gradient(90deg, ${barColor}cc, ${barColor})`,
                          borderRadius:5, transition:'width 0.3s',
                        }} />
                      </div>
                      <span style={{ fontSize:12, fontWeight:700, color: barColor, minWidth:38, textAlign:'right' }}>
                        {a.freq}%
                      </span>
                    </div>
                  </td>
                  <td style={{ padding:'5px 0 5px 8px', textAlign:'right', fontSize:11, color:'#8899bb' }}>
                    {a.count}/{total}
                  </td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>

      {/* Legend + note */}
      <div style={{ marginTop:16, display:'flex', gap:20, flexWrap:'wrap', fontSize:11, color:'#6688aa', borderTop:'1px solid #1a3060', paddingTop:14, alignItems:'center' }}>
        <span style={{ fontWeight:700, color:'#8899bb' }}>Legend:</span>
        <span><span style={{ color:'#ff6b6b', fontWeight:700, textDecoration:'line-through' }}>·· </span>Deletions (shown as strikethrough)</span>
        <span><span style={{ color:'#00d4aa', fontWeight:700, border:'1px solid rgba(0,212,170,0.4)', padding:'0 3px', borderRadius:3 }}>▲ </span>Insertions (boxed)</span>
        <span><span style={{ fontWeight:700 }}>A G T C </span>Colour-coded bases</span>
        <span style={{ marginLeft:'auto', fontSize:10, color:'#4d6a9a', fontStyle:'italic' }}>
          Sequences shown ±{WIN} bp around predicted cleavage site
        </span>
      </div>
    </div>
  );
}

// ── Indel Characterisation Chart component ────────────────────────────────────
function IndelChart({ results, cfg }) {
  const canvasRef = React.useRef(null);
  const cutPos = cfg?.cutPos || 137;
  const refLen = (cfg?.ref || '').length || 300;

  React.useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas || !results || results.length === 0) return;
    const ctx = canvas.getContext('2d');
    const W = canvas.width, H = canvas.height;
    const PAD = { top:30, right:20, bottom:40, left:60 };
    const plotW = W - PAD.left - PAD.right;
    const plotH = H - PAD.top - PAD.bottom;

    ctx.clearRect(0, 0, W, H);
    ctx.fillStyle = '#0a1220';
    ctx.fillRect(0, 0, W, H);

    // Build position frequency arrays (per base position, grouped by type)
    const buckets = Array.from({ length: refLen }, () => ({ ins: 0, del: 0, sub: 0 }));
    results.forEach(r => {
      if (!r.indelSize || r.indelSize === 0) return;
      const pos = r.sampleCutPosition ?? cutPos;
      const p = Math.min(Math.max(0, pos), refLen - 1);
      if (r.editType === 'Deletion')  buckets[p].del += 1 / results.length * 100;
      else if (r.editType === 'Insertion') buckets[p].ins += 1 / results.length * 100;
      else buckets[p].sub += 1 / results.length * 100;
    });

    const maxY = Math.max(...buckets.map(b => b.ins + b.del + b.sub), 1);

    // Grid lines
    ctx.strokeStyle = 'rgba(100,130,180,0.12)';
    ctx.lineWidth = 1;
    for (let i = 0; i <= 4; i++) {
      const y = PAD.top + plotH - (i / 4) * plotH;
      ctx.beginPath(); ctx.moveTo(PAD.left, y); ctx.lineTo(PAD.left + plotW, y); ctx.stroke();
      ctx.fillStyle = '#4d6a9a'; ctx.font = '9px monospace'; ctx.textAlign = 'right';
      ctx.fillText((maxY * i / 4).toFixed(1) + '%', PAD.left - 4, y + 3);
    }

    // Quantification window (±20bp around cut)
    const qx1 = PAD.left + Math.max(0, cutPos - 20) / refLen * plotW;
    const qx2 = PAD.left + Math.min(refLen, cutPos + 20) / refLen * plotW;
    ctx.fillStyle = 'rgba(100,150,255,0.06)';
    ctx.fillRect(qx1, PAD.top, qx2 - qx1, plotH);
    ctx.strokeStyle = 'rgba(100,150,255,0.25)'; ctx.lineWidth = 1;
    ctx.strokeRect(qx1, PAD.top, qx2 - qx1, plotH);

    // Draw lines for each type
    const types = [
      { key: 'sub', color: '#4d9fff', label: 'Substitutions' },
      { key: 'del', color: '#ff6b6b', label: 'Deletions' },
      { key: 'ins', color: '#00d4aa', label: 'Insertions' },
    ];
    types.forEach(t => {
      ctx.beginPath();
      ctx.strokeStyle = t.color; ctx.lineWidth = 2; ctx.lineJoin = 'round';
      let first = true;
      buckets.forEach((b, i) => {
        const x = PAD.left + i / refLen * plotW;
        const y = PAD.top + plotH - (b[t.key] / maxY) * plotH;
        if (first) { ctx.moveTo(x, y); first = false; } else ctx.lineTo(x, y);
      });
      ctx.stroke();
    });

    // Cut site marker
    const cx = PAD.left + cutPos / refLen * plotW;
    ctx.strokeStyle = '#ff4444'; ctx.lineWidth = 1.5;
    ctx.setLineDash([5, 4]);
    ctx.beginPath(); ctx.moveTo(cx, PAD.top); ctx.lineTo(cx, PAD.top + plotH); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = '#ff6666'; ctx.font = 'bold 9px monospace'; ctx.textAlign = 'center';
    ctx.fillText('cut', cx, PAD.top - 4);

    // X axis labels
    ctx.fillStyle = '#4d6a9a'; ctx.font = '9px monospace'; ctx.textAlign = 'center';
    [0, 50, 100, 150, 200, 250, 300].filter(v => v <= refLen).forEach(v => {
      const x = PAD.left + v / refLen * plotW;
      ctx.fillText(v, x, PAD.top + plotH + 14);
    });
    ctx.textAlign = 'center';
    ctx.fillText('Reference amplicon position (bp)', PAD.left + plotW / 2, H - 4);

    // Legend
    types.forEach((t, i) => {
      const lx = PAD.left + i * 120;
      const ly = PAD.top - 14;
      ctx.strokeStyle = t.color; ctx.lineWidth = 2;
      ctx.beginPath(); ctx.moveTo(lx, ly); ctx.lineTo(lx + 22, ly); ctx.stroke();
      ctx.fillStyle = '#c8d6f0'; ctx.font = '9px monospace'; ctx.textAlign = 'left';
      ctx.fillText(t.label, lx + 26, ly + 3);
    });

    // Y axis label
    ctx.save();
    ctx.translate(12, PAD.top + plotH / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillStyle = '#4d6a9a'; ctx.font = '9px monospace'; ctx.textAlign = 'center';
    ctx.fillText('% Sequences', 0, 0);
    ctx.restore();
  }, [results, cutPos, refLen]);

  return (
    <div style={{ background:'#0d1b35', border:'1px solid #1a3060', borderRadius:12, padding:'18px 20px', marginTop:16 }}>
      <div style={{ fontSize:13, fontWeight:700, color:'#00D4AA', marginBottom:12 }}>
        📊 Indel Characterisation — mutation position distribution
      </div>
      <canvas ref={canvasRef} width={680} height={220}
        style={{ width:'100%', height:220, borderRadius:8, display:'block' }} />
      <div style={{ fontSize:10, color:'#6688aa', marginTop:8, textAlign:'center' }}>
        Frequency of insertions, deletions, and substitutions across the reference amplicon.
        Shaded region = ±20bp quantification window around predicted cut site.
      </div>
    </div>
  );
}

function buildSyntheticTracesPDF(sequence, quality, step = 10) {
  const seq  = (sequence || '').replace(/[^ATGCN]/gi, 'N').toUpperCase();
  const w    = seq.length * step;
  const ch   = { A: new Float32Array(w), T: new Float32Array(w), G: new Float32Array(w), C: new Float32Array(w) };
  for (let i = 0; i < seq.length; i++) {
    const base = seq[i];
    const q    = quality ? (quality[i] || 30) : 30;
    const amp  = 180 + q * 1.5;
    const sigma = step * 0.52;
    for (let dx = Math.round(-step * 2); dx <= Math.round(step * 2); dx++) {
      const px = i * step + dx;
      if (px < 0 || px >= w) continue;
      const g = Math.exp(-(dx * dx) / (2 * sigma * sigma));
      Object.keys(ch).forEach(b => {
        ch[b][px] += b === base ? g * amp : g * 10 + Math.random() * 4;
      });
    }
  }
  return { channels: ch, width: w };
}

function renderChromatogramPNG(result, cfg, cutPosition, grnaList) {
  const W = 860, H = 190;
  const canvas = document.createElement('canvas');
  canvas.width = W; canvas.height = H;
  const ctx = canvas.getContext('2d');
  const BOTTOM = H - 30;

  ctx.fillStyle = '#0a1628';
  ctx.fillRect(0, 0, W, H);

  const seq  = result.sequence || '';
  const step = 10;
  const traces = (result.traces && result.traces.channels)
    ? result.traces
    : buildSyntheticTracesPDF(seq, result.quality, step);

  // Window: focus ±50bp around cut site
  const focusS = Math.max(0, cutPosition - 35);
  const focusE = Math.min(seq.length, cutPosition + 55);
  const winLen = focusE - focusS;
  const zoom   = winLen > 0 ? W / (winLen * step) : 1;
  const offset = focusS * step * zoom;

  // gRNA highlight(s) — draw all gRNAs from grnaList with distinct colours
  const ref = cfg.ref || '';
  const grnasToRender = (grnaList && grnaList.length > 0)
    ? grnaList.filter(g => g.sequence && g.sequence.length >= 18)
    : (cfg.grna ? [{ label: 'gRNA', sequence: cfg.grna }] : []);
  const grnaColors = ['#00D4AA', '#BF5AF2', '#4D9FFF', '#F0C040'];

  grnasToRender.forEach((g, gi) => {
    if (!ref) return;
    const gseq = g.sequence.toUpperCase();
    const idx  = ref.toUpperCase().indexOf(gseq);
    if (idx < 0) return;
    const col = grnaColors[gi % grnaColors.length];
    const x1  = Math.max(0, idx * step * zoom - offset);
    const x2  = Math.min(W, (idx + gseq.length) * step * zoom - offset);
    ctx.fillStyle = `${col}18`;
    ctx.fillRect(x1, 0, x2 - x1, BOTTOM);
    ctx.fillStyle = col + 'aa';
    ctx.font = '9px monospace';
    ctx.fillText(g.label || `gRNA${gi+1}`, x1 + 2, 11 + gi * 10);
  });

  // Indel highlight
  if (result.indelSize > 0) {
    const ix = cutPosition * step * zoom - offset;
    const iw = Math.max(8, result.indelSize * step * zoom * 0.6);
    ctx.fillStyle = result.editType === 'Deletion' ? 'rgba(255,77,109,0.22)' : 'rgba(0,212,170,0.22)';
    ctx.fillRect(ix - iw / 2, 0, iw, BOTTOM);
    ctx.fillStyle = result.editType === 'Deletion' ? '#ff4d6d' : '#00d4aa';
    ctx.font = 'bold 9px monospace';
    ctx.textAlign = 'center';
    ctx.fillText(result.editType === 'Deletion' ? `−${result.indelSize}bp` : `+${result.indelSize}bp`, ix, BOTTOM + 14);
    ctx.textAlign = 'left';
  }

  // Cut site
  const cx = cutPosition * step * zoom - offset;
  if (cx >= 0 && cx <= W) {
    ctx.save();
    ctx.setLineDash([4, 3]);
    ctx.strokeStyle = '#ff4d6d';
    ctx.lineWidth   = 1.5;
    ctx.beginPath(); ctx.moveTo(cx, 0); ctx.lineTo(cx, BOTTOM); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle   = '#ff4d6d';
    ctx.font = '12px serif';
    ctx.fillText('✂', cx - 5, 14);
    ctx.restore();
  }

  // Find max for normalisation
  const chArr = Object.values(traces.channels);
  let maxV = 1;
  chArr.forEach(ch => {
    for (let i = focusS * step; i < Math.min(focusE * step, ch.length); i++) {
      if (ch[i] > maxV) maxV = ch[i];
    }
  });

  // Trace channels
  Object.entries(traces.channels).forEach(([base, ch]) => {
    ctx.strokeStyle = BASE_COLORS_PDF[base] || '#888';
    ctx.lineWidth   = 1.5;
    ctx.globalAlpha = 0.92;
    ctx.beginPath();
    let started = false;
    for (let px = 0; px < ch.length; px++) {
      const sx = px * zoom - offset;
      if (sx < -1 || sx > W + 1) continue;
      const y = BOTTOM - Math.min((ch[px] / maxV) * (BOTTOM - 8), BOTTOM - 8);
      if (started) ctx.lineTo(sx, y); else { ctx.moveTo(sx, y); started = true; }
    }
    ctx.stroke();
    ctx.globalAlpha = 1;
  });

  // Base letters
  ctx.textAlign = 'center';
  for (let i = focusS; i < focusE && i < seq.length; i++) {
    const b  = (seq[i] || 'N').toUpperCase();
    const bx = i * step * zoom - offset + (step * zoom) / 2;
    if (bx < 0 || bx > W) continue;
    ctx.fillStyle = BASE_COLORS_PDF[b] || '#8899bb';
    ctx.font = `${Math.max(7, Math.min(11, Math.round(step * zoom * 0.8)))}px monospace`;
    ctx.fillText(b, bx, H - 14);
  }
  ctx.textAlign = 'left';

  // Baseline
  ctx.strokeStyle = '#1a3060';
  ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(0, BOTTOM); ctx.lineTo(W, BOTTOM); ctx.stroke();

  // Cut-site label at bottom
  if (cx >= 0 && cx <= W) {
    ctx.fillStyle = '#ff4d6d99';
    ctx.font = '8px monospace';
    ctx.textAlign = 'center';
    ctx.fillText(`cut site (+${cutPosition}bp)`, cx, BOTTOM + 22);
    ctx.textAlign = 'left';
  }

  // Legend bottom-right
  ctx.font = '8px monospace';
  const legend = [['A','#00d4aa'],['T','#ff4d6d'],['G','#4d9fff'],['C','#f0c040']];
  legend.forEach(([b, col], li) => {
    ctx.fillStyle = col;
    ctx.fillText(`■ ${b}`, W - 130 + li * 32, H - 14);
  });

  return canvas.toDataURL('image/png');
}

// ── PDF generator with chromatogram images ──────────────────
function generatePDF(results, cfg, stats, wtSampleResult, grnaList) {
  const date  = new Date().toLocaleDateString('en-US', { year:'numeric', month:'long', day:'numeric' });
  const ZYG_L = ['Wild-type Segregant','Low / Mosaic Edit','Heterozygous / Chimeric','Homozygous Edited'];
  const ZYG_C = ['#6c757d','#fd7e14','#6f42c1','#0d6efd'];

  // Pre-compute cut position using primary gRNA (first valid one)
  const _primaryGrna = grnaList?.find(g => g.sequence?.length >= 18) || { sequence: cfg.grna || '' };
  const cutInfo = deriveCutPosition(cfg.ref || '', _primaryGrna.sequence || '');
  const cutPos  = cutInfo.cutPosition;
  const grnaStrandInfo = cutInfo.grnaStrand; // eslint-disable-line no-unused-vars

  // Render chromatogram PNG for each sample
  const sampleSections = results.map(r => {
    let imgTag = '';
    try {
      const png = renderChromatogramPNG(r, cfg, cutPos, grnaList);
      imgTag = `<img src="${png}" style="width:100%;border-radius:6px;margin:8px 0;display:block" alt="Chromatogram ${r.id}"/>`;
    } catch(e) { imgTag = '<p style="color:#aaa;font-size:11px">Chromatogram not available</p>'; }

    const zClass = r.zygosityClass || 0;
    const zColor = ZYG_C[zClass] || '#6c757d';
    const zLabel = ZYG_L[zClass] || '—';
    const isHomo = zClass === 3;

    return `
    <div class="sample" style="border-left:4px solid ${zColor};${isHomo ? 'background:#f0f7ff;' : ''}">
      <div class="sample-header">
        <div>
          <span class="sample-id" style="color:${zColor}">${isHomo ? '★ ' : ''}${r.id}</span>
          <span class="badge" style="background:${zColor}22;color:${zColor};border:1px solid ${zColor}44">${zLabel}</span>
          ${r.flag ? `<span class="badge" style="background:#fff3cd;color:#856404;border:1px solid #ffc107">⚠ ${r.flag}</span>` : ''}
        </div>
        <div class="metrics-row">
          <span class="metric"><b>Edit:</b> ${r.editType}${r.indelSize > 0 ? ` (${r.indelSize} bp)` : ''}</span>
          <span class="metric"><b>ICE:</b> ${r.iceScore || 0}%</span>
          <span class="metric"><b>R²:</b> ${(r.r2 || 0).toFixed(3)}</span>
          <span class="metric"><b>Identity:</b> ${r.sequenceIdentity ? (r.sequenceIdentity * 100).toFixed(1) + '%' : '—'}</span>
          <span class="metric"><b>ML:</b> ${r.mlScore || 0}%</span>
          <span class="metric"><b>Frameshift:</b> ${r.isFrameshift ? '<span style="color:#dc3545">Yes</span>' : 'No'}</span>
        </div>
      </div>
      <div class="chroma-label">
        Chromatogram (${r.traces?.channels ? 'real .ab1 trace' : 'synthetic from sequence'}) 
        · Cut site at +${cutPos}bp ${grnaList.length > 1 ? '· ' + grnaList.length + ' gRNAs' : (grnaList[0]?.sequence ? '· gRNA: ' + grnaList[0].sequence : '')}
        ${r.indelSize > 0 ? '· <span style="color:' + (r.editType==='Deletion'?'#dc3545':'#198754') + '">Indel highlighted</span>' : ''}
      </div>
      ${imgTag}
      <div class="protein-row" style="color:${r.isFrameshift ? '#dc3545' : r.indelSize > 0 ? '#856404' : '#198754'}">
        ${r.proteinEffect || '—'}
      </div>
    </div>`;
  }).join('');

  // Summary table
  const tableRows = results.map(r => {
    const zClass = r.zygosityClass || 0;
    return `<tr>
      <td style="color:${ZYG_C[zClass]};font-weight:600">${zClass===3?'★ ':''}${r.id}</td>
      <td>${r.generation}</td>
      <td>${r.editType}${r.indelSize>0?` (${r.indelSize} bp)`:''}</td>
      <td style="color:${ZYG_C[zClass]}">${ZYG_L[zClass]}</td>
      <td>${r.iceScore||0}%</td>
      <td>${(r.r2||0).toFixed(3)}</td>
      <td>${r.sequenceIdentity?(r.sequenceIdentity*100).toFixed(1)+'%':'—'}</td>
      <td>${r.mlScore||0}%</td>
      <td style="color:${r.isFrameshift?'#dc3545':'#198754'}">${r.isFrameshift?'Yes':'No'}</td>
      <td style="font-size:11px;color:${r.flag?'#856404':'#6c757d'}">${r.flag||'OK'}</td>
    </tr>`;
  }).join('');

  const homo = results.filter(r=>r.zygosityClass===3);
  const het  = results.filter(r=>r.zygosityClass===2);
  const wt   = results.filter(r=>r.zygosityClass===0);
  const low  = results.filter(r=>r.zygosityClass===1); // eslint-disable-line no-unused-vars
  const flagged = results.filter(r=>r.flag);

  const html = `<!DOCTYPE html>
<html><head><meta charset="UTF-8">
<title>PlantaCRISPAI Report — ${cfg.geneName || 'Analysis'}</title>
<style>
* { box-sizing:border-box; margin:0; padding:0; }
body { font-family:'Segoe UI',Arial,sans-serif; font-size:12px; color:#212529; background:#fff; padding:28px; }
@media print { body { padding:10px; } .no-print { display:none; } }
.header { display:flex; justify-content:space-between; align-items:flex-start; margin-bottom:20px; border-bottom:2px solid #0d6efd; padding-bottom:14px; }
.logo { font-size:20px; font-weight:800; color:#0d6efd; }
.subtitle { font-size:11px; color:#6c757d; margin-top:2px; }
h2 { font-size:14px; color:#212529; margin:22px 0 10px; padding-bottom:5px; border-bottom:1px solid #dee2e6; }
h3 { font-size:12px; color:#343a40; margin:14px 0 7px; }
.meta-grid { display:grid; grid-template-columns:repeat(4,1fr); gap:10px; margin-bottom:20px; }
.meta-card { background:#f8f9fa; border:1px solid #dee2e6; border-radius:7px; padding:10px 12px; }
.meta-card .val { font-size:22px; font-weight:700; color:#0d6efd; line-height:1.1; }
.meta-card .lbl { font-size:9px; color:#6c757d; text-transform:uppercase; letter-spacing:0.5px; margin-top:2px; }
table { width:100%; border-collapse:collapse; font-size:11px; margin-bottom:16px; }
thead tr { background:#0d6efd; color:#fff; }
thead th { padding:8px 9px; text-align:left; font-weight:600; font-size:10px; }
tbody tr:nth-child(even) { background:#f8f9fa; }
td { padding:7px 9px; border-bottom:1px solid #dee2e6; vertical-align:middle; }
.sample { border-radius:8px; padding:14px; margin-bottom:20px; border:1px solid #dee2e6; page-break-inside:avoid; }
.sample-header { display:flex; justify-content:space-between; align-items:flex-start; margin-bottom:6px; flex-wrap:wrap; gap:8px; }
.sample-id { font-size:15px; font-weight:700; margin-right:8px; }
.badge { font-size:10px; padding:2px 8px; border-radius:20px; margin-right:6px; font-weight:600; }
.metrics-row { display:flex; flex-wrap:wrap; gap:10px; margin-top:4px; }
.metric { font-size:11px; color:#495057; }
.chroma-label { font-size:9px; color:#6c757d; margin:4px 0; font-style:italic; }
.protein-row { font-size:11px; margin-top:6px; font-weight:600; }
.methods { background:#f8f9fa; border-left:3px solid #0d6efd; padding:12px 16px; margin-top:10px; font-size:11px; line-height:1.75; color:#495057; }
.recs { margin-top:8px; }
.rec { padding:9px 12px; border-radius:6px; margin-bottom:7px; font-size:11px; line-height:1.55; }
.rec.high  { background:#fff3cd; border-left:3px solid #ffc107; }
.rec.prio  { background:#d1e7dd; border-left:3px solid #198754; }
.rec.info  { background:#cfe2ff; border-left:3px solid #0d6efd; }
.footer { margin-top:24px; font-size:9px; color:#adb5bd; border-top:1px solid #dee2e6; padding-top:10px; }
.print-btn { position:fixed; top:16px; right:16px; background:#0d6efd; color:#fff; border:none; padding:10px 20px; border-radius:8px; cursor:pointer; font-size:13px; font-weight:600; z-index:999; }
</style></head>
<body>
<button class="print-btn no-print" onclick="window.print()">🖨️ Save as PDF</button>

<div class="header">
  <div>
    <div class="logo">🧬 PlantaCRISPAI Analysis Report</div>
    <div class="subtitle">CRISPR Genome Editing Platform · v2.0 · Generated: ${date}</div>
  </div>
  <div style="text-align:right;font-size:11px;color:#6c757d;line-height:1.8">
    <b>Gene:</b> ${cfg.geneName||'—'}<br>
    <b>gRNA${grnaList.length > 1 ? 's' : ''}:</b> <code>${grnaList.filter(g=>g.sequence).map(g=>`${g.label}: ${g.sequence}`).join(' | ') || cfg.grna || '—'}</code><br>
    <b>Organism:</b> ${cfg.organism||'—'} · <b>Generation:</b> ${cfg.generation||'—'}<br>
    <b>Cut site:</b> +${cutPos} bp from amplicon start<br>
    ${wtSampleResult ? `<b>WT Reference:</b> ${wtSampleResult.id}<br>` : ''}
  </div>
</div>

<div class="meta-grid">
  <div class="meta-card"><div class="val">${stats.eff.toFixed(1)}%</div><div class="lbl">Editing Efficiency</div></div>
  <div class="meta-card"><div class="val">${stats.edited}/${stats.total}</div><div class="lbl">Edited / Total</div></div>
  <div class="meta-card"><div class="val">${stats.homo}</div><div class="lbl">Homozygous ★</div></div>
  <div class="meta-card"><div class="val">${stats.het}</div><div class="lbl">Heterozygous</div></div>
  <div class="meta-card"><div class="val">${stats.fs}</div><div class="lbl">Frameshifts</div></div>
  <div class="meta-card"><div class="val">${stats.wt}</div><div class="lbl">Wild-type Segregants</div></div>
  <div class="meta-card"><div class="val">${stats.low}</div><div class="lbl">Low / Mosaic</div></div>
  <div class="meta-card"><div class="val">${stats.flagged}</div><div class="lbl">Needs Re-sequencing</div></div>
</div>

<h2>Summary Table</h2>
<table>
  <thead><tr><th>Sample</th><th>Gen.</th><th>Edit</th><th>Zygosity</th><th>ICE</th><th>R²</th><th>Identity</th><th>ML%</th><th>Frameshift</th><th>Flag</th></tr></thead>
  <tbody>${tableRows}</tbody>
</table>

<h2>AI Recommendations</h2>
<div class="recs">
  ${homo.length>0 ? `<div class="rec prio">🏆 <b>Advance to T2 immediately:</b> ${homo.map(r=>'★ '+r.id).join(', ')} — confirmed homozygous edited. Self-pollinate to fix the edit in T2.</div>` : ''}
  ${het.filter(r=>!r.flag).length>0 ? `<div class="rec info">🌿 <b>Heterozygous — advance with self-pollination:</b> ${het.filter(r=>!r.flag).map(r=>r.id).join(', ')} — one allele edited. Screen ~8 T2 plants to recover homozygotes (~25%).</div>` : ''}
  ${flagged.length>0 ? `<div class="rec high">⚠️ <b>Re-sequence before advancing:</b> ${flagged.map(r=>r.id).join(', ')} — low alignment confidence. Clonal sequencing recommended.</div>` : ''}
  ${wt.length>0 ? `<div class="rec info">🌱 <b>Null-segregant controls:</b> ${wt.map(r=>r.id).join(', ')} — wild-type. Use for phenotyping comparisons.</div>` : ''}
  <div class="rec info">🔬 <b>Transgene screening:</b> Screen all edited plants for T-DNA (Cas9/Hygromycin PCR). Transgene-free edited plants (~37% of T1) are your cleanest lines.</div>
</div>

<h2>Per-Sample Chromatograms & Analysis</h2>
${sampleSections}

<h2>Allele Distribution — Around Cut Site</h2>
<div style="background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:20px;margin-bottom:24px;">
  <div style="font-size:13px;font-weight:700;color:#1a1a2e;margin-bottom:4px;">
    ${cfg.geneName||'Gene'} · gRNA: ${grnaList[0]?.sequence||cfg.grna||'—'} · Cut site: +${cutPos} bp
  </div>
  <div style="font-size:11px;color:#666;margin-bottom:16px;">
    Visualisation of identified alleles relative to the predicted Cas9 cleavage position
  </div>
  <!-- Reference row -->
  <table style="width:100%;border-collapse:collapse;font-family:monospace;font-size:12px;margin-bottom:12px;">
    <tr style="border-bottom:2px solid #dee2e6;">
      <td style="width:80px;text-align:right;padding:6px 12px 10px 0;color:#3182ce;font-weight:700;font-size:11px;">Reference</td>
      <td style="padding:6px 16px 10px;">
        ${(()=>{
          const BC = {A:'#00a87a',T:'#e53e3e',G:'#d69e2e',C:'#3182ce'};
          const refL = (cfg.ref||'').slice(Math.max(0,cutPos-20),cutPos);
          const refR = (cfg.ref||'').slice(cutPos, cutPos+20);
          const leftBases = refL.split('').map(b=>`<span style="color:${BC[b]||'#333'};font-weight:700;padding:0 2px;">${b}</span>`).join('');
          const rightBases = refR.split('').map(b=>`<span style="color:${BC[b]||'#333'};font-weight:700;padding:0 2px;">${b}</span>`).join('');
          return leftBases + '<span style="border-left:2px dashed #e53e3e;margin:0 4px;padding:0 2px;color:#e53e3e;font-size:14px;">✂</span>' + rightBases;
        })()}
      </td>
      <td style="width:160px;padding:6px 0 10px 16px;color:#3182ce;font-size:11px;">← 20bp · cut · 20bp →</td>
      <td style="width:50px;text-align:right;padding:6px 0 10px 8px;color:#666;font-size:11px;">—</td>
    </tr>
    ${(()=>{
      const BC = {A:'#00a87a',T:'#e53e3e',G:'#d69e2e',C:'#3182ce'};
      const groups = {};
      results.forEach(r => {
        const key = r.indelSize > 0 ? `${r.editType}_${r.indelSize}` : 'WT';
        if (!groups[key]) groups[key] = {type:r.editType||'WT',size:r.indelSize||0,count:0,samples:[]};
        groups[key].count++;
        groups[key].samples.push(r.id||r.name);
      });
      const total = results.length;
      const maxFreq = Math.max(...Object.values(groups).map(g=>g.count/total*100));
      return Object.entries(groups)
        .sort((a,b)=>b[1].count-a[1].count)
        .slice(0,12)
        .map(([k,g],ri)=>{
          const freq = (g.count/total*100).toFixed(1);
          const isWT = k==='WT';
          const isDel = g.type==='Deletion';
          const barCol = isWT?'#3182ce':isDel?'#e53e3e':'#00a87a';
          const label = isWT ? 'WT' : `${isDel?'−':'+'}${g.size} bp`;
          const refL = (cfg.ref||'').slice(Math.max(0,cutPos-20),cutPos);
          const refR = (cfg.ref||'').slice(cutPos,cutPos+20);
          let seqHtml = '';
          if(isWT){
            seqHtml = [...refL].map(b=>`<span style="color:${BC[b]||'#333'};padding:0 2px;">${b}</span>`).join('')
              + '<span style="border-left:1px dashed #ccc;margin:0 4px;"></span>'
              + [...refR].map(b=>`<span style="color:${BC[b]||'#333'};padding:0 2px;">${b}</span>`).join('');
          } else if(isDel){
            const dels = Math.min(g.size,10);
            seqHtml = [...refL].map(b=>`<span style="color:${BC[b]||'#333'};padding:0 2px;">${b}</span>`).join('')
              + `<span style="color:#e53e3e;text-decoration:line-through;padding:0 2px;">${'·'.repeat(dels)}</span>`
              + [...refR.slice(dels)].map(b=>`<span style="color:${BC[b]||'#333'};padding:0 2px;">${b}</span>`).join('');
          } else {
            const ins = Math.min(g.size,6);
            seqHtml = [...refL].map(b=>`<span style="color:${BC[b]||'#333'};padding:0 2px;">${b}</span>`).join('')
              + `<span style="color:#00a87a;border:1px solid #00a87a;border-radius:2px;padding:0 3px;font-weight:700;">${'▲'.repeat(ins)}</span>`
              + [...refR].map(b=>`<span style="color:${BC[b]||'#333'};padding:0 2px;">${b}</span>`).join('');
          }
          const barW = Math.round((parseFloat(freq)/maxFreq)*140);
          const bg = ri%2===0?'#f8f9fa':'#fff';
          return `<tr style="background:${bg};border-bottom:1px solid #eee;">
            <td style="text-align:right;padding:5px 12px 5px 0;"><span style="background:${barCol}18;border:1px solid ${barCol}44;color:${barCol};border-radius:4px;padding:2px 7px;font-size:11px;font-weight:700;">${label}</span></td>
            <td style="padding:5px 16px;font-size:12px;">${seqHtml}<br><span style="font-size:9px;color:#999;">${g.samples.slice(0,5).join(', ')}${g.samples.length>5?' +more':''}</span></td>
            <td style="padding:5px 0 5px 16px;">
              <div style="display:flex;align-items:center;gap:6px;">
                <div style="width:${barW}px;height:8px;background:${barCol};border-radius:4px;"></div>
                <span style="font-weight:700;color:${barCol};font-size:11px;">${freq}%</span>
              </div>
            </td>
            <td style="text-align:right;padding:5px 0 5px 8px;font-size:11px;color:#666;">${g.count}/${total}</td>
          </tr>`;
        }).join('');
    })()}
  </table>
  <div style="font-size:10px;color:#666;border-top:1px solid #dee2e6;padding-top:8px;display:flex;gap:16px;">
    <span><span style="color:#e53e3e;text-decoration:line-through;">·· </span>Deletions</span>
    <span><span style="color:#00a87a;border:1px solid #00a87a;padding:0 2px;">▲</span> Insertions</span>
    <span style="color:#3182ce;font-weight:700;">— </span>Wild-type
    <span style="margin-left:auto;">Dashed line = predicted cleavage position (between bp ${cutPos-1}–${cutPos})</span>
  </div>
</div>

<h2>Indel Position Distribution — Mutation Frequency Across Amplicon</h2>
<div style="background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:20px;margin-bottom:24px;">
  <div style="font-size:11px;color:#666;margin-bottom:12px;">
    Frequency of insertions, deletions, and substitutions across the reference amplicon.
    Shaded region = ±20bp quantification window around predicted cut site (+${cutPos} bp).
  </div>
  ${(()=>{
    const refLen = (cfg.ref||'').length || 300;
    const total = results.length;
    const W = 700, H = 180;
    const PAD = {top:24,right:16,bottom:36,left:52};
    const plotW = W-PAD.left-PAD.right, plotH = H-PAD.top-PAD.bottom;
    const buckets = Array.from({length:refLen},()=>({ins:0,del:0,sub:0}));
    results.forEach(r=>{
      if(!r.indelSize||r.indelSize===0) return;
      const p = Math.min(Math.max(0,(r.sampleCutPosition??cutPos)),refLen-1);
      if(r.editType==='Deletion') buckets[p].del+=1/total*100;
      else if(r.editType==='Insertion') buckets[p].ins+=1/total*100;
      else buckets[p].sub+=1/total*100;
    });
    const maxY = Math.max(...buckets.map(b=>b.ins+b.del+b.sub),1);
    const toX = i => PAD.left + i/refLen*plotW;
    const toY = v => PAD.top + plotH - (v/maxY)*plotH;
    // Quantification window
    const qx1 = toX(Math.max(0,cutPos-20));
    const qx2 = toX(Math.min(refLen,cutPos+20));
    // Build SVG paths
    const mkPath = (key,col)=>{
      const pts = buckets.map((b,i)=>`${toX(i)},${toY(b[key])}`).join(' ');
      // const first = ... (unused)
      // const last = ... (unused)
      return `<polyline points="${pts}" fill="none" stroke="${col}" stroke-width="1.5" stroke-linejoin="round"/>`;
    };
    // Y axis labels
    const yLabels = [0,25,50,75,100].filter(v=>v<=maxY*1.1).map(v=>{
      const y = toY(v/100*maxY);
      return `<text x="${PAD.left-4}" y="${y+3}" text-anchor="end" font-size="8" fill="#888">${(v/100*maxY).toFixed(1)}%</text>
              <line x1="${PAD.left}" y1="${y}" x2="${PAD.left+plotW}" y2="${y}" stroke="#eee" stroke-width="0.5"/>`;
    }).join('');
    // X axis labels
    const xTicks = [0,50,100,150,200,250,300].filter(v=>v<=refLen);
    const xLabels = xTicks.map(v=>{
      const x = toX(v);
      return `<text x="${x}" y="${PAD.top+plotH+14}" text-anchor="middle" font-size="8" fill="#888">${v}</text>
              <line x1="${x}" y1="${PAD.top}" x2="${x}" y2="${PAD.top+plotH}" stroke="#f0f0f0" stroke-width="0.5"/>`;
    }).join('');
    const cutX = toX(cutPos);
    return `<svg width="${W}" height="${H}" xmlns="http://www.w3.org/2000/svg" style="display:block;max-width:100%;">
      <!-- Quant window -->
      <rect x="${qx1}" y="${PAD.top}" width="${qx2-qx1}" height="${plotH}" fill="#3182ce" fill-opacity="0.06"/>
      <rect x="${qx1}" y="${PAD.top}" width="${qx2-qx1}" height="${plotH}" fill="none" stroke="#3182ce" stroke-width="0.5" stroke-opacity="0.3"/>
      <!-- Grid -->
      ${yLabels}${xLabels}
      <!-- Cut site -->
      <line x1="${cutX}" y1="${PAD.top}" x2="${cutX}" y2="${PAD.top+plotH}" stroke="#e53e3e" stroke-width="1.2" stroke-dasharray="4,3"/>
      <text x="${cutX}" y="${PAD.top-4}" text-anchor="middle" font-size="8" fill="#e53e3e" font-weight="700">cut</text>
      <!-- Data lines -->
      ${mkPath('sub','#3182ce')}
      ${mkPath('del','#e53e3e')}
      ${mkPath('ins','#00a87a')}
      <!-- Axes -->
      <line x1="${PAD.left}" y1="${PAD.top}" x2="${PAD.left}" y2="${PAD.top+plotH}" stroke="#ccc" stroke-width="1"/>
      <line x1="${PAD.left}" y1="${PAD.top+plotH}" x2="${PAD.left+plotW}" y2="${PAD.top+plotH}" stroke="#ccc" stroke-width="1"/>
      <!-- X axis title -->
      <text x="${PAD.left+plotW/2}" y="${H-4}" text-anchor="middle" font-size="9" fill="#666">Reference amplicon position (bp)</text>
      <!-- Y axis title -->
      <text transform="translate(10,${PAD.top+plotH/2}) rotate(-90)" text-anchor="middle" font-size="9" fill="#666">% Sequences</text>
      <!-- Legend -->
      <line x1="${PAD.left}" y1="10" x2="${PAD.left+18}" y2="10" stroke="#00a87a" stroke-width="2"/>
      <text x="${PAD.left+22}" y="13" font-size="8" fill="#333">Insertions</text>
      <line x1="${PAD.left+75}" y1="10" x2="${PAD.left+93}" y2="10" stroke="#e53e3e" stroke-width="2"/>
      <text x="${PAD.left+97}" y="13" font-size="8" fill="#333">Deletions</text>
      <line x1="${PAD.left+148}" y1="10" x2="${PAD.left+166}" y2="10" stroke="#3182ce" stroke-width="2"/>
      <text x="${PAD.left+170}" y="13" font-size="8" fill="#333">Substitutions</text>
      <rect x="${qx1}" y="5" width="${Math.max(10,qx2-qx1)}" height="9" fill="#3182ce" fill-opacity="0.15" stroke="#3182ce" stroke-width="0.5"/>
      <text x="${qx2+3}" y="13" font-size="7" fill="#3182ce">Quant. window</text>
    </svg>`;
  })()}
</div>

<h2>Indel Characterisation — Summary by Type & Size</h2>
<table style="width:100%;border-collapse:collapse;font-size:12px;margin-bottom:24px;">
  <thead>
    <tr style="background:#2d3748;color:#fff;">
      <th style="padding:8px 12px;text-align:left;">Edit Type</th>
      <th style="padding:8px 12px;text-align:center;">Indel Size (bp)</th>
      <th style="padding:8px 12px;text-align:center;">Count</th>
      <th style="padding:8px 12px;text-align:center;">Frequency (%)</th>
      <th style="padding:8px 12px;text-align:center;">Frameshift</th>
      <th style="padding:8px 12px;text-align:left;">Sample IDs</th>
    </tr>
  </thead>
  <tbody>
    ${(()=>{
      const rows = [];
      const edited = results.filter(r=>r.indelSize>0);
      const groups = {};
      edited.forEach(r=>{
        const key=`${r.editType}_${r.indelSize}`;
        if(!groups[key]) groups[key]={type:r.editType,size:r.indelSize,count:0,fs:0,samples:[]};
        groups[key].count++;
        if(r.isFrameshift) groups[key].fs++;
        groups[key].samples.push(r.id||r.name);
      });
      if(Object.keys(groups).length===0) return '<tr><td colspan="6" style="padding:12px;text-align:center;color:#666;">No indels detected — all samples are wild-type</td></tr>';
      Object.values(groups).sort((a,b)=>b.count-a.count).forEach((g,i)=>{
        const freq=(g.count/results.length*100).toFixed(1);
        const bg=i%2===0?'#f8f9fa':'#fff';
        const col=g.type==='Deletion'?'#e53e3e':'#00a87a';
        rows.push(`<tr style="background:${bg};">
          <td style="padding:7px 12px;"><span style="color:${col};font-weight:600;">${g.type}</span></td>
          <td style="padding:7px 12px;text-align:center;font-weight:700;">${g.size} bp</td>
          <td style="padding:7px 12px;text-align:center;">${g.count}</td>
          <td style="padding:7px 12px;text-align:center;">${freq}%</td>
          <td style="padding:7px 12px;text-align:center;">${g.fs>0?`<span style="color:#e53e3e;">${g.fs} (${(g.fs/g.count*100).toFixed(0)}%)</span>`:'—'}</td>
          <td style="padding:7px 12px;font-size:10px;color:#666;">${g.samples.slice(0,6).join(', ')}${g.samples.length>6?' ...':''}</td>
        </tr>`);
      });
      return rows.join('');
    })()}
    <tr style="background:#edf2f7;font-weight:700;">
      <td colspan="2" style="padding:7px 12px;">Total edited</td>
      <td style="padding:7px 12px;text-align:center;">${results.filter(r=>r.indelSize>0).length}</td>
      <td style="padding:7px 12px;text-align:center;">${(results.filter(r=>r.indelSize>0).length/results.length*100).toFixed(1)}%</td>
      <td style="padding:7px 12px;text-align:center;">${results.filter(r=>r.isFrameshift).length}</td>
      <td style="padding:7px 12px;"></td>
    </tr>
  </tbody>
</table>

<h2>Methods</h2>
<div class="methods">
Genome editing efficiency and mutation characterization were performed using PlantaCRISPAI (v2.0). Sanger sequencing reads (.ab1 format) were parsed from Applied Biosystems trace files. Sample sequences were aligned to the reference wild-type amplicon using the Smith-Waterman local alignment algorithm (match=+2, mismatch=−1, gap open=−2). Indels were identified from gap patterns in pairwise alignments. Zygosity classification employed two orthogonal signals: (1) the alignment-detected indel size from gap character counts, used as the primary homozygous indicator; and (2) the downstream non-gap mismatch rate (≥15% = double-peak Sanger trace, indicating heterozygosity). R² metrics were computed separately on upstream and downstream amplicon segments relative to the predicted Cas9 cut site (gRNA 3′ end − 3 bp), excluding gap columns from comparison. Editing efficiency was calculated as confirmed edited samples / total confirmed samples. Chromatogram images in this report are rendered from ${results.some(r=>r.traces?.channels) ? 'raw .ab1 fluorescence trace data' : 'sequence-derived synthetic traces'}.
<br><br>
<b>Reference amplicon:</b> ${cfg.ref ? cfg.ref.length + ' bp, GC ' + (((cfg.ref.match(/[GC]/g)||[]).length/cfg.ref.length)*100).toFixed(1) + '%' : '—'}<br>
<b>Guide RNA${grnaList.length > 1 ? 's' : ''}:</b> ${grnaList.filter(g=>g.sequence).map(g=>`${g.label}: ${g.sequence}`).join(' | ') || cfg.grna || '—'} · <b>Predicted cut site${grnaList.length > 1 ? 's' : ''}:</b> +${cutPos} bp${grnaList.length > 1 ? ' (primary)' : ''}
</div>

<div class="footer">
  PlantaCRISPAI v2.0 · For research use only · ${date} ·
  Gene: ${cfg.geneName||'—'} · ${grnaList.length > 1 ? grnaList.length + ' gRNAs' : 'gRNA: ' + (grnaList[0]?.sequence || cfg.grna || '—')} · ${cfg.organism||'—'} ${cfg.generation||'—'}
</div>
</body></html>`;

  const win = window.open('', '_blank');
  if (!win) { alert('Pop-up blocked — please allow pop-ups for this site to generate the PDF.'); return; }
  win.document.write(html);
  win.document.close();
}


// ── Excel CSV export (proper for Excel) ───────────────────
function downloadExcel(results, cfg, stats, wtControl, grnaList) {
  const BOM = '\uFEFF';
  const zyg = r => ZYG[r.zygosityClass]?.label || '—';
  const header = [
    'Sample ID','Generation','Edit Type','Indel Size (bp)',
    'Zygosity Class','Zygosity Label','ICE Score (%)','R² Overall',
    'R² Upstream','R² Downstream','Sequence Identity (%)','ML Score (%)',
    'Is Frameshift','Protein Effect',
    'Read Strand','gRNA Strand','Read Auto-Flipped',
    'Flag / Note'
  ];
  const rows = results.map(r => [
    r.id, r.generation,
    r.editType || '—',
    r.indelSize || 0,
    r.zygosityClass,
    zyg(r),
    r.iceScore || 0,
    (r.r2 || 0).toFixed(3),
    (r.r2up || 0).toFixed(3),
    (r.r2dn || 0).toFixed(3),
    r.sequenceIdentity ? (r.sequenceIdentity * 100).toFixed(1) : '—',
    r.mlScore || 0,
    r.isFrameshift ? 'Yes' : 'No',
    `"${r.proteinEffect || '—'}"`,
    r.readStrand || '+',
    r.grnaStrand || '+',
    r.wasFlipped ? 'Yes' : 'No',
    `"${r.flag || 'OK'}"`,
  ]);

  const csv = BOM + [header, ...rows].map(row => row.join('\t')).join('\n');
  const a = document.createElement('a');
  a.href = URL.createObjectURL(new Blob([csv], { type: 'text/tab-separated-values' }));
  a.download = `PlantaCRISPAI_${cfg.geneName || 'Analysis'}_Results.xls`;
  a.click();
}

// ── Main Component ─────────────────────────────────────────
function AnalysisApp({ onBack }) {
  const [tab, setTab]                   = useState(0);
  // Navigation helper — always scroll to top when changing tabs
  const navigateTo = useCallback((newTab, sample = null) => {
    setTab(newTab);
    if (sample) setSelectedSample(sample);
    window.scrollTo({ top: 0, behavior: 'smooth' });
  }, []);
  const [cfg, setCfg]                   = useState({
    geneName: '', ref: '', grna: '',
    generation: 'T1', ploidy: 'diploid', organism: 'rice'
  });
  // Phase 3: multi-gRNA list — each entry { label, sequence }
  // When length > 1, all gRNAs are analysed independently per sample
  const [grnaList, setGrnaList] = useState([{ label: 'gRNA 1', sequence: '' }]);
  const [samples, setSamples]           = useState([]);
  const [results, setResults]           = useState(null);
  const [wtControlResult, setWTResult]  = useState(null);  // WT reference result row
  const [selectedSample, setSelectedSample] = useState(null);
  const [busy, setBusy]                 = useState(false);
  const [progress, setProgress]         = useState(0);
  const [progressMsg, setProgressMsg]   = useState('');
  // Phase 3: trimming state — shared across Setup tab
  const [trimMode, setTrimMode]         = useState('standard'); // eslint-disable-line no-unused-vars
  const [trimOpts, setTrimOpts]         = useState({ phredThreshold: 20, window: 10 }); // eslint-disable-line no-unused-vars
  // Explorer filter — set by clicking Dashboard stat cards
  const [explorerFilter, setExplorerFilter] = useState('all'); // 'all'|'edited'|'homo'|'het'|'fs'|'wt'|'flagged'
  const fileRef                         = useRef(); // eslint-disable-line no-unused-vars

  // Auto-select first visible sample when Explorer filter changes or when navigating to Explorer
  React.useEffect(() => {
    if (tab !== 2 || !results) return;
    const filtered = results.filter(r => {
      if (explorerFilter === 'edited')  return r.zygosityClass > 0;
      if (explorerFilter === 'homo')    return r.zygosityClass === 3;
      if (explorerFilter === 'het')     return r.zygosityClass === 2;
      if (explorerFilter === 'wt')      return r.zygosityClass === 0;
      if (explorerFilter === 'fs')      return r.isFrameshift;
      if (explorerFilter === 'flagged') return !!r.flag;
      return true;
    });
    // If current selection is not in filtered list, auto-select first
    const currentInFilter = filtered.some(r => r.id === selectedSample?.id);
    if (!currentInFilter && filtered.length > 0) {
      setSelectedSample(filtered[0]);
    }
  }, [explorerFilter, tab, results, selectedSample?.id]); // eslint-disable-line react-hooks/exhaustive-deps

  // ── Demo data ──────────────────────────────────────────
  const loadDemo = () => {
    setCfg({
      geneName: DEMO_CONFIG.geneName,
      ref:      DEMO_CONFIG.ref,
      grna:     DEMO_CONFIG.grna,
      generation: DEMO_CONFIG.generation,
      ploidy:   DEMO_CONFIG.ploidy,
      organism: DEMO_CONFIG.organism,
    });
    setGrnaList([{ label: 'gRNA 1', sequence: DEMO_CONFIG.grna }]);
    setSamples(DEMO_SAMPLES.map(s => ({ ...s })));
    setResults(null);
  };

  // ── Run Analysis ────────────────────────────────────────
  const runAnalysis = async () => {
    // Require at least one valid gRNA and a reference
    const validGrnas = grnaList.filter(g => g.sequence && g.sequence.length >= 18);
    if (!validGrnas.length || !cfg.ref || !samples.length) return;
    setBusy(true);
    setProgress(0);

    const steps = [
      { p: 15, msg: 'Parsing sequences...' },
      { p: 35, msg: 'Running Smith-Waterman alignment...' },
      { p: 60, msg: `Detecting indels at ${validGrnas.length} cut site${validGrnas.length > 1 ? 's' : ''}...` },
      { p: 80, msg: 'Classifying zygosity...' },
      { p: 95, msg: 'Computing ML confidence scores...' },
      { p: 100, msg: 'Done.' },
    ];
    for (const step of steps) {
      await new Promise(r => setTimeout(r, 180));
      setProgress(step.p);
      setProgressMsg(step.msg);
    }

    // Primary gRNA for backwards-compat single-gRNA fields
    const primaryGrna = validGrnas[0];
    const cutInfo = deriveCutPosition(cfg.ref, primaryGrna.sequence);
    const cutPos  = cutInfo.cutPosition;

    const config = {
      grnaSeq:     primaryGrna.sequence,
      grnaList:    validGrnas,          // Phase 3: full multi-gRNA list
      cutPosition: cutPos,
      refAmplicon: cfg.ref,
      generation:  cfg.generation,
      ploidy:      cfg.ploidy,
      geneName:    cfg.geneName,
      organism:    cfg.organism,
    };

    // Run ALL samples (reference-less first pass to get sequences)
    // Then re-run using the best-quality WT sample as the alignment reference
    const allResults = samples.map(s => analyseSample(s, null, config));

    // Identify the WT reference sample: class 0, highest sequence identity, no flag
    const wtCandidates = allResults
      .filter(r => r.zygosityClass === 0 && !r.flag)
      .sort((a, b) => (b.sequenceIdentity || 0) - (a.sequenceIdentity || 0));
    const wtResult = wtCandidates.length > 0
      ? { ...wtCandidates[0], isWTControl: true }
      : null;

    // Find the matching sample object for the WT result (to use its traces)
    const wtSample = wtResult
      ? samples.find(s => (s.id || s.name) === (wtResult.id || wtResult.name))
      : null;

    // Re-run all samples aligned to the WT sample (better accuracy)
    const finalRes = wtSample
      ? samples.map(s => analyseSample(s, wtSample, config))
      : allResults;

    setResults(finalRes);
    setWTResult(finalRes.find(r => r.id === (wtResult?.id)) || null);
    setSelectedSample(finalRes[0] || null);
    setBusy(false);
    navigateTo(1);
  };

  // ── File upload — BINARY for .ab1, text for others ─────
  const refFileRef = useRef();

  // ── Reference file upload (FASTA or .ab1) → cfg.ref ─────
  const handleRefFile = useCallback((e) => {
    const file = Array.from(e.target?.files || e.dataTransfer?.files || [])[0];
    if (!file) return;
    const isAB1 = file.name.toLowerCase().endsWith('.ab1');
    const reader = new FileReader();
    if (isAB1) {
      reader.onload = ev => {
        const parsed = parseAB1(ev.target.result, file.name);
        if (parsed.success && parsed.trimmedSeq) {
          setCfg(c => ({ ...c, ref: (parsed.trimmedSeq || parsed.sequence || '').toUpperCase().replace(/[^ATGCN]/g, '') }));
        }
      };
      reader.readAsArrayBuffer(file);
    } else {
      reader.onload = ev => {
        const text = ev.target.result || '';
        const seq = text.split('\n')
          .filter(l => !l.startsWith('>') && !l.startsWith(';') && l.trim())
          .join('').replace(/[^ATGCN]/gi, '').toUpperCase();
        if (seq.length > 20) setCfg(c => ({ ...c, ref: seq }));
      };
      reader.readAsText(file);
    }
  }, []);

  const handleFiles = useCallback((e) => { // eslint-disable-line no-unused-vars
    Array.from(e.target.files).forEach(file => {
      const isAB1 = file.name.toLowerCase().endsWith('.ab1');
      const reader = new FileReader();

      if (isAB1) {
        reader.onload = ev => {
          const parsed = parseAB1(ev.target.result, file.name);
          if (parsed.success) {
            const s = ab1ToSample(parsed, cfg.generation, cfg.ploidy, trimMode, trimOpts);
            if (s) setSamples(p => [...p, s]);
          } else {
            setSamples(p => [...p, {
              id: file.name.replace(/\.ab1$/i, ''),
              generation: cfg.generation, ploidy: cfg.ploidy,
              sequence: '', parseError: parsed.error || 'Failed to parse .ab1 file',
            }]);
          }
        };
        reader.readAsArrayBuffer(file);
      } else {
        reader.onload = ev => {
          const seq = (ev.target.result || '')
            .split('\n').filter(l => !l.startsWith('>') && !l.startsWith(';'))
            .join('').replace(/[^ATGCN]/gi, '').toUpperCase();
          setSamples(p => [...p, {
            id: file.name.replace(/\.[^.]+$/, ''),
            generation: cfg.generation, ploidy: cfg.ploidy, sequence: seq,
          }]);
        };
        reader.readAsText(file);
      }
    });
  }, [cfg.generation, cfg.ploidy, trimMode, trimOpts]);

  const addSample  = () => setSamples(p => [...p, { id: `Sample_${p.length + 1}`, generation: cfg.generation, ploidy: cfg.ploidy, sequence: '' }]);
  const delSample  = i => setSamples(p => p.filter((_, j) => j !== i));
  const updSample  = (i, k, v) => setSamples(p => p.map((s, j) => j === i ? { ...s, [k]: v } : s));

  // ── Stats ───────────────────────────────────────────────
  const stats = results ? {
    total:     results.length,
    edited:    results.filter(r => r.zygosityClass > 0).length,
    wt:        results.filter(r => r.zygosityClass === 0).length,
    homo:      results.filter(r => r.zygosityClass === 3).length,
    het:       results.filter(r => r.zygosityClass === 2).length,
    low:       results.filter(r => r.zygosityClass === 1).length,
    flagged:   results.filter(r => r.flag).length,
    eff: (() => {
      const conf = results.filter(r => !r.flag).length;
      const confEd = results.filter(r => r.zygosityClass > 0 && !r.flag).length;
      return conf > 0 ? Math.round(confEd / conf * 1000) / 10 : 0;
    })(),
    fs:        results.filter(r => r.isFrameshift).length,
    deletions: results.filter(r => r.editType === 'Deletion').length,
    insertions:results.filter(r => r.editType === 'Insertion').length,
  } : null;

  const TABS = ['🔬 Setup','📊 Dashboard','🧬 Explorer','🤖 AI Insights','🧠 ML Model','📄 Report','📈 Chromatogram','🔍 CAPS Designer'];

  // ── Shared styles ────────────────────────────────────────
  const card   = { background: C.card, border: `1px solid ${C.border}`, borderRadius: 16, padding: 24, marginBottom: 20 };
  const inp    = { background: C.card2, border: `1px solid ${C.border}`, borderRadius: 8, color: C.text, padding: '9px 13px', fontSize: 13, width: '100%', fontFamily: 'inherit', outline: 'none' };
  const lbl    = { fontSize: 10, color: C.muted, letterSpacing: 1, textTransform: 'uppercase', marginBottom: 5, display: 'block' };
  const btn    = { background: 'linear-gradient(135deg, #00d4aa, #00a882)', color: '#000', border: 'none', borderRadius: 10, padding: '11px 26px', fontSize: 13, fontWeight: 700, fontFamily: 'inherit', cursor: 'pointer' };
  const btnSec = { background: 'transparent', color: C.cyan, border: `1px solid ${C.cyan}`, borderRadius: 10, padding: '9px 18px', fontSize: 12, fontWeight: 600, fontFamily: 'inherit', cursor: 'pointer' };
  const th     = { background: C.card2, color: C.muted, padding: '9px 13px', textAlign: 'left', fontSize: 10, letterSpacing: 1, textTransform: 'uppercase', borderBottom: `1px solid ${C.border}` };
  const td     = { padding: '9px 13px', borderBottom: `1px solid rgba(30,58,95,0.4)`, verticalAlign: 'middle' };

  // Helper — is this the confirmed star sample?
  const isStarSample = r => r.zygosityClass === 3;
  const isWTSample   = r => r.id === (wtControlResult?.id);

  return (
    <div style={{ paddingTop: 70, minHeight: '100vh', background: C.bg, fontFamily: "'Segoe UI', system-ui, sans-serif" }}>
      <div style={{ maxWidth: 1320, margin: '0 auto', padding: '28px 28px' }}>

        {/* Header */}
        <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', marginBottom: 24, flexWrap: 'wrap', gap: 12 }}>
          <button onClick={onBack} style={btnSec}>← Back to Home</button>
          <div style={{ display: 'flex', alignItems: 'center', gap: 12 }}>
            <img src="/logo.png" alt="PlantaCRISPAI" style={{ width: 36, height: 36, objectFit: 'contain', filter: 'drop-shadow(0 0 6px rgba(0,212,170,0.4))' }} />
            <div>
              <div style={{ fontSize: 18, fontWeight: 800, background: 'linear-gradient(135deg, #00d4aa, #4d9fff)', WebkitBackgroundClip: 'text', WebkitTextFillColor: 'transparent' }}>PlantaCRISPAI Analysis</div>
              <div style={{ fontSize: 10, color: C.muted, letterSpacing: 1 }}>CRISPR GENOME EDITING PLATFORM</div>
            </div>
          </div>
          <div style={{ background: 'rgba(0,212,170,0.1)', border: '1px solid rgba(0,212,170,0.3)', color: C.cyan, padding: '4px 14px', borderRadius: 20, fontSize: 11, letterSpacing: 1 }}>v2.0</div>
        </div>

        {/* Tab bar */}
        <div style={{ display: 'flex', gap: 3, background: C.card, padding: 4, borderRadius: 12, marginBottom: 24, border: `1px solid ${C.border}`, overflowX: 'auto' }}>
          {TABS.map((t, i) => (
            <button key={i} onClick={() => navigateTo(i)} style={{
              flex: 1, minWidth: 80, padding: '10px 4px',
              background: tab === i ? C.cyan : 'transparent',
              color: tab === i ? '#000' : C.muted,
              border: 'none', borderRadius: 8, fontSize: 11,
              fontWeight: tab === i ? 700 : 500, fontFamily: 'inherit', cursor: 'pointer', transition: 'all 0.2s'
            }}>{t}</button>
          ))}
        </div>

        {/* ═══════════ TAB 0 — SETUP ═══════════ */}
        {tab === 0 && (
          <div>
            {/* Demo banner */}
            <div style={{ background: 'linear-gradient(135deg, rgba(0,212,170,0.08), rgba(77,159,255,0.08))', border: '1px solid rgba(0,212,170,0.3)', borderRadius: 14, padding: '18px 24px', marginBottom: 22, display: 'flex', justifyContent: 'space-between', alignItems: 'center', flexWrap: 'wrap', gap: 12 }}>
              <div>
                <div style={{ color: C.cyan, fontWeight: 700, marginBottom: 3 }}>🧬 Quick Demo — Rice OsAAP3 Editing Data</div>
                <div style={{ color: C.muted, fontSize: 12 }}>10 real OsAAP3 T1 AB1 reads · WT, Low, Het, Homo★ · validated by ICE+TIDE</div>
              </div>
              <button style={btn} onClick={loadDemo}>Load Demo Data →</button>
            </div>

            <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: 22 }}>

              {/* Config */}
              <div style={card}>
                <SectionTitle>⚙️ Experiment Configuration</SectionTitle>
                <div style={{ display: 'flex', flexDirection: 'column', gap: 14 }}>
                  <div>
                    <label style={lbl}>Target Gene Name</label>
                    <input style={inp} value={cfg.geneName} onChange={e => setCfg(c => ({ ...c, geneName: e.target.value }))} placeholder="e.g., OsAAP3" />
                  </div>
                  <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: 12 }}>
                    <div>
                      <label style={lbl}>Organism</label>
                      <select style={inp} value={cfg.organism} onChange={e => setCfg(c => ({ ...c, organism: e.target.value }))}>
                        {ORGANISM_OPTIONS.map(o => <option key={o.value} value={o.value}>{o.label}</option>)}
                      </select>
                    </div>
                    <div>
                      <label style={lbl}>Ploidy Level</label>
                      <select style={inp} value={cfg.ploidy} onChange={e => setCfg(c => ({ ...c, ploidy: e.target.value }))}>
                        {PLOIDY_OPTIONS.map(o => <option key={o.value} value={o.value}>{o.label}</option>)}
                      </select>
                    </div>
                  </div>
                  <div>
                    <label style={lbl}>Plant Generation</label>
                    <div style={{ display: 'flex', gap: 6 }}>
                      {['T0','T1','T2','T3'].map(g => (
                        <button key={g} onClick={() => setCfg(c => ({ ...c, generation: g }))} style={{ flex: 1, padding: '8px 0', background: cfg.generation === g ? C.cyan : 'transparent', color: cfg.generation === g ? '#000' : C.muted, border: `1px solid ${cfg.generation === g ? C.cyan : C.border}`, borderRadius: 8, fontSize: 13, fontWeight: 600, fontFamily: 'inherit', cursor: 'pointer' }}>{g}</button>
                      ))}
                    </div>
                  </div>
                  <div>
                    {/* ── Multi-gRNA Panel (Phase 3) ── */}
                    <div style={{ display:'flex', justifyContent:'space-between', alignItems:'center', marginBottom:8 }}>
                      <label style={lbl}>Guide RNA Sequence{grnaList.length > 1 ? 's' : ''} (20 bp, no PAM)</label>
                      <div style={{ display:'flex', gap:6 }}>
                        {grnaList.length > 1 && (
                          <span style={{ fontSize:10, color:'#BF5AF2', background:'rgba(191,90,242,0.1)', border:'1px solid rgba(191,90,242,0.3)', borderRadius:4, padding:'2px 8px', fontWeight:600 }}>
                            ✦ Multiplex — {grnaList.length} gRNAs
                          </span>
                        )}
                        <button
                          onClick={() => setGrnaList(l => [...l, { label:`gRNA ${l.length+1}`, sequence:'' }])}
                          style={{ fontSize:11, color:C.cyan, background:'transparent', border:`1px solid ${C.border}`, borderRadius:6, padding:'3px 10px', cursor:'pointer' }}>
                          + Add gRNA
                        </button>
                      </div>
                    </div>

                    {/* gRNA rows */}
                    {grnaList.map((g, gi) => {
                      const rc   = s => s.split('').reverse().map(b=>({A:'T',T:'A',G:'C',C:'G'}[b]||b)).join('');
                      const ci   = cfg.ref && g.sequence.length >= 18 ? deriveCutPosition(cfg.ref, g.sequence) : null;
                      const fwd  = ci && cfg.ref.toUpperCase().includes(g.sequence);
                      const rev  = ci && !fwd && cfg.ref.toUpperCase().includes(rc(g.sequence));
                      const valid = g.sequence.length >= 18 && g.sequence.length <= 22;
                      return (
                        <div key={gi} style={{ display:'flex', gap:8, alignItems:'flex-start', marginBottom:8 }}>
                          {/* Label pill */}
                          <input
                            value={g.label}
                            onChange={e => setGrnaList(l => l.map((x,i)=>i===gi?{...x,label:e.target.value}:x))}
                            style={{ ...inp, width:80, flexShrink:0, fontSize:11, padding:'8px 8px', fontWeight:600, color:C.cyan }}
                            placeholder={`gRNA ${gi+1}`}
                          />
                          {/* Sequence input */}
                          <div style={{ flex:1 }}>
                            <input
                              style={{ ...inp, fontFamily:'monospace', fontSize:12,
                                borderColor: g.sequence.length > 0 ? (valid ? (fwd?'rgba(0,212,170,0.5)':rev?'rgba(191,90,242,0.5)':'rgba(240,192,64,0.5)') : 'rgba(255,77,109,0.4)') : C.border }}
                              value={g.sequence}
                              onChange={e => {
                                const seq = e.target.value.toUpperCase().replace(/[^ATGCN]/g,'');
                                setGrnaList(l => l.map((x,i)=>i===gi?{...x,sequence:seq}:x));
                                // keep cfg.grna in sync with first gRNA for backwards compat
                                if (gi === 0) setCfg(c => ({...c, grna: seq}));
                              }}
                              placeholder="e.g., CGGCGTGCTCTCTCTCGCAT"
                            />
                            {g.sequence.length > 0 && (
                              <div style={{ marginTop:3, fontSize:10, display:'flex', gap:8, flexWrap:'wrap' }}>
                                <span style={{ color: valid ? C.cyan : C.red }}>
                                  {valid ? '✓ Valid length' : `⚠ ${g.sequence.length} bp`}
                                </span>
                                {ci && fwd && <span style={{ color:C.cyan }}>✓ (+) strand · cut ~{ci.cutPosition} bp</span>}
                                {ci && rev && <span style={{ color:'#BF5AF2' }}>✓ (−) strand · cut ~{ci.cutPosition} bp</span>}
                                {ci && !fwd && !rev && g.sequence.length >= 18 && <span style={{ color:C.gold }}>⚠ Not found in amplicon</span>}
                                {/* Show inter-cut distance when ≥2 gRNAs found */}
                                {gi > 0 && ci && (fwd||rev) && (() => {
                                  const prev = grnaList.slice(0,gi).findLast(x => x.sequence.length >= 18);
                                  if (!prev || !cfg.ref) return null;
                                  const ci0 = deriveCutPosition(cfg.ref, prev.sequence);
                                  if (!ci0.cutPosition) return null;
                                  const dist = Math.abs(ci.cutPosition - ci0.cutPosition);
                                  return <span style={{ color:'#BF5AF2', fontWeight:600 }}>↔ {dist} bp from {prev.label||`gRNA ${gi}`}</span>;
                                })()}
                              </div>
                            )}
                          </div>
                          {/* Remove button (only if > 1 gRNA) */}
                          {grnaList.length > 1 && (
                            <button
                              onClick={() => {
                                const next = grnaList.filter((_,i)=>i!==gi);
                                setGrnaList(next);
                                if (gi === 0) setCfg(c=>({...c, grna: next[0]?.sequence||''}));
                              }}
                              style={{ background:'rgba(255,77,109,0.1)', border:'1px solid rgba(255,77,109,0.3)', borderRadius:6, color:C.red, fontSize:14, width:30, height:30, cursor:'pointer', flexShrink:0, marginTop:2 }}
                              title="Remove this gRNA">
                              ✕
                            </button>
                          )}
                        </div>
                      );
                    })}

                    {/* Multiplex summary */}
                    {grnaList.length > 1 && cfg.ref && (() => {
                      const found = grnaList.filter(g => g.sequence.length >= 18 && cfg.ref.toUpperCase().includes(g.sequence));
                      if (found.length < 2) return null;
                      const cuts  = found.map(g => deriveCutPosition(cfg.ref, g.sequence).cutPosition).sort((a,b)=>a-b);
                      const span  = cuts[cuts.length-1] - cuts[0];
                      return (
                        <div style={{ background:'rgba(191,90,242,0.08)', border:'1px solid rgba(191,90,242,0.25)', borderRadius:8, padding:'8px 12px', fontSize:11, color:'#BF5AF2', marginTop:4 }}>
                          ✦ <strong>Dual-cut configuration:</strong> {found.length} gRNAs active · cut sites span <strong>{span} bp</strong> ·
                          potential inter-cut deletion of ~{span} bp will be auto-detected per sample
                        </div>
                      );
                    })()}
                  </div>
                  <div>
                    <label style={lbl}>Reference Amplicon (wild-type sequence)</label>
                    {/* Upload zone for FASTA or .ab1 */}
                    <div style={{ display: 'flex', gap: 8, marginBottom: 8 }}>
                      <div
                        style={{ flex: 1, border: `1px dashed ${C.border}`, borderRadius: 8, padding: '8px 12px', textAlign: 'center', cursor: 'pointer', background: C.card2, fontSize: 11, color: C.muted }}
                        onClick={() => refFileRef.current?.click()}
                        onDragOver={e => e.preventDefault()}
                        onDrop={e => { e.preventDefault(); handleRefFile({ dataTransfer: e.dataTransfer }); }}>
                        <input ref={refFileRef} type="file" accept=".fasta,.fa,.txt,.ab1,.seq" onChange={handleRefFile} style={{ display: 'none' }} />
                        📁 Upload reference .fasta or .ab1
                      </div>
                      <div style={{ fontSize: 10, color: C.dim, display: 'flex', alignItems: 'center', paddingLeft: 4 }}>or paste below ↓</div>
                    </div>
                    <textarea style={{ ...inp, minHeight: 80, resize: 'vertical', fontFamily: 'monospace', fontSize: 11, lineHeight: 1.6 }}
                      value={cfg.ref} onChange={e => setCfg(c => ({ ...c, ref: e.target.value.toUpperCase().replace(/\s/g,'') }))}
                      placeholder="Paste wild-type reference amplicon sequence here (or upload FASTA above)..." />
                    {cfg.ref && <div style={{ marginTop: 5, fontSize: 11, color: C.muted }}>{cfg.ref.length} bp · GC: {((cfg.ref.match(/[GC]/g)||[]).length / cfg.ref.length * 100).toFixed(1)}%{grnaList.filter(g=>g.sequence.length>=18 && cfg.ref.toUpperCase().includes(g.sequence)).length > 0 ? ` · ✓ ${grnaList.filter(g=>g.sequence.length>=18&&cfg.ref.toUpperCase().includes(g.sequence)).length}/${grnaList.filter(g=>g.sequence.length>=18).length} gRNA(s) found` : ''}</div>}
                  </div>
                </div>
              </div>

              {/* Sample Manager */}
              <div style={card}>
                <div style={{ fontSize: 14, fontWeight: 700, color: C.cyan, marginBottom: 14, display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                  <span>🔬 Sample Manager</span>
                  <span style={{ fontSize: 11, color: C.muted }}>{samples.length} sample{samples.length !== 1 ? 's' : ''}</span>
                </div>

                {/* AB1 Upload Panel — with 4-mode trimming controls (Phase 2) */}
                <AB1UploadPanel
                  generation={cfg.generation}
                  ploidy={cfg.ploidy}
                  cutSite={cfg.ref && grnaList[0]?.sequence ? deriveCutPosition(cfg.ref, grnaList[0].sequence).cutPosition : null}
                  onSamplesLoaded={(newSamples) => {
                    setSamples(prev => {
                      // Merge: avoid duplicates by id, new samples override existing
                      const existingIds = new Set(prev.map(s => s.id));
                      const fresh = newSamples.filter(s => !existingIds.has(s.id));
                      return [...prev, ...fresh];
                    });
                  }}
                />

                <div style={{ background: 'rgba(77,159,255,0.06)', border: '1px solid rgba(77,159,255,0.2)', borderRadius: 8, padding: '8px 12px', marginBottom: 10, marginTop: 10, fontSize: 11, color: C.blue }}>
                  💡 Upload ALL samples here including your WT control (e.g. MTU1010). Any sample whose sequence matches the reference amplicon will automatically be identified as wild-type.
                </div>

                <button style={{ ...btnSec, width: '100%', marginBottom: 12 }} onClick={addSample}>+ Add Sample Manually</button>

                <div style={{ maxHeight: 340, overflowY: 'auto', display: 'flex', flexDirection: 'column', gap: 8 }}>
                  {samples.length === 0 && (
                    <div style={{ textAlign: 'center', color: C.dim, padding: 28, fontSize: 12 }}>No samples yet. Load demo or upload files above.</div>
                  )}
                  {samples.map((s, i) => {
                    const isWT = /MTU|WT|wildtype|control/i.test(s.id);
                    return (
                      <div key={i} style={{ background: isWT ? 'rgba(0,212,170,0.06)' : C.card2, border: `1px solid ${isWT ? 'rgba(0,212,170,0.4)' : C.border}`, borderRadius: 10, padding: 10 }}>
                        <div style={{ display: 'flex', gap: 8, marginBottom: 6, alignItems: 'center' }}>
                          {isWT && <span style={{ fontSize: 14 }}>🧬</span>}
                          <input style={{ ...inp, flex: 1, padding: '5px 10px', fontSize: 12 }} value={s.id} onChange={e => updSample(i, 'id', e.target.value)} />
                          <select style={{ ...inp, width: 68, padding: '5px 6px', fontSize: 12 }} value={s.generation} onChange={e => updSample(i, 'generation', e.target.value)}>
                            <option>T0</option><option>T1</option><option>T2</option><option>T3</option>
                          </select>
                          <button onClick={() => delSample(i)} style={{ background: 'rgba(255,77,109,0.15)', border: '1px solid rgba(255,77,109,0.3)', color: C.red, borderRadius: 7, padding: '5px 9px', cursor: 'pointer', fontFamily: 'inherit' }}>✕</button>
                        </div>
                        {isWT && <div style={{ fontSize: 10, color: C.cyan, marginBottom: 5, paddingLeft: 2 }}>⭐ WT Reference Control — used for all alignments</div>}
                        {s.parseError && <div style={{ fontSize: 10, color: C.red, marginBottom: 5 }}>⚠ {s.parseError}</div>}
                        {/* Trim mode indicator for AB1-sourced samples */}
                        {s.source === 'ab1' && s.trimMode && s.trimMode !== 'standard' && (
                          <div style={{ fontSize: 10, color: '#BF5AF2', marginBottom: 4 }}>
                            ✂️ Trim: {s.trimMode}{s.trimMode === 'manual' ? ` (${s.trimStart}–${s.trimEnd})` : ''}{s.trimWarning ? ' ⚠' : ''}
                          </div>
                        )}
                        <textarea style={{ ...inp, minHeight: 40, fontSize: 11, fontFamily: 'monospace', lineHeight: 1.5, resize: 'none' }}
                          value={s.sequence} onChange={e => updSample(i, 'sequence', e.target.value.toUpperCase().replace(/\s/g, ''))}
                          placeholder="Sequence (auto-filled from .ab1)" rows={2} />
                        {s.sequence && <div style={{ fontSize: 10, color: C.dim, marginTop: 3 }}>{s.sequence.length} bp{s.source === 'ab1' ? ` · from .ab1` : ''}{s.avgQuality ? ` · avg Q${s.avgQuality.toFixed(0)}` : ''}{s.rawSeqLength && s.rawSeqLength !== s.sequence.length ? ` · raw ${s.rawSeqLength} bp` : ''}</div>}
                      </div>
                    );
                  })}
                </div>
              </div>
            </div>

            {/* Run button */}
            <div style={{ textAlign: 'center', marginTop: 8 }}>
              {busy ? (
                <div style={{ ...card, display: 'inline-block', minWidth: 400, textAlign: 'left' }}>
                  <div style={{ color: C.cyan, marginBottom: 12, textAlign: 'center', fontWeight: 700 }}>🤖 AI Analysis Running...</div>
                  <ProgressBar value={progress} />
                  <div style={{ display: 'flex', justifyContent: 'space-between', marginTop: 8, fontSize: 11, color: C.muted }}>
                    <span>{progressMsg}</span><span>{progress}%</span>
                  </div>
                </div>
              ) : (
                <button style={{ ...btn, fontSize: 15, padding: '14px 52px', opacity: (!grnaList.some(g=>g.sequence.length>=18) || !cfg.ref || !samples.length) ? 0.35 : 1 }}
                  onClick={runAnalysis} disabled={!grnaList.some(g=>g.sequence.length>=18) || !cfg.ref || !samples.length}>
                  🚀 Run Analysis ({samples.length} sample{samples.length !== 1 ? 's' : ''}{grnaList.filter(g=>g.sequence.length>=18).length > 1 ? ` · ${grnaList.filter(g=>g.sequence.length>=18).length} gRNAs` : ''})
                </button>
              )}
              {!busy && (!grnaList.some(g=>g.sequence.length>=18) || !cfg.ref || !samples.length) && (
                <div style={{ marginTop: 8, color: C.muted, fontSize: 11 }}>Required: gRNA sequence + reference amplicon + at least 1 sample</div>
              )}
            </div>
          </div>
        )}

        {/* ═══════════ TAB 1 — DASHBOARD ═══════════ */}
        {tab === 1 && results && (
          <div>
            {/* WT Control banner */}
            {wtControlResult && (
              <div style={{ background: 'rgba(0,212,170,0.06)', border: '1px solid rgba(0,212,170,0.3)', borderRadius: 12, padding: '12px 20px', marginBottom: 18, display: 'flex', alignItems: 'center', gap: 16 }}>
                <div style={{ fontSize: 20 }}>🧬</div>
                <div>
                  <div style={{ color: C.cyan, fontWeight: 700, fontSize: 13 }}>WT Reference Control: {wtControlResult.id}</div>
                  <div style={{ color: C.muted, fontSize: 11 }}>Identified as wild-type by alignment · R²={(wtControlResult.r2||0).toFixed(3)} · Identity={wtControlResult.sequenceIdentity ? (wtControlResult.sequenceIdentity*100).toFixed(1)+'%' : '—'} · ICE={wtControlResult.iceScore||0}%</div>
                </div>
              </div>
            )}

            {/* Stat cards */}
            <div style={{ display: 'grid', gridTemplateColumns: 'repeat(4, 1fr)', gap: 16, marginBottom: 22 }}>
              <StatCard label="Editing Efficiency" value={`${stats.eff}%`} color={C.cyan}
                onClick={() => { setExplorerFilter('edited'); navigateTo(2); }} />
              <StatCard label="Edited Plants" value={`${stats.edited}/${stats.total}`} color={C.blue}
                onClick={() => { setExplorerFilter('edited'); navigateTo(2); }} />
              <StatCard label="Frameshift Mutations" value={stats.fs} color={C.red}
                onClick={() => { setExplorerFilter('fs'); navigateTo(2); }} />
              <StatCard label="Wild-type Controls" value={stats.wt} color={C.muted}
                onClick={() => { setExplorerFilter('wt'); navigateTo(2); }} />
            </div>

            {/* Second row */}
            <div style={{ display: 'grid', gridTemplateColumns: 'repeat(3, 1fr)', gap: 16, marginBottom: 22 }}>
              <StatCard label="Homozygous ★" value={stats.homo} color={C.cyan} sub="Target plants for advancement"
                onClick={() => { setExplorerFilter('homo'); navigateTo(2); }} />
              <StatCard label="Heterozygous" value={stats.het} color={C.purple} sub="Self-pollinate → T2"
                onClick={() => { setExplorerFilter('het'); navigateTo(2); }} />
              <StatCard label="Needs Verification" value={stats.flagged} color={C.gold} sub="Re-sequence recommended"
                onClick={() => { setExplorerFilter('flagged'); navigateTo(2); }} />
            </div>

            {/* Charts row */}
            <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: 22, marginBottom: 22 }}>
              <div style={card}>
                <SectionTitle>🎯 Edit Type Distribution</SectionTitle>
                {[
                  { label: 'Homozygous Edited ★', val: stats.homo,      color: C.cyan   },
                  { label: 'Heterozygous',         val: stats.het,       color: C.purple },
                  { label: 'Low / Mosaic Edit',    val: stats.low,       color: C.gold   },
                  { label: 'Wild-type Segregant',  val: stats.wt,        color: C.dim    },
                  { label: 'Deletions',            val: stats.deletions, color: C.red    },
                  { label: 'Insertions',           val: stats.insertions,color: C.green  },
                ].map(item => (
                  <div key={item.label} style={{ marginBottom: 12 }}>
                    <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: 4 }}>
                      <span style={{ color: item.color, fontSize: 12 }}>● {item.label}</span>
                      <span style={{ color: C.muted, fontSize: 11 }}>{item.val} ({stats.total > 0 ? Math.round(item.val/stats.total*100) : 0}%)</span>
                    </div>
                    <ProgressBar value={stats.total > 0 ? item.val/stats.total*100 : 0} color={item.color} />
                  </div>
                ))}
              </div>

              <div style={card}>
                <SectionTitle>🔬 ML Confidence Scores</SectionTitle>
                <div style={{ maxHeight: 320, overflowY: 'auto' }}>
                  {results.map((r, i) => (
                    <div key={i} style={{ marginBottom: 8, display: 'flex', alignItems: 'center', gap: 10 }}>
                      <div style={{ width: 100, fontSize: 11, color: isStarSample(r) ? C.cyan : C.muted, fontWeight: isStarSample(r) ? 700 : 400, flexShrink: 0 }}>
                        {isStarSample(r) ? '★ ' : ''}{r.id}
                      </div>
                      <ZygBadge cls={r.zygosityClass} />
                      <ScoreBar score={r.mlScore} width={80} />
                    </div>
                  ))}
                </div>
              </div>
            </div>

            {/* Full results table */}
            <div style={card}>
              <SectionTitle>📋 All Samples</SectionTitle>
              <div style={{ overflowX: 'auto' }}>
                <table style={{ width: '100%', borderCollapse: 'collapse', fontSize: 12 }}>
                  <thead>
                    <tr>{['Sample ID','Gen.','Edit Type','Size','Zygosity','Strand','R²','ICE','Protein Effect','ML Score','Flag',''].map(h => <th key={h} style={th}>{h}</th>)}</tr>
                  </thead>
                  <tbody>
                    {results.map((r, i) => (
                      <tr key={i} style={{ cursor: 'pointer', background: selectedSample?.id === r.id ? 'rgba(0,212,170,0.04)' : 'transparent' }}
                        onClick={() => navigateTo(2, r)}>
                        <td style={td}>
                          <code style={{ color: isStarSample(r) ? C.cyan : C.blue, fontWeight: isStarSample(r) ? 700 : 400 }}>
                            {isStarSample(r) ? '★ ' : ''}{r.id}
                          </code>
                        </td>
                        <td style={td}><Chip color={C.purple} small>{r.generation}</Chip></td>
                        <td style={td}><EditChip editType={r.editType} indelSize={r.indelSize} /></td>
                        <td style={td}>{r.indelSize > 0 ? <code style={{ color: C.gold }}>{r.indelSize} bp</code> : <span style={{ color: C.dim }}>—</span>}</td>
                        <td style={td}><ZygBadge cls={r.zygosityClass} /></td>
                        <td style={td}>
                          <span title={r.wasFlipped ? 'Reverse-primer read — auto-corrected to forward orientation' : r.grnaStrand === '-' ? 'gRNA targets minus strand' : 'Forward primer read'}
                            style={{
                              fontSize: 10, fontFamily: 'monospace', fontWeight: 700,
                              color: (r.readStrand === '-' || r.grnaStrand === '-') ? '#BF5AF2' : C.cyan,
                              background: (r.readStrand === '-' || r.grnaStrand === '-') ? 'rgba(191,90,242,0.1)' : 'rgba(0,212,170,0.1)',
                              border: `1px solid ${(r.readStrand === '-' || r.grnaStrand === '-') ? 'rgba(191,90,242,0.35)' : 'rgba(0,212,170,0.35)'}`,
                              borderRadius: 3, padding: '1px 5px', cursor: 'help',
                            }}>
                            {r.wasFlipped ? '↺ RC' : r.grnaStrand === '-' ? '(−)' : '(+)'}
                          </span>
                        </td>
                        <td style={td}><span style={{ color: (r.r2||0) > 0.85 ? C.cyan : C.gold, fontSize: 11, fontFamily: 'monospace' }}>{(r.r2||0).toFixed(2)}</span></td>
                        <td style={td}><span style={{ color: C.muted, fontSize: 11 }}>{r.iceScore||0}%</span></td>
                        <td style={{ ...td, maxWidth: 200 }}><span style={{ fontSize: 11, color: r.isFrameshift ? C.red : r.indelSize > 0 ? C.gold : C.dim }}>{r.proteinEffect}</span></td>
                        <td style={td}><ScoreBar score={r.mlScore} /></td>
                        <td style={td}>
                          {r.flag ? (
                            <span title={r.flag} style={{ fontSize: 10, color: r.flag.includes('Borderline') ? C.blue : C.gold, cursor: 'help' }}>
                              {r.flag.includes('Borderline') ? 'ℹ️ Borderline' : r.flag.includes('Poor') ? '⚠️ Quality' : r.flag.includes('strand') ? '⚠️ Strand' : '⚠️ Verify'}
                            </span>
                          ) : <span style={{ color: C.dim, fontSize: 10 }}>—</span>}
                        </td>
                        <td style={td}><button style={{ ...btnSec, padding: '4px 12px', fontSize: 11 }} onClick={e => { e.stopPropagation(); navigateTo(2, r); }}>View →</button></td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>
          </div>
        )}

        {/* ═══════════ TAB 2 — EXPLORER ═══════════ */}
        {tab === 2 && results && (
          <div style={{ display: 'grid', gridTemplateColumns: '280px 1fr', gap: 22 }}>
            {/* Sample list */}
            <div style={{ ...card, padding: 14 }}>
              <SectionTitle>🔎 Select Sample</SectionTitle>

              {/* Filter bar — shown when navigated from Dashboard */}
              <div style={{ marginBottom: 10 }}>
                <div style={{ display: 'flex', flexWrap: 'wrap', gap: 4 }}>
                  {[
                    { id:'all',     label:'All',          color: C.muted },
                    { id:'edited',  label:'Edited',       color: C.cyan  },
                    { id:'homo',    label:'Homo',         color: C.cyan  },
                    { id:'het',     label:'Het',          color: C.purple},
                    { id:'wt',      label:'WT',           color: C.blue  },
                    { id:'fs',      label:'Frameshift',   color: C.red   },
                    { id:'flagged', label:'⚠ Flagged',    color: C.gold  },
                  ].map(f => (
                    <button key={f.id} onClick={() => setExplorerFilter(f.id)} style={{
                      background: explorerFilter === f.id ? `${f.color}22` : 'transparent',
                      border: `1px solid ${explorerFilter === f.id ? f.color : C.border}`,
                      color: explorerFilter === f.id ? f.color : C.muted,
                      borderRadius: 6, padding: '3px 8px', fontSize: 10,
                      cursor: 'pointer', fontWeight: explorerFilter === f.id ? 700 : 400,
                      transition: 'all 0.15s',
                    }}>{f.label}</button>
                  ))}
                </div>
                {explorerFilter !== 'all' && (
                  <div style={{ fontSize: 9, color: C.muted, marginTop: 4 }}>
                    Filtering by: <b style={{ color: C.cyan }}>{explorerFilter}</b>
                    {' '}— showing {results.filter(r => {
                      if (explorerFilter === 'edited')  return r.zygosityClass > 0;
                      if (explorerFilter === 'homo')    return r.zygosityClass === 3;
                      if (explorerFilter === 'het')     return r.zygosityClass === 2;
                      if (explorerFilter === 'wt')      return r.zygosityClass === 0;
                      if (explorerFilter === 'fs')      return r.isFrameshift;
                      if (explorerFilter === 'flagged') return !!r.flag;
                      return true;
                    }).length} of {results.length} samples
                  </div>
                )}
              </div>

              {results.filter(r => {
                if (explorerFilter === 'edited')  return r.zygosityClass > 0;
                if (explorerFilter === 'homo')    return r.zygosityClass === 3;
                if (explorerFilter === 'het')     return r.zygosityClass === 2;
                if (explorerFilter === 'wt')      return r.zygosityClass === 0;
                if (explorerFilter === 'fs')      return r.isFrameshift;
                if (explorerFilter === 'flagged') return !!r.flag;
                return true;
              }).map((r, i) => (
                <div key={i} onClick={() => setSelectedSample(r)} style={{ padding: '10px 12px', borderRadius: 10, cursor: 'pointer', marginBottom: 6, background: selectedSample?.id === r.id ? 'rgba(0,212,170,0.1)' : C.card2, border: `1px solid ${selectedSample?.id === r.id ? C.cyan : C.border}`, transition: 'all 0.15s' }}>
                  <div style={{ fontWeight: 600, fontSize: 12, color: isStarSample(r) ? C.cyan : C.text }}>
                    {isStarSample(r) ? '★ ' : ''}{r.id}
                  </div>
                  <div style={{ fontSize: 10, color: C.muted, marginTop: 2 }}>{r.generation} · <EditChip editType={r.editType} indelSize={r.indelSize} /></div>
                  <div style={{ marginTop: 6, display: 'flex', gap: 6, alignItems: 'center', flexWrap: 'wrap' }}>
                    <ZygBadge cls={r.zygosityClass} />
                    {isWTSample(r) && <Chip color={C.cyan} small>Ref Control</Chip>}
                    <ScoreBar score={r.mlScore} width={50} />
                  </div>
                  {r.flag && <div style={{ fontSize: 9, color: C.gold, marginTop: 4 }}>⚠ {r.flag}</div>}
                </div>
              ))}
            </div>

            {/* Detail panel */}
            {selectedSample ? (
              <div style={{ ...card, border: `1px solid ${isStarSample(selectedSample) ? 'rgba(0,212,170,0.5)' : C.border}` }}>
                {/* Header */}
                <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: 20, flexWrap: 'wrap', gap: 12 }}>
                  <div>
                    <div style={{ fontSize: 22, fontWeight: 800, color: isStarSample(selectedSample) ? C.cyan : C.text, marginBottom: 6 }}>
                      {isStarSample(selectedSample) ? '★ ' : ''}{selectedSample.id}
                    </div>
                    <div style={{ display: 'flex', gap: 8, flexWrap: 'wrap' }}>
                      <Chip color={C.blue}>{selectedSample.generation}</Chip>
                      <Chip color={C.purple}>{selectedSample.ploidy}</Chip>
                      <EditChip editType={selectedSample.editType} indelSize={selectedSample.indelSize} />
                      <ZygBadge cls={selectedSample.zygosityClass} />
                    </div>
                  </div>
                  <div style={{ textAlign: 'right' }}>
                    <div style={{ fontSize: 10, color: C.muted, marginBottom: 5 }}>ML CONFIDENCE</div>
                    <ScoreBar score={selectedSample.mlScore} width={80} />
                  </div>
                </div>

                {/* Star banner for confirmed homo */}
                {isStarSample(selectedSample) && (
                  <div style={{ background: 'rgba(0,212,170,0.1)', border: '1px solid rgba(0,212,170,0.4)', borderRadius: 10, padding: '10px 16px', marginBottom: 16, color: C.cyan, fontSize: 12, fontWeight: 600 }}>
                    ⭐ Confirmed Homozygous Edited — Priority advance to T2. Self-pollinate immediately.
                  </div>
                )}

                {selectedSample.flag && (
                  <div style={{
                    background: selectedSample.flag.includes('Borderline') ? 'rgba(77,159,255,0.08)' : 'rgba(240,192,64,0.08)',
                    border: `1px solid ${selectedSample.flag.includes('Borderline') ? 'rgba(77,159,255,0.3)' : 'rgba(240,192,64,0.3)'}`,
                    borderRadius: 10, padding: '10px 16px', marginBottom: 16,
                    color: selectedSample.flag.includes('Borderline') ? C.blue : C.gold, fontSize: 12
                  }}>
                    {selectedSample.flag.includes('Borderline') ? 'ℹ️' : '⚠️'} {selectedSample.flag}
                  </div>
                )}

                {/* Zygosity + Effect */}
                <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: 14, marginBottom: 20 }}>
                  <div style={{ background: C.card2, border: `1px solid ${C.border}`, borderRadius: 12, padding: '16px 18px', textAlign: 'center' }}>
                    <div style={{ fontSize: 16, fontWeight: 700, color: ZYG[selectedSample.zygosityClass]?.color || C.muted, marginBottom: 4 }}>{ZYG[selectedSample.zygosityClass]?.label}</div>
                    <div style={{ fontSize: 10, color: C.muted, letterSpacing: 1 }}>ZYGOSITY</div>
                    <div style={{ marginTop: 10 }}><ScoreBar score={selectedSample.mlScore} width={100} /></div>
                  </div>
                  <div style={{ background: C.card2, border: `1px solid ${C.border}`, borderRadius: 12, padding: '16px 18px', textAlign: 'center' }}>
                    <div style={{ fontSize: 13, fontWeight: 700, color: selectedSample.isFrameshift ? C.red : selectedSample.indelSize > 0 ? C.gold : C.muted, marginBottom: 4, lineHeight: 1.5 }}>{selectedSample.proteinEffect}</div>
                    <div style={{ fontSize: 10, color: C.muted, letterSpacing: 1 }}>PREDICTED PROTEIN EFFECT</div>
                  </div>
                </div>

                {/* Metrics grid */}
                <div style={{ display: 'grid', gridTemplateColumns: 'repeat(4, 1fr)', gap: 10, marginBottom: 12 }}>
                  {[
                    { l: 'ICE Score',    v: `${selectedSample.iceScore||0}%`,                       c: C.blue   },
                    { l: 'R² Overall',   v: (selectedSample.r2||0).toFixed(3),                       c: (selectedSample.r2||0) > 0.85 ? C.cyan : C.gold },
                    { l: 'R² Upstream',  v: (selectedSample.r2up||0).toFixed(3),                     c: C.muted  },
                    { l: 'R² Downstream',v: (selectedSample.r2dn||0).toFixed(3),                     c: C.muted  },
                    { l: 'Seq. Identity',v: selectedSample.sequenceIdentity ? (selectedSample.sequenceIdentity*100).toFixed(1)+'%' : '—', c: C.muted },
                    { l: 'Indel Size',   v: selectedSample.indelSize > 0 ? `${selectedSample.indelSize} bp` : '—', c: C.gold },
                    { l: 'Frameshift',   v: selectedSample.isFrameshift ? 'Yes' : 'No',              c: selectedSample.isFrameshift ? C.red : C.muted },
                    { l: 'Edit Direction',v: selectedSample.lenDiff > 0 ? `+${selectedSample.lenDiff} bp ins.` : selectedSample.lenDiff < 0 ? `${selectedSample.lenDiff} bp del.` : 'None', c: C.muted },
                  ].map((m, i) => (
                    <div key={i} style={{ background: C.card2, border: `1px solid ${C.border}`, borderRadius: 8, padding: '10px 12px' }}>
                      <div style={{ fontSize: 10, color: C.muted, marginBottom: 4 }}>{m.l}</div>
                      <div style={{ fontSize: 14, fontWeight: 700, color: m.c, fontFamily: 'monospace' }}>{m.v}</div>
                    </div>
                  ))}
                </div>

                {/* Strand information row */}
                <div style={{ display: 'flex', gap: 8, marginBottom: 18, flexWrap: 'wrap' }}>
                  {/* Read strand badge */}
                  <div style={{ display: 'flex', alignItems: 'center', gap: 6, background: C.card2, border: `1px solid ${C.border}`, borderRadius: 8, padding: '6px 12px' }}>
                    <span style={{ fontSize: 10, color: C.muted }}>READ STRAND</span>
                    <span style={{
                      fontFamily: 'monospace', fontSize: 12, fontWeight: 700,
                      color: selectedSample.readStrand === '-' ? '#BF5AF2' : C.cyan,
                      background: selectedSample.readStrand === '-' ? 'rgba(191,90,242,0.12)' : 'rgba(0,212,170,0.12)',
                      border: `1px solid ${selectedSample.readStrand === '-' ? 'rgba(191,90,242,0.4)' : 'rgba(0,212,170,0.4)'}`,
                      borderRadius: 4, padding: '2px 7px',
                    }}>
                      {selectedSample.readStrand === '-' ? '(−) Reverse primer' : '(+) Forward primer'}
                    </span>
                    {selectedSample.wasFlipped && (
                      <span style={{ fontSize: 10, color: '#BF5AF2', fontStyle: 'italic' }}>auto-corrected ✓</span>
                    )}
                  </div>
                  {/* gRNA strand badge */}
                  <div style={{ display: 'flex', alignItems: 'center', gap: 6, background: C.card2, border: `1px solid ${C.border}`, borderRadius: 8, padding: '6px 12px' }}>
                    <span style={{ fontSize: 10, color: C.muted }}>gRNA TARGET</span>
                    <span style={{
                      fontFamily: 'monospace', fontSize: 12, fontWeight: 700,
                      color: selectedSample.grnaStrand === '-' ? '#BF5AF2' : C.cyan,
                      background: selectedSample.grnaStrand === '-' ? 'rgba(191,90,242,0.12)' : 'rgba(0,212,170,0.12)',
                      border: `1px solid ${selectedSample.grnaStrand === '-' ? 'rgba(191,90,242,0.4)' : 'rgba(0,212,170,0.4)'}`,
                      borderRadius: 4, padding: '2px 7px',
                    }}>
                      {selectedSample.grnaStrand === '-' ? '(−) Minus strand' : '(+) Plus strand'}
                    </span>
                  </div>
                  {/* Ambiguous strand warning */}
                  {selectedSample.strandAmbiguous && (
                    <div style={{ display: 'flex', alignItems: 'center', gap: 6, background: 'rgba(240,192,64,0.08)', border: '1px solid rgba(240,192,64,0.3)', borderRadius: 8, padding: '6px 12px' }}>
                      <span style={{ fontSize: 11, color: C.gold }}>⚠ Strand ambiguous — verify primer orientation</span>
                    </div>
                  )}
                </div>

                {/* Multi-gRNA breakdown panel — only shown when multiple gRNAs */}
                {selectedSample.isMultiGrna && selectedSample.perGrnaResults?.length > 0 && (
                  <div style={{ background: 'rgba(191,90,242,0.06)', border: '1px solid rgba(191,90,242,0.25)', borderRadius: 12, padding: '14px 18px', marginBottom: 18 }}>
                    <div style={{ fontSize: 11, color: '#BF5AF2', letterSpacing: 1, fontWeight: 700, marginBottom: 12, textTransform: 'uppercase', display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                      <span>✦ Multi-gRNA Analysis — {selectedSample.perGrnaResults.length} guide RNAs</span>
                      {selectedSample.isDualCut && (
                        <span style={{ background: 'rgba(191,90,242,0.15)', border: '1px solid rgba(191,90,242,0.4)', borderRadius: 6, padding: '2px 10px', fontSize: 10 }}>
                          ✦ Dual cut · ~{selectedSample.interCutDeletion} bp inter-cut deletion possible
                        </span>
                      )}
                    </div>
                    <div style={{ display: 'grid', gridTemplateColumns: `repeat(${Math.min(selectedSample.perGrnaResults.length, 3)}, 1fr)`, gap: 10 }}>
                      {selectedSample.perGrnaResults.map((pg, pgi) => {
                        const PCOLS = ['#00D4AA','#BF5AF2','#4D9FFF','#F0C040','#FF4D6D'];
                        const col   = PCOLS[pgi % PCOLS.length];
                        const zyg   = {3:'Homozygous',2:'Heterozygous',1:'Low/Mosaic',0:'Wild-type'}[pg.zygosity?.class ?? 0];
                        const zyCol = {3:C.cyan,2:C.purple,1:C.gold,0:C.muted}[pg.zygosity?.class ?? 0];
                        return (
                          <div key={pgi} style={{ background: C.card2, border: `1px solid ${col}30`, borderTop: `2px solid ${col}`, borderRadius: 8, padding: '10px 12px' }}>
                            <div style={{ fontSize: 12, fontWeight: 700, color: col, marginBottom: 8 }}>{pg.label}</div>
                            <div style={{ fontSize: 10, color: C.muted, marginBottom: 2 }}>ZYGOSITY</div>
                            <div style={{ fontSize: 13, fontWeight: 700, color: zyCol, marginBottom: 6 }}>{zyg}</div>
                            <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: 4, fontSize: 10 }}>
                              <div><span style={{ color: C.muted }}>ICE: </span><span style={{ color: col, fontWeight: 600 }}>{pg.iceScore}%</span></div>
                              <div><span style={{ color: C.muted }}>R²↓: </span><span style={{ color: col, fontWeight: 600 }}>{(pg.r2dn||0).toFixed(3)}</span></div>
                              <div><span style={{ color: C.muted }}>dMM: </span><span style={{ fontFamily: 'monospace' }}>{(pg.downMM||0).toFixed(3)}</span></div>
                              <div><span style={{ color: C.muted }}>cut: </span><span style={{ fontFamily: 'monospace' }}>+{pg.cutPosition}bp</span></div>
                            </div>
                            {pg.grnaStrand && (
                              <div style={{ marginTop: 6, fontSize: 10, color: pg.grnaStrand==='-'?'#BF5AF2':C.cyan }}>
                                {pg.grnaStrand==='-'?'(−) minus strand':'(+) plus strand'}
                              </div>
                            )}
                          </div>
                        );
                      })}
                    </div>
                    {selectedSample.multiGrnaSummary && (
                      <div style={{ marginTop: 10, fontSize: 11, color: '#BF5AF2', fontStyle: 'italic' }}>
                        {selectedSample.multiGrnaSummary}
                      </div>
                    )}
                  </div>
                )}

                {/* Indel details */}
                {selectedSample.indels?.length > 0 ? (
                  <div>
                    <div style={{ fontSize: 11, color: C.muted, letterSpacing: 1, marginBottom: 10, textTransform: 'uppercase' }}>Indel Details</div>
                    <table style={{ width: '100%', borderCollapse: 'collapse', fontSize: 12 }}>
                      <thead>
                        <tr>{['Position','Type','Size','Rel. to Cut Site','Frame Effect'].map(h => <th key={h} style={th}>{h}</th>)}</tr>
                      </thead>
                      <tbody>
                        {selectedSample.indels.map((ind, i) => (
                          <tr key={i}>
                            <td style={td}><code style={{ color: C.blue }}>{ind.position}</code></td>
                            <td style={td}><Chip color={ind.type === 'deletion' ? C.red : ind.type === 'insertion' ? C.cyan : C.gold}>{ind.type}</Chip></td>
                            <td style={td}><code>{ind.length} bp</code></td>
                            <td style={td}><span style={{ color: Math.abs(ind.relToCutSite) <= 5 ? C.cyan : C.muted }}>{ind.relToCutSite > 0 ? '+' : ''}{ind.relToCutSite} bp</span></td>
                            <td style={td}><span style={{ color: ind.isFrameshift ? C.red : C.gold, fontWeight: 600 }}>{ind.isFrameshift ? '🔴 Frameshift' : '🟡 In-frame'}</span></td>
                          </tr>
                        ))}
                      </tbody>
                    </table>
                  </div>
                ) : (
                  <div style={{ background: 'rgba(0,212,170,0.07)', border: '1px solid rgba(0,212,170,0.2)', borderRadius: 10, padding: '12px 18px', color: C.cyan, fontSize: 13 }}>
                    ✓ No indels detected — wild-type at target site.
                  </div>
                )}
              </div>
            ) : (
              <div style={{ ...card, display: 'flex', flexDirection: 'column', alignItems: 'center', justifyContent: 'center', minHeight: 400, color: C.muted }}>
                <div style={{ fontSize: 48, marginBottom: 16 }}>🧬</div>
                <div style={{ fontSize: 16, fontWeight: 700, color: C.text, marginBottom: 8 }}>Select a sample</div>
                <div style={{ fontSize: 12, color: C.muted, textAlign: 'center', maxWidth: 260 }}>
                  Click any sample from the list on the left to view its full analysis details, metrics, and indel information.
                </div>
                {explorerFilter !== 'all' && (
                  <div style={{ marginTop: 16, fontSize: 11, color: C.cyan, background: 'rgba(0,212,170,0.08)', border: '1px solid rgba(0,212,170,0.2)', borderRadius: 8, padding: '8px 16px' }}>
                    Filtered by: <b>{explorerFilter}</b> — {results.filter(r => {
                      if (explorerFilter === 'edited')  return r.zygosityClass > 0;
                      if (explorerFilter === 'homo')    return r.zygosityClass === 3;
                      if (explorerFilter === 'het')     return r.zygosityClass === 2;
                      if (explorerFilter === 'wt')      return r.zygosityClass === 0;
                      if (explorerFilter === 'fs')      return r.isFrameshift;
                      if (explorerFilter === 'flagged') return !!r.flag;
                      return true;
                    }).length} samples shown
                  </div>
                )}
              </div>
            )}
          </div>
        )}

        {/* ═══════════ TAB 3 — AI INSIGHTS ═══════════ */}
        {tab === 3 && results && (() => {
          const homo   = results.filter(r => r.zygosityClass === 3);
          const het    = results.filter(r => r.zygosityClass === 2);
          const wt     = results.filter(r => r.zygosityClass === 0);
          const low    = results.filter(r => r.zygosityClass === 1);
          const flagged= results.filter(r => r.flag);
          const ready   = homo.concat(het).filter(r => !r.flag); // eslint-disable-line no-unused-vars
          return (
            <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: 22 }}>
              {/* Recommendations — data-driven */}
              <div style={card}>
                <SectionTitle>🤖 AI Recommendations</SectionTitle>
                {[
                  homo.length > 0 && {
                    priority: 'high', emoji: '🏆',
                    text: `${homo.length} homozygous edited plant${homo.length>1?'s':''} confirmed: ${homo.map(r=>r.id).join(', ')}. Advance to T2 immediately — self-pollinate to fix the edit.`
                  },
                  het.filter(r=>!r.flag).length > 0 && {
                    priority: 'high', emoji: '🌾',
                    text: `${het.filter(r=>!r.flag).length} heterozygous plant${het.filter(r=>!r.flag).length>1?'s':''} ready for T2 advancement: ${het.filter(r=>!r.flag).map(r=>r.id).join(', ')}. Self-pollinate → screen ~8 T2 plants to recover homozygotes.`
                  },
                  flagged.length > 0 && {
                    priority: 'medium', emoji: '⚠️',
                    text: `${flagged.length} plant${flagged.length>1?'s':''} need re-sequencing before advancing: ${flagged.map(r=>r.id).join(', ')}. Reason: ${[...new Set(flagged.map(r=>r.flag))].join('; ')}.`
                  },
                  {
                    priority: 'high', emoji: '🔬',
                    text: 'Screen all edited lines for T-DNA (Cas9/Hygromycin PCR). Transgene-free edited plants (~37% of T1) are your final target material.'
                  },
                  wt.length > 0 && {
                    priority: 'medium', emoji: '🌱',
                    text: `${wt.length} wild-type segregant${wt.length>1?'s':''} identified: ${wt.map(r=>r.id).join(', ')}. Use as null-segregant controls for phenotyping.`
                  },
                  low.length > 0 && {
                    priority: 'medium', emoji: '🔄',
                    text: `${low.length} plant${low.length>1?'s':''} with low/mosaic edits: ${low.map(r=>r.id).join(', ')}. Confirm by clonal sequencing (TA-clone ≥8 colonies) to identify allele composition.`
                  },
                  {
                    priority: 'low', emoji: '📋',
                    text: `Screen ~${Math.max(16, results.length)} T1 plants per line to find ≥3 homozygous edited plants with 95% probability.`
                  },
                ].filter(Boolean).map((rec, i) => (
                  <div key={i} style={{ background: rec.priority==='high' ? 'rgba(255,77,109,0.06)' : rec.priority==='medium' ? 'rgba(240,192,64,0.06)' : 'rgba(0,212,170,0.06)', border: `1px solid ${rec.priority==='high' ? 'rgba(255,77,109,0.25)' : rec.priority==='medium' ? 'rgba(240,192,64,0.25)' : 'rgba(0,212,170,0.25)'}`, borderRadius: 12, padding: '12px 16px', display: 'flex', gap: 12, marginBottom: 10 }}>
                    <div style={{ fontSize: 18, flexShrink: 0 }}>{rec.emoji}</div>
                    <div>
                      <div style={{ fontSize: 9, letterSpacing: 1, color: rec.priority==='high' ? C.red : rec.priority==='medium' ? C.gold : C.cyan, marginBottom: 3 }}>{rec.priority.toUpperCase()} PRIORITY</div>
                      <div style={{ fontSize: 12, color: C.text, lineHeight: 1.6 }}>{rec.text}</div>
                    </div>
                  </div>
                ))}
              </div>

              {/* Plant advancement strategy */}
              <div style={card}>
                <SectionTitle>🌾 Plant Advancement Strategy</SectionTitle>
                {[
                  { icon: '🏆', title: 'Priority Advance', color: C.cyan,   plants: homo, desc: 'Confirmed homozygous — advance to T2 now' },
                  { icon: '🌿', title: 'Advance (het)',     color: C.purple, plants: het.filter(r=>!r.flag), desc: 'Heterozygous — self-pollinate to recover homo T2' },
                  { icon: '⚠️', title: 'Verify First',     color: C.gold,   plants: flagged, desc: 'Low confidence — re-sequence before advancing' },
                  { icon: '🌱', title: 'Use as Controls',  color: C.dim,    plants: wt, desc: 'Wild-type — null-segregant controls' },
                ].map((g, i) => (
                  <div key={i} style={{ marginBottom: 12, padding: '12px 14px', background: `${g.color}0d`, border: `1px solid ${g.color}30`, borderRadius: 10 }}>
                    <div style={{ display: 'flex', gap: 10, alignItems: 'center', marginBottom: g.plants.length ? 6 : 0 }}>
                      <span style={{ fontSize: 16 }}>{g.icon}</span>
                      <div style={{ flex: 1 }}>
                        <div style={{ color: g.color, fontWeight: 700, fontSize: 12 }}>{g.title} ({g.plants.length} plant{g.plants.length!==1?'s':''})</div>
                        <div style={{ color: C.muted, fontSize: 10 }}>{g.desc}</div>
                      </div>
                    </div>
                    {g.plants.length > 0 && (
                      <div style={{ background: 'rgba(0,0,0,0.3)', borderRadius: 6, padding: '5px 10px', fontSize: 11, color: g.color, fontFamily: 'monospace' }}>
                        {g.plants.map(p => p.id).join(', ')}
                      </div>
                    )}
                  </div>
                ))}

                {/* T1 segregation guide */}
                <div style={{ marginTop: 18 }}>
                  <SectionTitle>📊 Expected T1 Segregation</SectionTitle>
                  {[
                    { label: 'Homozygous edited', color: C.cyan,   ratio: '~25%', desc: 'Both alleles carry the edit — your target' },
                    { label: 'Heterozygous',       color: C.purple, ratio: '~50%', desc: 'One allele edited — self to get homozygotes' },
                    { label: 'Wild-type',           color: C.muted,  ratio: '~25%', desc: 'No editing — use as null-segregant control' },
                    { label: 'Transgene-free edited',color: C.gold,  ratio: '~37%', desc: 'Edited plants without T-DNA — cleanest lines' },
                  ].map((p, i) => (
                    <div key={i} style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', padding: '7px 0', borderBottom: i < 3 ? `1px solid ${C.border}` : 'none' }}>
                      <div>
                        <div style={{ color: p.color, fontSize: 12, fontWeight: 600 }}>{p.label}</div>
                        <div style={{ color: C.muted, fontSize: 10 }}>{p.desc}</div>
                      </div>
                      <span style={{ color: p.color, fontWeight: 700, fontSize: 16, marginLeft: 16 }}>{p.ratio}</span>
                    </div>
                  ))}
                  <div style={{ marginTop: 10, fontSize: 11, color: C.cyan }}>💡 Screen ~16 T1 plants per line to find ≥3 homozygous edited plants (95% probability)</div>
                </div>
              </div>
            </div>
          );
        })()}

        {/* ═══════════ TAB 4 — ML MODEL ═══════════ */}
        {tab === 4 && (
          <MLQualityPanel
            results={results || []}
            selected={selectedSample}
            onSelect={s => setSelectedSample(s)}
          />
        )}

        {/* ═══════════ TAB 5 — REPORT ═══════════ */}
        {tab === 5 && results && (() => {
          const reportCutPos = cfg.ref && grnaList[0]?.sequence
            ? deriveCutPosition(cfg.ref, grnaList[0].sequence).cutPosition
            : 137;
          return (
          <div>
            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 20, flexWrap: 'wrap', gap: 12 }}>
              <div>
                <div style={{ fontSize: 18, fontWeight: 700 }}>📄 Publication-Ready Report</div>
                <div style={{ fontSize: 12, color: C.muted, marginTop: 4 }}>Export data for records and manuscript preparation</div>
              </div>
              <div style={{ display: 'flex', gap: 10, flexWrap: 'wrap' }}>
                <button style={btnSec} onClick={() => downloadExcel(results, cfg, stats, wtControlResult, grnaList)}>⬇ Download Excel (.xls)</button>
                <button style={btn} onClick={() => generatePDF(results, cfg, stats, wtControlResult, grnaList)}>⬇ Download PDF Report</button>
              </div>
            </div>

            <div style={{ display: 'grid', gridTemplateColumns: '200px 1fr', gap: 22 }}>
              {/* Summary sidebar */}
              <div style={{ display: 'flex', flexDirection: 'column', gap: 10 }}>
                {[
                  { icon: '🧬', l: 'Gene',       v: cfg.geneName || '—' },
                  { icon: '🌾', l: 'Organism',   v: cfg.organism },
                  { icon: '🌱', l: 'Generation', v: cfg.generation },
                  { icon: '🔬', l: 'Samples',    v: stats.total },
                  { icon: '✂️', l: 'Efficiency', v: `${stats.eff}%`, hi: true },
                  { icon: '📐', l: 'Frameshift', v: `${stats.fs} (${stats.total>0?Math.round(stats.fs/stats.total*100):0}%)` },
                  { icon: '★',  l: 'Homozygous', v: stats.homo, hi: true },
                  { icon: '◑',  l: 'Heterozygous',v: stats.het },
                ].map((s, i) => (
                  <div key={i} style={{ background: s.hi ? 'rgba(0,212,170,0.08)' : C.card2, border: `1px solid ${s.hi ? 'rgba(0,212,170,0.4)' : C.border}`, borderRadius: 10, padding: '10px 14px', display: 'flex', alignItems: 'center', gap: 10 }}>
                    <span style={{ fontSize: 16 }}>{s.icon}</span>
                    <div>
                      <div style={{ fontSize: 9, color: C.muted, letterSpacing: 1, textTransform: 'uppercase' }}>{s.l}</div>
                      <div style={{ fontSize: 13, fontWeight: 700, color: s.hi ? C.cyan : C.text }}>{s.v}</div>
                    </div>
                  </div>
                ))}
                {wtControlResult && (
                  <div style={{ background: 'rgba(0,212,170,0.06)', border: '1px solid rgba(0,212,170,0.3)', borderRadius: 10, padding: '10px 14px' }}>
                    <div style={{ fontSize: 9, color: C.muted, letterSpacing: 1 }}>WT CONTROL</div>
                    <div style={{ fontSize: 12, fontWeight: 700, color: C.cyan }}>{wtControlResult.id}</div>
                  </div>
                )}
              </div>

              {/* Report preview */}
              <div style={{ background: C.card2, border: `1px solid ${C.border}`, borderRadius: 14, padding: 22, fontSize: 11, lineHeight: 1.9, color: C.muted, fontFamily: 'monospace', maxHeight: 680, overflowY: 'auto' }}>
                <div style={{ color: C.cyan, fontWeight: 700, fontSize: 13, marginBottom: 4 }}>PlantaCRISPAI — CRISPR Analysis Report</div>
                <div style={{ color: C.blue, marginBottom: 16 }}>
                  Gene: {cfg.geneName||'—'} | {cfg.organism} | {cfg.generation}
                  {grnaList.length === 1
                    ? ` | gRNA: ${grnaList[0].sequence || cfg.grna || '—'}`
                    : ` | ${grnaList.length} gRNAs: ${grnaList.map(g=>g.label||g.sequence.slice(0,8)+'…').join(', ')}`
                  }
                </div>
                {wtControlResult && <div style={{ color: C.cyan, marginBottom: 8 }}>WT Control: {wtControlResult.id}</div>}
                <div style={{ color: C.text, marginBottom: 4 }}>Editing Efficiency: <span style={{ color: C.cyan, fontWeight: 700 }}>{stats.eff}%</span></div>
                <div style={{ color: C.text, marginBottom: 4 }}>Edited: {stats.edited}/{stats.total} | Homo: {stats.homo} | Het: {stats.het} | Low: {stats.low} | WT: {stats.wt}</div>
                <div style={{ color: C.text, marginBottom: 16 }}>Frameshift: {stats.fs} | Flagged: {stats.flagged}</div>
                <div style={{ borderTop: `1px solid ${C.border}`, paddingTop: 16 }}>
                  {results.map((r, i) => (
                    <div key={i} style={{ marginBottom: 16, paddingBottom: 14, borderBottom: `1px solid rgba(30,58,95,0.4)` }}>
                      <div style={{ color: isStarSample(r) ? C.cyan : C.blue, fontWeight: 700, fontSize: 12 }}>
                        {isStarSample(r) ? '★ ' : ''}{r.id}
                      </div>
                      <div>Edit: {r.editType}{r.indelSize>0?` (${r.indelSize} bp)`:''} | Direction: {r.lenDiff>0?`+${r.lenDiff} bp insertion`:r.lenDiff<0?`${r.lenDiff} bp deletion`:'none'}</div>
                      <div>Zygosity: <span style={{ color: ZYG[r.zygosityClass]?.color }}>{ZYG[r.zygosityClass]?.label}</span></div>
                      <div>ICE: {r.iceScore}% | R²: {(r.r2||0).toFixed(3)} (up: {(r.r2up||0).toFixed(3)}, dn: {(r.r2dn||0).toFixed(3)}) | Identity: {r.sequenceIdentity?(r.sequenceIdentity*100).toFixed(1)+'%':'—'}</div>
                      <div>ML Score: <span style={{ color: C.cyan }}>{r.mlScore}%</span> | Frameshift: {r.isFrameshift?'Yes':'No'}</div>
                      <div>Protein: {r.proteinEffect}</div>
                      {r.flag && <div style={{ color: C.gold }}>⚠ FLAG: {r.flag}</div>}
                    </div>
                  ))}
                </div>
                <div style={{ borderTop: `1px solid ${C.border}`, paddingTop: 14, color: C.dim, fontSize: 10, lineHeight: 1.7 }}>
                  METHODS: Sanger sequencing reads aligned to reference amplicon using Smith-Waterman local alignment. Indels identified from pairwise gap patterns in aligned strings. Zygosity classified by three independent signals: (1) gap-aware downstream R² ≥0.97 + downstream mismatch &lt;3% = Homozygous (single-allele fingerprint); (2) downstream non-gap mismatch rate ≥8% = Heterozygous (double-peak signal from two alleles); (3) indel present with mismatch above WT noise floor (0.034) = Low/Mosaic. R² computed on non-gap base pairs only upstream and downstream of cut site. ICE-equivalent score: Homo = 85 + r2dn×14%; Het/Low = mismatch rate × 270%. Validated on 19 OsAAP3 T1 samples confirmed by ICE+TIDE (17/19 = 89.5% accuracy). PlantaCRISPAI v2.0 — analysisEngine v5.0.
                </div>
              </div>
            </div>{/* end 200px/1fr grid */}

            {/* ── Full-width visualisation panels ── */}
            {/* AllelePlot first — then IndelChart below */}
            <AllelePlot results={results} cfg={{ ref: cfg.ref, cutPos: reportCutPos, grnaSeq: grnaList[0]?.sequence || cfg.grna }} />
            <IndelChart results={results} cfg={{ ref: cfg.ref, cutPos: reportCutPos }} />
          </div>
          );
        })()}

        {/* ═══════════ TAB 6 — CHROMATOGRAM ═══════════ */}
        {tab === 6 && (
          <div style={{ display: 'flex', flexDirection: 'column', gap: 20 }}>
            <div style={{ background: 'rgba(77,159,255,0.06)', border: `1px solid rgba(77,159,255,0.25)`, borderRadius: 12, padding: '12px 20px', fontSize: 12, color: C.blue }}>
              ℹ️ Upload .ab1 files in Setup to view real fluorescence traces. Synthetic traces are shown below for samples loaded from text sequences.
            </div>
            <div style={{ ...card, padding: 0, overflow: 'hidden' }}>
              <ChromatogramViewer
                sample={selectedSample || (results && results[0]) || null}
                wtSample={wtControlResult}
                allSamples={results || []}
                onSelect={s => setSelectedSample(s)}
                cutPosition={cfg.ref && grnaList[0]?.sequence ? deriveCutPosition(cfg.ref, grnaList[0].sequence).cutPosition : 137}
                grnaSeq={grnaList[0]?.sequence || cfg.grna || 'CGGCGTGCTCTCTCTCGCAT'}
                grnaList={grnaList}
              />
            </div>
          </div>
        )}

        {/* ═══════════ TAB 7 — CAPS DESIGNER ═══════════ */}
        {tab === 7 && (
          <CAPSDesignerPanel results={results} cfg={cfg} refSequence={cfg.ref} />
        )}

        {/* Empty state */}
        {tab > 0 && tab < 7 && !results && (
          <div style={{ textAlign: 'center', padding: '80px 32px' }}>
            <div style={{ fontSize: 64, marginBottom: 18 }}>🧬</div>
            <div style={{ fontSize: 20, fontWeight: 700, marginBottom: 10 }}>No Analysis Run Yet</div>
            <div style={{ color: C.muted, marginBottom: 24 }}>Go to Setup, configure your experiment, and click Run Analysis.</div>
            <button style={btn} onClick={() => navigateTo(0)}>← Go to Setup</button>
          </div>
        )}

      </div>
    </div>
  );
}

export default AnalysisApp;
