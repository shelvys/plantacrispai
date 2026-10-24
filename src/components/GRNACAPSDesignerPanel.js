// ============================================================
// FILE: src/components/GRNACAPSDesignerPanel.js
// PURPOSE: ML-powered gRNA + CAPS designer UI (Mode C upgrade)
//
// Shows per-gRNA:
//   • Ensemble ML score with confidence interval
//   • 4-model breakdown bar chart
//   • Thermodynamic profile (ΔG, Tm, seed stability)
//   • Position-specific feature heatmap
//   • CAPS compatibility with best enzyme
//   • Pareto frontier visualization
//   • Off-target risk assessment
//   • Exportable methods text for publication
// ============================================================

import React, { useState } from 'react';
import { designGRNAsWithCAPS } from '../engine/grnaDesigner';

// ── Theme ──────────────────────────────────────────────────
const T = {
  bg:      '#060D1F', card: '#0D1B35', border: '#1a3060',
  cyan:    '#00D4AA', blue: '#4D9FFF', purple: '#B060FF',
  gold:    '#F0C040', red:  '#FF4D6D', green:  '#00D4AA',
  orange:  '#FF8C42', text: '#E8F0FE', muted:  '#8899BB',
  mutedDk: '#3A4A6A',
};

const card  = { background: T.card, border: `1px solid ${T.border}`, borderRadius: 14, padding: 20 };
const inp   = { background: T.bg, border: `1px solid ${T.border}`, borderRadius: 8, color: T.text, padding: '9px 13px', fontSize: 13, fontFamily: 'Calibri, sans-serif', outline: 'none', boxSizing: 'border-box', width: '100%' };
const lbl   = { fontSize: 10, color: T.muted, letterSpacing: 1, textTransform: 'uppercase', marginBottom: 5, display: 'block' };
const btn   = { background: `linear-gradient(135deg, ${T.cyan}, #00a882)`, color: '#000', border: 'none', borderRadius: 10, padding: '10px 24px', fontSize: 13, fontWeight: 700, cursor: 'pointer' };
const btnSm = { background: 'transparent', color: T.cyan, border: `1px solid ${T.cyan}55`, borderRadius: 8, padding: '5px 14px', fontSize: 11, fontWeight: 600, cursor: 'pointer' };

// ── Reusable mini components ───────────────────────────────

function ScoreRing({ score, size = 56, color }) {
  const c     = color || (score >= 75 ? T.cyan : score >= 55 ? T.blue : score >= 35 ? T.gold : T.red);
  const r     = (size - 8) / 2;
  const circ  = 2 * Math.PI * r;
  const dash  = (score / 100) * circ;
  return (
    <svg width={size} height={size} style={{ display: 'block' }}>
      <circle cx={size/2} cy={size/2} r={r} fill="none"
        stroke={`${c}22`} strokeWidth={5} />
      <circle cx={size/2} cy={size/2} r={r} fill="none"
        stroke={c} strokeWidth={5}
        strokeDasharray={`${dash} ${circ}`}
        strokeLinecap="round"
        transform={`rotate(-90 ${size/2} ${size/2})`} />
      <text x={size/2} y={size/2 + 1} textAnchor="middle"
        dominantBaseline="middle"
        fill={c} fontSize={size > 50 ? 13 : 10} fontWeight="800"
        fontFamily="monospace">{score}%</text>
    </svg>
  );
}

function BarBreakdown({ breakdown, size = 'md' }) {
  const bars = [
    { label: 'Position', key: 'positionWeight', color: T.cyan   },
    { label: 'Thermo',   key: 'thermodynamic',  color: T.blue   },
    { label: 'Structure',key: 'structural',      color: T.purple },
    { label: 'Plant',    key: 'plantSpecific',   color: T.gold   },
  ];
  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: 5 }}>
      {bars.map(({ label, key, color }) => {
        const val = breakdown?.[key] || 0;
        return (
          <div key={key}>
            <div style={{ display: 'flex', justifyContent: 'space-between',
                          marginBottom: 2 }}>
              <span style={{ fontSize: 10, color: T.muted }}>{label}</span>
              <span style={{ fontSize: 10, color, fontWeight: 700 }}>{val}%</span>
            </div>
            <div style={{ background: `${color}18`, borderRadius: 3, height: 5 }}>
              <div style={{ width: `${val}%`, height: '100%',
                            background: color, borderRadius: 3,
                            transition: 'width 0.4s ease' }} />
            </div>
          </div>
        );
      })}
    </div>
  );
}

function ThermoCard({ thermo }) {
  if (!thermo) return null;
  const dG       = thermo.fullDuplex?.dG;
  const seedDG   = thermo.seedRegion?.dG;
  const tm       = thermo.meltingTemp;
  const selfFold = thermo.selfFolding;
  const interp   = thermo.interpretation;

  return (
    <div style={{ background: T.bg, borderRadius: 10, padding: '12px 14px' }}>
      <div style={{ fontSize: 10, color: T.muted, letterSpacing: 1,
                    textTransform: 'uppercase', marginBottom: 10 }}>
        Thermodynamic Profile
      </div>
      <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: 8 }}>
        {[
          { label: 'Full ΔG', val: `${dG} kcal/mol`,  color: dG >= -28 && dG <= -18 ? T.cyan : T.gold },
          { label: 'Seed ΔG', val: `${seedDG} kcal/mol`, color: seedDG <= -10 ? T.cyan : T.gold },
          { label: 'Tm',      val: `${tm}°C`,           color: T.blue },
          { label: 'Self-fold', val: selfFold?.risk || '—', color: selfFold?.color || T.muted },
        ].map(({ label, val, color }) => (
          <div key={label} style={{ background: `${color}0A`,
            border: `1px solid ${color}33`, borderRadius: 8, padding: '6px 10px' }}>
            <div style={{ fontSize: 9, color: T.muted, letterSpacing: 1 }}>{label}</div>
            <div style={{ fontSize: 12, fontWeight: 700, color, fontFamily: 'monospace' }}>{val}</div>
          </div>
        ))}
      </div>
      {interp && (
        <div style={{ marginTop: 8, fontSize: 11, color: interp.color }}>
          → {interp.text}
        </div>
      )}
    </div>
  );
}

function PositionHeatmap({ grnaSeq }) {
  if (!grnaSeq) return null;
  const seq = grnaSeq.toUpperCase();

  // Base color map
  const baseColors = { G: T.cyan, C: T.blue, A: T.gold, T: T.red, N: T.muted };

  // Importance weights per position (from Doench 2016)
  const posImportance = [
    0.3,0.3,0.4,0.5,0.4,0.5,0.4,0.6,
    0.5,0.6,0.6,0.8,0.8,0.7,0.7,0.8,
    0.8,0.9,0.9,1.0  // position 20 is most important
  ];

  return (
    <div>
      <div style={{ fontSize: 10, color: T.muted, letterSpacing: 1,
                    textTransform: 'uppercase', marginBottom: 8 }}>
        Position Map (1 → 20, PAM-proximal at right)
      </div>
      <div style={{ display: 'flex', gap: 2 }}>
        {seq.split('').map((base, i) => {
          const color   = baseColors[base] || T.muted;
          const opacity = 0.3 + posImportance[i] * 0.7;
          const isKey   = [0,3,5,11,12,17,18,19].includes(i);
          return (
            <div key={i} title={`Position ${i+1}: ${base}`}
              style={{ flex: 1, display: 'flex', flexDirection: 'column',
                       alignItems: 'center', gap: 2 }}>
              {/* Base letter */}
              <div style={{
                width: '100%', textAlign: 'center',
                fontSize: 11, fontWeight: isKey ? 800 : 500,
                fontFamily: 'monospace', color,
                opacity,
              }}>{base}</div>
              {/* Importance bar */}
              <div style={{
                width: '80%', height: Math.round(posImportance[i] * 24),
                background: color, opacity, borderRadius: 2,
                minHeight: 4,
              }} />
              {/* Position number for key positions */}
              {isKey && (
                <div style={{ fontSize: 7, color: T.muted }}>{i+1}</div>
              )}
            </div>
          );
        })}
      </div>
      {/* Legend */}
      <div style={{ display: 'flex', gap: 12, marginTop: 8, flexWrap: 'wrap' }}>
        {[['G', T.cyan], ['C', T.blue], ['A', T.gold], ['T', T.red]].map(([b, c]) => (
          <div key={b} style={{ display: 'flex', alignItems: 'center', gap: 4 }}>
            <div style={{ width: 8, height: 8, background: c, borderRadius: 2 }} />
            <span style={{ fontSize: 10, color: T.muted }}>{b}</span>
          </div>
        ))}
        <span style={{ fontSize: 10, color: T.muted, marginLeft: 8 }}>
          Bar height = position importance
        </span>
      </div>
    </div>
  );
}

function ParetoChart({ candidates }) {
  if (!candidates || candidates.length === 0) return null;

  const W = 280, H = 200;
  const PAD = 32;

  // Normalize scores to plot area
  const toX = v => PAD + (v / 100) * (W - PAD * 2);
  const toY = v => H - PAD - (v / 100) * (H - PAD * 2);

  const pareto  = candidates.filter(c => c.isPareto);
  const nonPar  = candidates.filter(c => !c.isPareto);

  return (
    <div>
      <div style={{ fontSize: 10, color: T.muted, letterSpacing: 1,
                    textTransform: 'uppercase', marginBottom: 8 }}>
        Pareto Frontier — Efficiency vs CAPS
      </div>
      <svg width={W} height={H} style={{ display: 'block', overflow: 'visible' }}>
        {/* Grid */}
        {[0, 25, 50, 75, 100].map(v => (
          <g key={v}>
            <line x1={toX(v)} y1={PAD} x2={toX(v)} y2={H - PAD}
              stroke={T.mutedDk} strokeWidth={0.5} strokeDasharray="3,4" />
            <line x1={PAD} y1={toY(v)} x2={W - PAD} y2={toY(v)}
              stroke={T.mutedDk} strokeWidth={0.5} strokeDasharray="3,4" />
            <text x={toX(v)} y={H - PAD + 12} textAnchor="middle"
              fill={T.muted} fontSize={8}>{v}</text>
            <text x={PAD - 4} y={toY(v) + 3} textAnchor="end"
              fill={T.muted} fontSize={8}>{v}</text>
          </g>
        ))}

        {/* Axis labels */}
        <text x={W / 2} y={H - 4} textAnchor="middle"
          fill={T.muted} fontSize={9}>Efficiency Score</text>
        <text x={10} y={H / 2} textAnchor="middle"
          fill={T.muted} fontSize={9}
          transform={`rotate(-90, 10, ${H/2})`}>CAPS Score</text>

        {/* Non-pareto points */}
        {nonPar.map((c, i) => (
          <circle key={i}
            cx={toX(c.scores.efficiency)}
            cy={toY(c.scores.capsCompat)}
            r={4} fill={T.mutedDk} opacity={0.7} />
        ))}

        {/* Pareto front line */}
        {pareto.length > 1 && (() => {
          const sorted = [...pareto].sort((a,b) =>
            a.scores.efficiency - b.scores.efficiency);
          const pts = sorted.map(c =>
            `${toX(c.scores.efficiency)},${toY(c.scores.capsCompat)}`
          ).join(' ');
          return <polyline points={pts} fill="none"
            stroke={`${T.cyan}55`} strokeWidth={1.5}
            strokeDasharray="4,3" />;
        })()}

        {/* Pareto points */}
        {pareto.map((c, i) => (
          <g key={i}>
            <circle
              cx={toX(c.scores.efficiency)}
              cy={toY(c.scores.capsCompat)}
              r={6} fill={T.cyan} opacity={0.9} />
            <text
              x={toX(c.scores.efficiency) + 8}
              y={toY(c.scores.capsCompat) + 3}
              fill={T.cyan} fontSize={7}>
              #{candidates.indexOf(c) + 1}
            </text>
          </g>
        ))}

        {/* Top candidate star */}
        {candidates[0] && (
          <text
            x={toX(candidates[0].scores.efficiency)}
            y={toY(candidates[0].scores.capsCompat) - 9}
            textAnchor="middle" fontSize={12}>⭐</text>
        )}
      </svg>
      <div style={{ fontSize: 11, color: T.muted, marginTop: 4 }}>
        <span style={{ color: T.cyan }}>● Pareto-optimal</span>
        <span style={{ color: T.mutedDk, marginLeft: 12 }}>● Dominated</span>
      </div>
    </div>
  );
}

// ── gRNA Row (expanded detail) ─────────────────────────────
function GRNARow({ cand, rank, isExpanded, onToggle }) {
  const s = cand.scores;
  const c = s.combined >= 75 ? T.cyan
          : s.combined >= 55 ? T.blue
          : s.combined >= 35 ? T.gold : T.red;

  return (
    <div style={{ marginBottom: 10 }}>
      {/* Collapsed row */}
      <div
        onClick={onToggle}
        style={{
          background: isExpanded ? `${c}0F` : '#060D1F',
          border: `1px solid ${isExpanded ? c : T.mutedDk}`,
          borderRadius: isExpanded ? '10px 10px 0 0' : 10,
          padding: '10px 14px',
          cursor: 'pointer',
          display: 'grid',
          gridTemplateColumns: '36px 1fr 60px 60px 60px 100px 80px 40px',
          alignItems: 'center',
          gap: 10,
          transition: 'all 0.15s',
        }}
      >
        {/* Rank */}
        <div style={{ textAlign: 'center' }}>
          <span style={{ fontSize: 13, fontWeight: 800, color: rank === 1 ? T.cyan : T.muted }}>
            #{rank}
          </span>
          {cand.isPareto && (
            <div style={{ fontSize: 8, color: T.cyan, marginTop: 1 }}>PARETO</div>
          )}
        </div>

        {/* Sequence */}
        <div>
          <code style={{ fontSize: 12, color: T.blue, letterSpacing: '1.2px',
                          fontFamily: 'monospace', wordBreak: 'break-all' }}>
            {cand.sequence}
          </code>
          <span style={{ marginLeft: 6, fontSize: 10, color: T.purple,
                          fontFamily: 'monospace' }}>{cand.pam}</span>
          <div style={{ fontSize: 10, color: T.muted, marginTop: 2 }}>
            {cand.strand === '+' ? '→' : '←'} pos {cand.position} · GC {cand.gcContent}%
          </div>
        </div>

        {/* Efficiency */}
        <div style={{ textAlign: 'center' }}>
          <ScoreRing score={s.efficiency} size={44}
            color={s.efficiency >= 70 ? T.cyan : s.efficiency >= 50 ? T.blue : T.gold} />
          <div style={{ fontSize: 8, color: T.muted, marginTop: 2 }}>EFFIC.</div>
        </div>

        {/* CAPS */}
        <div style={{ textAlign: 'center' }}>
          <ScoreRing score={s.capsCompat} size={44}
            color={s.capsCompat >= 70 ? T.cyan : s.capsCompat >= 50 ? T.blue : T.gold} />
          <div style={{ fontSize: 8, color: T.muted, marginTop: 2 }}>CAPS</div>
        </div>

        {/* Combined */}
        <div style={{ textAlign: 'center' }}>
          <ScoreRing score={s.combined} size={44} color={c} />
          <div style={{ fontSize: 8, color: T.muted, marginTop: 2 }}>TOTAL</div>
        </div>

        {/* Best CAPS enzyme */}
        <div>
          {cand.capsData?.hasCAPSCompatibility ? (
            <div>
              <span style={{ background: `${T.cyan}18`, border: `1px solid ${T.cyan}44`,
                borderRadius: 20, padding: '2px 10px', color: T.cyan,
                fontSize: 11, fontWeight: 700 }}>
                {cand.capsData.bestEnzyme}
              </span>
              <div style={{ fontSize: 9, color: T.muted, marginTop: 3 }}>
                {cand.capsData.capsCategory}
              </div>
            </div>
          ) : (
            <span style={{ fontSize: 10, color: T.muted }}>dCAPS needed</span>
          )}
        </div>

        {/* Off-target */}
        <div>
          <span style={{ fontSize: 11, fontWeight: 700,
                          color: cand.offTarget?.color || T.muted }}>
            {cand.offTarget?.level || '—'}
          </span>
          <div style={{ fontSize: 9, color: T.muted }}>off-tgt</div>
        </div>

        {/* Expand arrow */}
        <div style={{ textAlign: 'center', fontSize: 14, color: T.muted,
                      transform: isExpanded ? 'rotate(180deg)' : 'none',
                      transition: 'transform 0.2s' }}>▼</div>
      </div>

      {/* Expanded detail */}
      {isExpanded && (
        <div style={{
          background: `${c}06`,
          border: `1px solid ${c}`,
          borderTop: 'none',
          borderRadius: '0 0 10px 10px',
          padding: 16,
        }}>
          <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr 1fr', gap: 14 }}>

            {/* Model breakdown */}
            <div style={{ background: T.bg, borderRadius: 10, padding: 14 }}>
              <div style={{ fontSize: 10, color: T.muted, letterSpacing: 1,
                            textTransform: 'uppercase', marginBottom: 10 }}>
                ML Model Breakdown
              </div>
              <BarBreakdown breakdown={s.breakdown} />
              {s.ci && (
                <div style={{ marginTop: 10, fontSize: 11, color: T.muted }}>
                  95% CI: {s.ci.low}% – {s.ci.high}%
                </div>
              )}
            </div>

            {/* Thermodynamics */}
            <ThermoCard thermo={cand.mlResult?.thermo} />

            {/* Position heatmap */}
            <div style={{ background: T.bg, borderRadius: 10, padding: 14 }}>
              <PositionHeatmap grnaSeq={cand.sequence} />
            </div>
          </div>

          {/* Recommendations */}
          {cand.mlResult?.recommendation?.length > 0 && (
            <div style={{ marginTop: 14 }}>
              {cand.mlResult.recommendation.map((r, i) => (
                <div key={i} style={{ fontSize: 12, color: T.muted,
                                      marginBottom: 4, lineHeight: 1.5 }}>
                  {r}
                </div>
              ))}
            </div>
          )}

          {/* Recommendation badge */}
          <div style={{ marginTop: 12, display: 'inline-flex',
                        alignItems: 'center', gap: 8,
                        background: `${c}18`, border: `1px solid ${c}44`,
                        borderRadius: 8, padding: '6px 14px' }}>
            <span style={{ fontSize: 12, color: c, fontWeight: 700 }}>
              {cand.recommendation}
            </span>
          </div>
        </div>
      )}
    </div>
  );
}

// ══════════════════════════════════════════════════════════
// MAIN COMPONENT
// ══════════════════════════════════════════════════════════

export default function GRNACAPSDesignerPanel({ refSequence }) {
  const [customSeq, setCustomSeq]     = useState('');
  const [species, setSpecies]         = useState('rice');
  const [effWeight, setEffWeight]     = useState(50);
  const [capsWeight, setCapsWeight]   = useState(35);
  const [report, setReport]           = useState(null);
  const [running, setRunning]         = useState(false);
  const [expandedIdx, setExpandedIdx] = useState(0);
  const [showMethods, setShowMethods] = useState(false);
  const [activeView, setActiveView]   = useState('table'); // 'table'|'pareto'

  const seq = customSeq.trim() || refSequence || '';

  const runDesign = async () => {
    if (!seq || seq.length < 30) return;
    setRunning(true);
    setReport(null);
    setExpandedIdx(0);
    await new Promise(r => setTimeout(r, 200));
    try {
      const result = designGRNAsWithCAPS(seq, {
        maxResults:   20,
        effWeight:    effWeight / 100,
        capsWeight:   capsWeight / 100,
        offTgtWeight: (100 - effWeight - capsWeight) / 100,
        species,
        minEfficiency: 25,
      });
      setReport(result);
    } catch (e) {
      setReport({ success: false, error: e.message });
    }
    setRunning(false);
  };

  const offWeight = Math.max(0, 100 - effWeight - capsWeight);

  return (
    <div style={{ fontFamily: 'Calibri, sans-serif' }}>

      {/* Header */}
      <div style={{ ...card, marginBottom: 20,
        background: 'linear-gradient(135deg, rgba(0,212,170,0.06), rgba(176,96,255,0.06))' }}>
        <div style={{ display: 'flex', justifyContent: 'space-between',
                      alignItems: 'flex-start', flexWrap: 'wrap', gap: 10 }}>
          <div>
            <div style={{ fontSize: 20, fontWeight: 800, color: T.text }}>
              🌱 ML-Powered gRNA + CAPS Designer
            </div>
            <div style={{ fontSize: 12, color: T.muted, marginTop: 4, lineHeight: 1.6 }}>
              Feature-engineered ensemble ML · Nearest-neighbor thermodynamics ·
              Position-specific weights (Doench 2016) · Plant-species corrections ·
              Pareto-optimal ranking
            </div>
          </div>
          <div style={{ display: 'flex', gap: 8 }}>
            {['Position\nWeight', 'Thermo\ndynamics', 'Structural\nPenalty', 'Plant\nSpecific'].map((m, i) => (
              <div key={i} style={{ textAlign: 'center', background: `${[T.cyan,T.blue,T.purple,T.gold][i]}18`,
                border: `1px solid ${[T.cyan,T.blue,T.purple,T.gold][i]}44`,
                borderRadius: 8, padding: '6px 10px', minWidth: 60 }}>
                <div style={{ fontSize: 11, fontWeight: 700, color: [T.cyan,T.blue,T.purple,T.gold][i] }}>
                  {['25%','35%','20%','20%'][i]}
                </div>
                <div style={{ fontSize: 9, color: T.muted, whiteSpace: 'pre-line' }}>{m}</div>
              </div>
            ))}
          </div>
        </div>
      </div>

      {/* Input panel */}
      {!report && (
        <div style={{ ...card, marginBottom: 20 }}>
          <div style={{ display: 'grid', gridTemplateColumns: '1fr 240px', gap: 20 }}>

            {/* Sequence input */}
            <div>
              <label style={lbl}>Target Gene / Amplicon Sequence</label>
              <textarea
                style={{ ...inp, minHeight: 120, resize: 'vertical',
                          fontFamily: 'monospace', fontSize: 12, lineHeight: 1.7 }}
                value={customSeq || refSequence || ''}
                onChange={e => setCustomSeq(e.target.value.toUpperCase().replace(/\s/g,''))}
                placeholder="Paste your target gene or amplicon sequence here (minimum 30bp)&#10;&#10;Example: ATGGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC..."
              />
              {seq.length > 0 && (
                <div style={{ marginTop: 5, fontSize: 11, color: T.muted }}>
                  {seq.length} bp · GC: {Math.round((seq.match(/[GC]/g)||[]).length / seq.length * 100)}% ·
                  Est. candidates: ~{Math.floor(seq.length / 8)}
                </div>
              )}
            </div>

            {/* Settings */}
            <div style={{ display: 'flex', flexDirection: 'column', gap: 14 }}>

              <div>
                <label style={lbl}>Target Species</label>
                <select style={{ ...inp }}
                  value={species} onChange={e => setSpecies(e.target.value)}>
                  <option value="rice">🌾 Rice (Oryza sativa)</option>
                  <option value="wheat">🌿 Wheat (T. aestivum)</option>
                  <option value="arabidopsis">🌱 Arabidopsis thaliana</option>
                  <option value="maize">🌽 Maize (Zea mays)</option>
                  <option value="tomato">🍅 Tomato (S. lycopersicum)</option>
                  <option value="other">🌿 Other plant species</option>
                </select>
              </div>

              <div>
                <label style={lbl}>Scoring Weights</label>
                <div style={{ background: T.bg, borderRadius: 8, padding: 12 }}>
                  {[
                    { label: 'Efficiency', val: effWeight,  set: setEffWeight,  color: T.cyan },
                    { label: 'CAPS Compat.', val: capsWeight, set: setCapsWeight, color: T.blue },
                  ].map(({ label, val, set, color }) => (
                    <div key={label} style={{ marginBottom: 10 }}>
                      <div style={{ display: 'flex', justifyContent: 'space-between',
                                    marginBottom: 3 }}>
                        <span style={{ fontSize: 11, color: T.muted }}>{label}</span>
                        <span style={{ fontSize: 11, color, fontWeight: 700 }}>{val}%</span>
                      </div>
                      <input type="range" min={10} max={80} value={val}
                        onChange={e => set(Number(e.target.value))}
                        style={{ width: '100%', accentColor: color }} />
                    </div>
                  ))}
                  <div style={{ fontSize: 10, color: T.muted, textAlign: 'right' }}>
                    Off-target risk: {offWeight}%
                  </div>
                </div>
              </div>
            </div>
          </div>

          <div style={{ marginTop: 16 }}>
            <button style={{ ...btn, opacity: running || seq.length < 30 ? 0.6 : 1 }}
              onClick={runDesign} disabled={running || seq.length < 30}>
              {running
                ? '⚙ Running ML analysis...'
                : `🧬 Design gRNA + CAPS Candidates`}
            </button>
            {seq.length < 30 && seq.length > 0 && (
              <span style={{ marginLeft: 12, fontSize: 11, color: T.gold }}>
                ⚠ Need at least 30bp
              </span>
            )}
          </div>
        </div>
      )}

      {/* Results */}
      {report && report.success && (
        <div>
          {/* Summary row */}
          <div style={{ display: 'grid', gridTemplateColumns: 'repeat(5,1fr)',
                        gap: 12, marginBottom: 20 }}>
            {[
              { label: 'Candidates',    val: report.summary.scoredCandidates, color: T.cyan   },
              { label: 'Excellent (≥75%)', val: report.summary.excellentGRNAs, color: T.cyan },
              { label: 'With CAPS',     val: report.summary.withCAPSCompat,   color: T.blue   },
              { label: 'Pareto-optimal',val: report.summary.paretoOptimal,    color: T.purple },
              { label: 'Top Score',     val: `${report.summary.topGRNA?.scores.combined || 0}%`,
                color: T.gold },
            ].map((s,i) => (
              <div key={i} style={{ ...card, textAlign: 'center', padding: 14 }}>
                <div style={{ fontSize: 26, fontWeight: 800, color: s.color }}>{s.val}</div>
                <div style={{ fontSize: 9, color: T.muted, letterSpacing: 1,
                              textTransform: 'uppercase', marginTop: 4 }}>{s.label}</div>
              </div>
            ))}
          </div>

          {/* Species + methods info */}
          <div style={{ ...card, marginBottom: 16,
            background: `${T.cyan}07`, border: `1px solid ${T.cyan}33` }}>
            <div style={{ display: 'flex', justifyContent: 'space-between',
                          alignItems: 'center', flexWrap: 'wrap', gap: 10 }}>
              <div>
                <span style={{ fontSize: 12, color: T.cyan, fontWeight: 700 }}>
                  🌾 Species: {species.charAt(0).toUpperCase() + species.slice(1)}
                </span>
                <span style={{ fontSize: 12, color: T.muted, marginLeft: 16 }}>
                  Model: PWM + Thermodynamics + Structure + Plant correction
                </span>
                <span style={{ fontSize: 12, color: T.muted, marginLeft: 16 }}>
                  Weights: Eff {Math.round(effWeight)}% · CAPS {Math.round(capsWeight)}% · Off-tgt {offWeight}%
                </span>
              </div>
              <button style={btnSm} onClick={() => setShowMethods(v => !v)}>
                {showMethods ? 'Hide' : '📄 Methods Text'}
              </button>
            </div>

            {showMethods && (
              <div style={{ marginTop: 12, background: T.bg, borderRadius: 8,
                            padding: 14, fontSize: 11, color: T.muted,
                            lineHeight: 1.8, fontStyle: 'italic' }}>
                {report.methodology}
              </div>
            )}
          </div>

          {/* Recommendations */}
          {report.recommendations?.map((rec, i) => {
            const color = rec.priority === 'high' ? T.cyan
                        : rec.priority === 'medium' ? T.blue : T.muted;
            return (
              <div key={i} style={{ ...card, marginBottom: 10,
                background: `${color}08`, border: `1px solid ${color}33`,
                padding: '10px 16px' }}>
                <span style={{ fontSize: 13, color: T.text }}>
                  {rec.icon} {rec.text}
                </span>
              </div>
            );
          })}

          {/* View toggle */}
          <div style={{ display: 'flex', gap: 8, marginBottom: 16, marginTop: 16 }}>
            {['table', 'pareto'].map(v => (
              <button key={v}
                onClick={() => setActiveView(v)}
                style={{ ...btnSm,
                  background: activeView === v ? `${T.cyan}22` : 'transparent',
                  borderColor: activeView === v ? T.cyan : T.mutedDk,
                  color: activeView === v ? T.cyan : T.muted,
                }}>
                {v === 'table' ? '📋 Ranked Table' : '📊 Pareto Chart'}
              </button>
            ))}
            <button style={{ ...btnSm, marginLeft: 'auto' }}
              onClick={() => { setReport(null); setCustomSeq(''); }}>
              ← New Search
            </button>
          </div>

          {/* Pareto chart view */}
          {activeView === 'pareto' && (
            <div style={{ ...card, marginBottom: 16 }}>
              <div style={{ display: 'grid', gridTemplateColumns: '300px 1fr', gap: 20 }}>
                <ParetoChart candidates={report.candidates} />
                <div>
                  <div style={{ fontSize: 12, color: T.muted, lineHeight: 1.8 }}>
                    <strong style={{ color: T.cyan }}>How to read this chart:</strong><br />
                    Each dot = one gRNA candidate. X-axis = cutting efficiency. Y-axis = CAPS compatibility.<br /><br />
                    <strong style={{ color: T.cyan }}>Pareto-optimal</strong> (cyan dots) = no other gRNA is better on BOTH axes simultaneously.<br /><br />
                    The best choice depends on your experimental priority:
                    <br />• T1 screening is critical → pick the rightmost pareto point (high CAPS)
                    <br />• Editing efficiency is critical → pick the topmost pareto point (high efficiency)
                    <br />• Balanced → pick the pareto point closest to the top-right corner
                  </div>
                </div>
              </div>
            </div>
          )}

          {/* Table view */}
          {activeView === 'table' && (
            <div>
              {/* Column headers */}
              <div style={{
                display: 'grid',
                gridTemplateColumns: '36px 1fr 60px 60px 60px 100px 80px 40px',
                gap: 10, padding: '6px 14px',
                fontSize: 9, color: T.muted, letterSpacing: 1, textTransform: 'uppercase',
              }}>
                <div>#</div>
                <div>gRNA Sequence</div>
                <div style={{ textAlign: 'center' }}>Effic.</div>
                <div style={{ textAlign: 'center' }}>CAPS</div>
                <div style={{ textAlign: 'center' }}>Total</div>
                <div>Best CAPS</div>
                <div>Off-target</div>
                <div></div>
              </div>

              {report.candidates.map((cand, i) => (
                <GRNARow
                  key={i}
                  cand={cand}
                  rank={i + 1}
                  isExpanded={expandedIdx === i}
                  onToggle={() => setExpandedIdx(expandedIdx === i ? -1 : i)}
                />
              ))}
            </div>
          )}
        </div>
      )}

      {/* Error state */}
      {report && !report.success && (
        <div style={{ ...card, background: `${T.red}0C`, border: `1px solid ${T.red}44` }}>
          <div style={{ color: T.red, fontSize: 14 }}>⚠ {report.error}</div>
          <button style={{ ...btnSm, marginTop: 12, borderColor: T.red, color: T.red }}
            onClick={() => setReport(null)}>← Try Again</button>
        </div>
      )}
    </div>
  );
}
