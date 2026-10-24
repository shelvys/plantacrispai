// ============================================================
// FILE: src/components/MLQualityPanel.js
//
// REPLACES the old ML Model tab.
// Shows: per-sample confidence, quality indicators, R² summary,
// ensemble agreement — WITHOUT exposing training accuracy internals.
//
// Props:
//   results     — array of analysed samples from analysisEngine
//   onSelect    — callback(sample)
//   selected    — currently selected sample
// ============================================================

import React, { useState, useMemo } from 'react';

const DARK_BG  = '#0a1628';
const PANEL_BG = '#0d1b35';
const BORDER   = '#1a3060';
const TEXT      = '#c8d6f0';
const ACCENT    = '#00d4aa';

const CLASS_CONFIG = {
  0: { label: 'Wild-type',          color: '#4d9fff', icon: '○' },
  1: { label: 'Low / Mosaic Edit',  color: '#f59e0b', icon: '◑' },
  2: { label: 'Heterozygous',       color: '#a855f7', icon: '◐' },
  3: { label: 'Homozygous Edited',  color: '#00d4aa', icon: '●' }
};

// ─── Confidence bar ────────────────────────────────────────────────────────────
function ConfBar({ value, color = ACCENT, height = 6 }) {
  return (
    <div style={{
      background: '#111c33', borderRadius: height / 2,
      height, width: '100%', overflow: 'hidden'
    }}>
      <div style={{
        width: `${value}%`, height: '100%',
        background: color, borderRadius: height / 2,
        transition: 'width 0.4s ease'
      }} />
    </div>
  );
}

// ─── Quality indicator dot ─────────────────────────────────────────────────────
function QualityDot({ r2 }) {
  const color = r2 >= 0.90 ? ACCENT : r2 >= 0.75 ? '#f59e0b' : '#ff6b6b';
  const label = r2 >= 0.90 ? 'High' : r2 >= 0.75 ? 'Med' : 'Low';
  return (
    <span style={{
      display: 'inline-flex', alignItems: 'center', gap: 4,
      color, fontSize: 10, fontWeight: 600
    }}>
      <span style={{
        width: 7, height: 7, borderRadius: '50%', background: color, display: 'inline-block'
      }} />
      {label}
    </span>
  );
}

// ─── Sample card ───────────────────────────────────────────────────────────────
function SampleCard({ sample, isSelected, onClick }) {
  const cls   = CLASS_CONFIG[sample.zygosityClass] || CLASS_CONFIG[2];
  const r2    = sample.r2 ?? 0;
  const score = sample.mlScore ?? 0;

  return (
    <div
      onClick={onClick}
      style={{
        background: isSelected ? `${cls.color}14` : DARK_BG,
        border: `1px solid ${isSelected ? cls.color + '55' : BORDER}`,
        borderRadius: 8, padding: '10px 14px', cursor: 'pointer',
        transition: 'all 0.15s', marginBottom: 6
      }}
    >
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 6 }}>
        <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
          <span style={{ color: cls.color, fontSize: 14 }}>{cls.icon}</span>
          <span style={{ color: TEXT, fontSize: 13, fontWeight: 600 }}>
            {sample.name || sample.id}
          </span>
          <span style={{
            background: '#111c33', color: '#4d6a9a',
            fontSize: 10, padding: '1px 6px', borderRadius: 10
          }}>
            {sample.generation}
          </span>
        </div>
        <div style={{ display: 'flex', alignItems: 'center', gap: 10 }}>
          <QualityDot r2={r2} />
          <span style={{
            color: cls.color, fontSize: 13, fontWeight: 700, minWidth: 36, textAlign: 'right'
          }}>
            {score}%
          </span>
        </div>
      </div>
      <ConfBar value={score} color={cls.color} />
      <div style={{ display: 'flex', justifyContent: 'space-between', marginTop: 5 }}>
        <span style={{ color: '#4d6a9a', fontSize: 10 }}>{cls.label}</span>
        <span style={{ color: '#4d6a9a', fontSize: 10 }}>
          R² {r2.toFixed(2)} · ICE {sample.iceScore ?? '--'}%
          {sample.flag ? ' · ⚠' : ''}
        </span>
      </div>
    </div>
  );
}

// ─── Detail panel ─────────────────────────────────────────────────────────────
function SampleDetail({ sample }) {
  if (!sample) return (
    <div style={{
      display: 'flex', alignItems: 'center', justifyContent: 'center',
      height: 300, color: '#4d6a9a', fontSize: 13, flexDirection: 'column', gap: 8
    }}>
      <span style={{ fontSize: 28 }}>🔬</span>
      <span>Select a sample to view analysis details</span>
    </div>
  );

  const cls   = CLASS_CONFIG[sample.zygosityClass] || CLASS_CONFIG[2];
  const r2    = sample.r2 ?? 0;
  const score = sample.mlScore ?? 0;

  // Confidence breakdown bars
  const features = [
    { label: 'Sequence Quality (R²)', value: Math.round(r2 * 100), color: r2 > 0.8 ? ACCENT : '#f59e0b' },
    { label: 'Edit Signal (ICE)',     value: sample.iceScore ?? 0,  color: '#4d9fff' },
    { label: 'Indel Clarity',         value: sample.indelSize > 0 ? Math.min(100, sample.indelSize * 8 + 40) : 0, color: '#a855f7' },
    { label: 'Alignment Score',       value: Math.round((sample.sequenceIdentity ?? 0.5) * 100), color: '#f59e0b' },
  ];

  return (
    <div>
      {/* Sample header */}
      <div style={{
        background: `${cls.color}14`, border: `1px solid ${cls.color}44`,
        borderRadius: 10, padding: '14px 18px', marginBottom: 14
      }}>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
          <div>
            <div style={{ color: cls.color, fontSize: 16, fontWeight: 700 }}>
              {cls.icon} {sample.name || sample.id}
            </div>
            <div style={{ color: cls.color + 'aa', fontSize: 12, marginTop: 3 }}>
              {cls.label}
            </div>
          </div>
          <div style={{ textAlign: 'right' }}>
            <div style={{ color: cls.color, fontSize: 28, fontWeight: 800 }}>{score}%</div>
            <div style={{ color: '#4d6a9a', fontSize: 10 }}>ML confidence</div>
          </div>
        </div>

        <div style={{ marginTop: 12 }}>
          <ConfBar value={score} color={cls.color} height={8} />
        </div>
      </div>

      {/* Feature contributions */}
      <div style={{
        background: DARK_BG, borderRadius: 10, padding: '14px 18px',
        border: `1px solid ${BORDER}`, marginBottom: 14
      }}>
        <div style={{ color: '#4d6a9a', fontSize: 11, fontWeight: 700, marginBottom: 12, letterSpacing: 1 }}>
          SIGNAL CONTRIBUTIONS
        </div>
        {features.map(f => (
          <div key={f.label} style={{ marginBottom: 10 }}>
            <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: 4 }}>
              <span style={{ color: TEXT, fontSize: 12 }}>{f.label}</span>
              <span style={{ color: f.color, fontSize: 12, fontWeight: 600 }}>{f.value}%</span>
            </div>
            <ConfBar value={f.value} color={f.color} />
          </div>
        ))}
      </div>

      {/* Edit summary */}
      <div style={{
        background: DARK_BG, borderRadius: 10, padding: '14px 18px',
        border: `1px solid ${BORDER}`
      }}>
        <div style={{ color: '#4d6a9a', fontSize: 11, fontWeight: 700, marginBottom: 12, letterSpacing: 1 }}>
          EDIT SUMMARY
        </div>
        <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: 10 }}>
          {[
            { k: 'Edit Type',      v: sample.editType || 'None' },
            { k: 'Indel Size',     v: sample.indelSize ? `${Math.abs(sample.indelSize)} bp` : '—' },
            { k: 'R² Quality',     v: r2.toFixed(3) },
            { k: 'ICE Score',      v: `${sample.iceScore ?? '--'}%` },
            { k: 'Frameshift',     v: sample.isFrameshift ? 'Yes' : 'No' },
            { k: 'Sequence ID',    v: sample.sequenceIdentity ? `${(sample.sequenceIdentity*100).toFixed(1)}%` : '--' }
          ].map(({ k, v }) => (
            <div key={k}>
              <div style={{ color: '#4d6a9a', fontSize: 10 }}>{k}</div>
              <div style={{ color: TEXT, fontSize: 13, fontWeight: 600 }}>{v}</div>
            </div>
          ))}
        </div>

        <div style={{ marginTop: 12, paddingTop: 12, borderTop: `1px solid ${BORDER}` }}>
          <div style={{ color: '#4d6a9a', fontSize: 10, marginBottom: 4 }}>PROTEIN EFFECT</div>
          <div style={{
            color: sample.isFrameshift ? '#ff6b6b' : sample.indelSize > 0 ? '#f59e0b' : ACCENT,
            fontSize: 12, fontWeight: 600
          }}>
            {sample.proteinEffect || 'No prediction available'}
          </div>
        </div>

        {sample.flag && (
          <div style={{
            marginTop: 12, background: '#ff6b6b12', borderRadius: 8, padding: '10px 12px',
            border: '1px solid #ff6b6b44'
          }}>
            <div style={{ color: '#ff6b6b', fontSize: 11, fontWeight: 600 }}>⚠ Re-sequencing recommended</div>
            <div style={{ color: '#ff6b6b88', fontSize: 10, marginTop: 2 }}>{sample.flag}</div>
          </div>
        )}
      </div>
    </div>
  );
}

// ─── Summary stats row ────────────────────────────────────────────────────────
function SummaryRow({ results }) {
  const total     = results.length;
  const highConf  = results.filter(r => (r.mlScore ?? 0) >= 80).length;
  const needsVerif= results.filter(r => r.flag).length;
  const avgR2     = total > 0 ? results.reduce((a, r) => a + (r.r2 ?? 0), 0) / total : 0;

  return (
    <div style={{
      display: 'grid', gridTemplateColumns: 'repeat(4, 1fr)', gap: 10, marginBottom: 18
    }}>
      {[
        { label: 'Total Samples',    value: total,              color: TEXT,    sub: 'analysed' },
        { label: 'High Confidence',  value: `${highConf}/${total}`, color: ACCENT, sub: '≥ 80% ML score' },
        { label: 'Verify First',     value: needsVerif,         color: needsVerif > 0 ? '#f59e0b' : ACCENT, sub: 'needs resequencing' },
        { label: 'Avg R² Quality',   value: avgR2.toFixed(3),   color: avgR2 > 0.8 ? ACCENT : '#f59e0b', sub: 'sequence quality' }
      ].map(stat => (
        <div key={stat.label} style={{
          background: DARK_BG, border: `1px solid ${BORDER}`,
          borderRadius: 10, padding: '12px 16px'
        }}>
          <div style={{ color: '#4d6a9a', fontSize: 10, marginBottom: 4 }}>{stat.label}</div>
          <div style={{ color: stat.color, fontSize: 22, fontWeight: 800 }}>{stat.value}</div>
          <div style={{ color: '#4d6a9a', fontSize: 10 }}>{stat.sub}</div>
        </div>
      ))}
    </div>
  );
}

// ─── Main export ──────────────────────────────────────────────────────────────
export default function MLQualityPanel({ results = [], onSelect, selected }) {
  const [sortBy, setSortBy] = useState('mlScore');

  const sorted = useMemo(() => {
    const arr = [...results];
    if (sortBy === 'mlScore')   arr.sort((a, b) => (b.mlScore ?? 0) - (a.mlScore ?? 0));
    if (sortBy === 'r2')        arr.sort((a, b) => (b.r2 ?? 0) - (a.r2 ?? 0));
    if (sortBy === 'ice')       arr.sort((a, b) => (b.iceScore ?? 0) - (a.iceScore ?? 0));
    if (sortBy === 'name')      arr.sort((a, b) => (a.name || '').localeCompare(b.name || ''));
    if (sortBy === 'zygosity')  arr.sort((a, b) => (b.zygosityClass ?? 0) - (a.zygosityClass ?? 0));
    return arr;
  }, [results, sortBy]);

  if (!results || results.length === 0) {
    return (
      <div style={{
        background: PANEL_BG, border: `1px solid ${BORDER}`, borderRadius: 12,
        padding: 40, textAlign: 'center', color: '#4d6a9a'
      }}>
        <div style={{ fontSize: 36, marginBottom: 12 }}>🧬</div>
        <div style={{ fontSize: 16, color: TEXT, marginBottom: 8 }}>No analysis results yet</div>
        <div style={{ fontSize: 13 }}>Run analysis in the Setup tab to see ML quality scores here</div>
      </div>
    );
  }

  return (
    <div style={{ display: 'flex', gap: 16 }}>
      {/* Left — sample list */}
      <div style={{ flex: 1, minWidth: 0 }}>
        <div style={{
          background: PANEL_BG, border: `1px solid ${BORDER}`, borderRadius: 12,
          padding: 18, height: '100%'
        }}>
          {/* Header */}
          <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', marginBottom: 16 }}>
            <div>
              <div style={{ color: ACCENT, fontWeight: 700, fontSize: 15 }}>
                🧠 Analysis Quality
              </div>
              <div style={{ color: '#4d6a9a', fontSize: 11, marginTop: 2 }}>
                Ensemble classifier · ICE + R² + alignment signals
              </div>
            </div>
            <select
              value={sortBy}
              onChange={e => setSortBy(e.target.value)}
              style={{
                background: DARK_BG, border: `1px solid ${BORDER}`,
                color: TEXT, borderRadius: 6, padding: '4px 8px', fontSize: 11
              }}
            >
              <option value="mlScore">Sort: ML Score</option>
              <option value="r2">Sort: R² Quality</option>
              <option value="ice">Sort: ICE Score</option>
              <option value="zygosity">Sort: Zygosity</option>
              <option value="name">Sort: Name</option>
            </select>
          </div>

          <SummaryRow results={results} />

          {/* Sample cards */}
          <div style={{ maxHeight: 440, overflowY: 'auto', paddingRight: 4 }}>
            {sorted.map(s => (
              <SampleCard
                key={s.id || s.name}
                sample={s}
                isSelected={selected?.id === s.id || selected?.name === s.name}
                onClick={() => onSelect && onSelect(s)}
              />
            ))}
          </div>
        </div>
      </div>

      {/* Right — detail panel */}
      <div style={{ width: 320, flexShrink: 0 }}>
        <div style={{
          background: PANEL_BG, border: `1px solid ${BORDER}`, borderRadius: 12, padding: 18
        }}>
          <SampleDetail sample={selected} />
        </div>
      </div>
    </div>
  );
}
