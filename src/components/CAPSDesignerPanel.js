// ============================================================
// FILE: src/components/CAPSDesignerPanel.js
// PURPOSE: Complete CAPS Designer UI — all three modes
//
// Mode A: Post-editing, no CAPS planned → discover/design
// Mode B: Post-editing, CAPS planned    → validate
// Mode C: Pre-editing                   → gRNA + CAPS design
// ============================================================

import React, { useState, useMemo } from 'react';
import { runCAPSAnalysis, generateGelData } from '../engine/capsEngine';
import GelSimulator from './GelSimulator';
import GRNACAPSDesignerPanel from './GRNACAPSDesignerPanel';

// ── Theme ──────────────────────────────────────────────────
const T = {
  bg:         '#060D1F',
  card:       '#0D1B35',
  cardBorder: '#1a3060',
  cyan:       '#00D4AA',
  blue:       '#4D9FFF',
  purple:     '#B060FF',
  gold:       '#F0C040',
  red:        '#FF4D6D',
  green:      '#00D4AA',
  text:       '#E8F0FE',
  muted:      '#8899BB',
  mutedDark:  '#3A4A6A',
};

const card   = { background: T.card, border: `1px solid ${T.cardBorder}`, borderRadius: 14, padding: 20 };
const input  = { background: '#060D1F', border: `1px solid ${T.cardBorder}`, borderRadius: 8, color: T.text, padding: '9px 13px', fontSize: 13, width: '100%', fontFamily: 'Calibri, sans-serif', outline: 'none', boxSizing: 'border-box' };
const btn    = { background: `linear-gradient(135deg, ${T.cyan}, #00a882)`, color: '#000', border: 'none', borderRadius: 10, padding: '10px 24px', fontSize: 13, fontWeight: 700, fontFamily: 'Calibri, sans-serif', cursor: 'pointer' };
const btnSec = { background: 'transparent', color: T.cyan, border: `1px solid ${T.cyan}`, borderRadius: 10, padding: '8px 18px', fontSize: 12, fontWeight: 600, fontFamily: 'Calibri, sans-serif', cursor: 'pointer' };
const label  = { fontSize: 10, color: T.muted, letterSpacing: 1, textTransform: 'uppercase', marginBottom: 5, display: 'block' };

// ── Score Badge ────────────────────────────────────────────
function ScoreBadge({ score, size = 'md' }) {
  const color = score >= 75 ? T.cyan : score >= 55 ? T.blue : score >= 35 ? T.gold : T.red;
  const fs    = size === 'lg' ? 22 : 14;
  return (
    <div style={{ display: 'inline-flex', flexDirection: 'column', alignItems: 'center',
                  background: `${color}18`, border: `1px solid ${color}55`,
                  borderRadius: 10, padding: size === 'lg' ? '10px 16px' : '4px 10px' }}>
      <span style={{ fontSize: fs, fontWeight: 800, color, fontFamily: 'monospace' }}>{score}%</span>
      {size === 'lg' && (
        <span style={{ fontSize: 9, color: T.muted, letterSpacing: 1 }}>
          {score >= 75 ? 'EXCELLENT' : score >= 55 ? 'GOOD' : score >= 35 ? 'POSSIBLE' : 'POOR'}
        </span>
      )}
    </div>
  );
}

// ── Enzyme Chip ────────────────────────────────────────────
function EnzymeChip({ name, score, size = 'sm' }) {
  const color = score >= 75 ? T.cyan : score >= 55 ? T.blue : T.gold;
  return (
    <span style={{ background: `${color}18`, border: `1px solid ${color}55`,
                   borderRadius: 20, padding: size === 'lg' ? '5px 14px' : '2px 10px',
                   color, fontSize: size === 'lg' ? 14 : 11, fontWeight: 700 }}>
      {name}
    </span>
  );
}

// ── Main Component ─────────────────────────────────────────
export default function CAPSDesignerPanel({ results, cfg, refSequence }) {
  const [mode, setMode]                     = useState(null);
  const [capsReport, setCapsReport]         = useState(null);
  const [selectedSample, setSelectedSample] = useState(null);
  const [selectedCAPS, setSelectedCAPS]     = useState(null);
  const [running, setRunning]               = useState(false);
  const [ampliconLen, setAmpliconLen]       = useState(500);

  const hasResults = results && results.length > 0;
  const hasRef     = refSequence && refSequence.length > 20;

  // ── Run CAPS Analysis (Modes A & B) ──────────────────────
  const runCAPS = async () => {
    if (!hasRef || !hasResults) return;
    setRunning(true);
    await new Promise(r => setTimeout(r, 300));

    const report = runCAPSAnalysis(results, refSequence, {
      ...cfg,
      ampliconLength: ampliconLen,
    });
    setCapsReport(report);
    setSelectedSample(report.samples[0] || null);
    setSelectedCAPS(report.samples[0]?.bestOption || null);
    setRunning(false);
  };

  // ── Gel data ──────────────────────────────────────────────
  const gelData = useMemo(() => {
    if (!selectedCAPS || !capsReport) return [];
    return generateGelData(selectedCAPS, capsReport.samples, ampliconLen);
  }, [selectedCAPS, capsReport, ampliconLen]);

  // ── Mode Selector ─────────────────────────────────────────
  if (!mode) {
    return (
      <div style={{ fontFamily: 'Calibri, sans-serif' }}>
        <div style={{ ...card, marginBottom: 20, background:
          'linear-gradient(135deg, rgba(0,212,170,0.06), rgba(176,96,255,0.06))' }}>
          <div style={{ fontSize: 22, fontWeight: 800, color: T.text, marginBottom: 6 }}>
            🔍 CAPS Marker Designer
          </div>
          <div style={{ color: T.muted, fontSize: 13, lineHeight: 1.6 }}>
            Cleaved Amplified Polymorphic Sequence (CAPS) markers allow rapid,
            low-cost screening of hundreds of T1 plants without sequencing each one.
            Select your situation below to get started.
          </div>
        </div>

        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(3, 1fr)', gap: 16 }}>
          {[
            {
              id: 'A',
              icon: '🔬',
              title: 'I have sequencing results',
              subtitle: 'Find or design CAPS from my indels',
              desc: 'Already edited and sequenced T0 plants. No CAPS was planned in advance. PlantaCRISPAI will scan all RE sites and design dCAPS primers if needed.',
              color: T.cyan,
              available: hasResults && hasRef,
              unavailMsg: 'Requires: run analysis first (Setup tab)',
            },
            {
              id: 'B',
              icon: '✅',
              title: 'I planned a CAPS site',
              subtitle: 'Validate it worked after sequencing',
              desc: 'Deliberately designed gRNA near an RE site. Now verify the editing disrupted the CAPS marker as expected, with gel band predictions per sample.',
              color: T.blue,
              available: hasResults && hasRef,
              unavailMsg: 'Requires: run analysis first (Setup tab)',
            },
            {
              id: 'C',
              icon: '🌱',
              title: "I haven't edited yet",
              subtitle: 'Design gRNA + CAPS together',
              desc: 'Pre-editing planning. Input your target gene sequence and get ranked gRNA candidates, each scored for both cutting efficiency AND built-in CAPS compatibility.',
              color: T.purple,
              available: true,
              unavailMsg: '',
            },
          ].map(m => (
            <div
              key={m.id}
              onClick={() => m.available && setMode(m.id)}
              style={{
                ...card,
                cursor: m.available ? 'pointer' : 'not-allowed',
                opacity: m.available ? 1 : 0.5,
                border: `1px solid ${m.color}44`,
                transition: 'all 0.2s',
                position: 'relative',
                overflow: 'hidden',
              }}
            >
              {/* Color accent top bar */}
              <div style={{ position: 'absolute', top: 0, left: 0, right: 0,
                            height: 3, background: m.color }} />

              <div style={{ fontSize: 32, marginBottom: 10 }}>{m.icon}</div>
              <div style={{ fontSize: 14, fontWeight: 800, color: m.color, marginBottom: 4 }}>
                Mode {m.id} — {m.title}
              </div>
              <div style={{ fontSize: 12, color: T.muted, marginBottom: 10, fontStyle: 'italic' }}>
                {m.subtitle}
              </div>
              <div style={{ fontSize: 12, color: T.text, lineHeight: 1.6, marginBottom: 14 }}>
                {m.desc}
              </div>
              {!m.available && (
                <div style={{ fontSize: 11, color: T.gold, background: `${T.gold}18`,
                              border: `1px solid ${T.gold}44`, borderRadius: 6,
                              padding: '4px 10px' }}>
                  ⚠ {m.unavailMsg}
                </div>
              )}
              {m.available && (
                <div style={{ fontSize: 12, color: m.color, fontWeight: 700 }}>
                  Select → 
                </div>
              )}
            </div>
          ))}
        </div>
      </div>
    );
  }

  // ══════════════════════════════════════════════════════════
  // MODE A & B — Post-editing CAPS analysis
  // ══════════════════════════════════════════════════════════
  if (mode === 'A' || mode === 'B') {
    return (
      <div style={{ fontFamily: 'Calibri, sans-serif' }}>

        {/* Header */}
        <div style={{ display:'flex', justifyContent:'space-between', alignItems:'center',
                      marginBottom: 20, flexWrap:'wrap', gap: 10 }}>
          <div>
            <div style={{ fontSize: 18, fontWeight: 800, color: T.text }}>
              {mode === 'A' ? '🔬 CAPS Discovery' : '✅ CAPS Validation'}
            </div>
            <div style={{ fontSize: 12, color: T.muted, marginTop: 3 }}>
              {mode === 'A'
                ? 'Finding restriction sites disrupted by your edits'
                : 'Validating your planned CAPS marker worked correctly'}
            </div>
          </div>
          <button onClick={() => { setMode(null); setCapsReport(null); }} style={btnSec}>
            ← Change Mode
          </button>
        </div>

        {/* Config + Run */}
        {!capsReport && (
          <div style={{ ...card, marginBottom: 20 }}>
            <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: 16, marginBottom: 16 }}>
              <div>
                <label style={label}>Expected PCR Amplicon Size (bp)</label>
                <input style={input} type="number" value={ampliconLen}
                  onChange={e => setAmpliconLen(Number(e.target.value))}
                  placeholder="500" />
                <div style={{ fontSize: 11, color: T.muted, marginTop: 4 }}>
                  Size of your diagnostic PCR product (not the sequencing amplicon)
                </div>
              </div>
              <div>
                <label style={label}>Reference Sequence</label>
                <div style={{ fontSize: 12, color: T.cyan, padding: '9px 0' }}>
                  ✓ {refSequence?.length || 0} bp loaded from Setup tab
                </div>
              </div>
            </div>

            <button style={{ ...btn, opacity: running ? 0.6 : 1 }}
              onClick={runCAPS} disabled={running}>
              {running ? '⚙ Scanning RE sites...' : `🔍 Run CAPS Analysis (${results?.length || 0} samples)`}
            </button>
          </div>
        )}

        {/* Results */}
        {capsReport && (
          <div>
            {/* Summary cards */}
            <div style={{ display:'grid', gridTemplateColumns:'repeat(4,1fr)', gap:14, marginBottom:20 }}>
              {[
                { label:'Best Enzyme',      val: capsReport.summary.bestEnzyme,   color: T.cyan },
                { label:'CAPS Score',       val: `${capsReport.summary.bestScore}%`, color: T.cyan },
                { label:'Natural CAPS',     val: capsReport.summary.withNaturalCAPS, color: T.blue },
                { label:'dCAPS Designed',   val: capsReport.summary.withDCAPS,    color: T.purple },
              ].map((s,i) => (
                <div key={i} style={{ ...card, textAlign:'center' }}>
                  <div style={{ fontSize: 28, fontWeight:800, color: s.color, lineHeight:1.1 }}>{s.val}</div>
                  <div style={{ fontSize:10, color: T.muted, letterSpacing:1, textTransform:'uppercase', marginTop:4 }}>{s.label}</div>
                </div>
              ))}
            </div>

            <div style={{ display:'grid', gridTemplateColumns:'280px 1fr', gap:20, marginBottom:20 }}>

              {/* Sample list */}
              <div style={{ ...card }}>
                <div style={{ fontSize:13, fontWeight:700, color:T.cyan, marginBottom:12 }}>
                  🔬 Samples
                </div>
                {capsReport.samples.map((s,i) => (
                  <div key={i}
                    onClick={() => { setSelectedSample(s); setSelectedCAPS(s.bestOption); }}
                    style={{
                      padding:'10px 12px', borderRadius:10, cursor:'pointer', marginBottom:6,
                      background: selectedSample?.sampleId === s.sampleId
                        ? `${T.cyan}18` : '#060D1F',
                      border: `1px solid ${selectedSample?.sampleId === s.sampleId ? T.cyan : T.cardBorder}`,
                    }}>
                    <div style={{ fontWeight:700, fontSize:13, color:T.text }}>{s.sampleId}</div>
                    <div style={{ fontSize:11, color:T.muted, marginTop:2 }}>
                      {s.editType} · {s.zygosity}
                    </div>
                    <div style={{ marginTop:6, display:'flex', gap:6, flexWrap:'wrap' }}>
                      {s.hasCaps && <EnzymeChip name={s.caps[0]?.enzyme} score={s.caps[0]?.capsScore || 0} />}
                      {!s.hasCaps && s.hasDcaps && (
                        <span style={{ fontSize:10, color:T.purple }}>dCAPS: {s.dcaps[0]?.enzyme}</span>
                      )}
                      {!s.hasCaps && !s.hasDcaps && (
                        <span style={{ fontSize:10, color:T.muted }}>No CAPS found</span>
                      )}
                    </div>
                  </div>
                ))}
              </div>

              {/* Sample detail */}
              {selectedSample && (
                <div>
                  <SampleCAPSDetail
                    sample={selectedSample}
                    onSelectCAPS={setSelectedCAPS}
                    selectedCAPS={selectedCAPS}
                  />
                </div>
              )}
            </div>

            {/* Gel Simulator */}
            {selectedCAPS && (
              <div style={{ ...card }}>
                <div style={{ fontSize:14, fontWeight:700, color:T.cyan, marginBottom:14 }}>
                  🧫 Gel Band Simulator — {selectedCAPS.enzyme}
                </div>
                <GelSimulator
                  lanes={gelData}
                  enzyme={selectedCAPS.enzyme}
                  ampliconLength={ampliconLen}
                />
              </div>
            )}

            <div style={{ marginTop:14, textAlign:'center' }}>
              <button style={btnSec} onClick={() => setCapsReport(null)}>
                ← Reconfigure Analysis
              </button>
            </div>
          </div>
        )}
      </div>
    );
  }

  // ══════════════════════════════════════════════════════════
  // MODE C — Pre-editing gRNA + CAPS design (ML-powered)
  // Delegated to GRNACAPSDesignerPanel
  // ══════════════════════════════════════════════════════════
  if (mode === 'C') {
    return (
      <div style={{ fontFamily: 'Calibri, sans-serif' }}>
        <div style={{ display:'flex', justifyContent:'space-between',
                      alignItems:'center', marginBottom: 20, flexWrap:'wrap', gap:10 }}>
          <div>
            <div style={{ fontSize: 18, fontWeight: 800, color: T.text }}>
              🌱 ML-Powered gRNA + CAPS Designer
            </div>
            <div style={{ fontSize: 12, color: T.muted, marginTop: 3 }}>
              Feature-engineered ensemble · Nearest-neighbor thermodynamics ·
              Doench 2016 position weights · Plant-specific corrections
            </div>
          </div>
          <button onClick={() => setMode(null)} style={btnSec}>
            ← Change Mode
          </button>
        </div>
        <GRNACAPSDesignerPanel refSequence={refSequence} />
      </div>
    );
  }

  return null;
}

// ── Sample CAPS Detail Sub-component ──────────────────────

function SampleCAPSDetail({ sample, onSelectCAPS, selectedCAPS }) {
  return (
    <div style={{ display:'flex', flexDirection:'column', gap:14 }}>

      {/* Natural CAPS */}
      {sample.hasCaps && (
        <div style={{ background: T.card, border:`1px solid ${T.cardBorder}`,
                      borderRadius:12, padding:16 }}>
          <div style={{ fontSize:13, fontWeight:700, color:T.cyan, marginBottom:12 }}>
            ✂ Natural CAPS Markers Found
          </div>
          {sample.caps.map((c,i) => (
            <CAPSCard key={i} caps={c} type="natural"
              isSelected={selectedCAPS?.enzyme === c.enzyme && selectedCAPS?.sitePosition === c.sitePosition}
              onClick={() => onSelectCAPS(c)} />
          ))}
        </div>
      )}

      {/* dCAPS */}
      {!sample.hasCaps && sample.hasDcaps && (
        <div style={{ background: T.card, border:`1px solid ${T.cardBorder}`,
                      borderRadius:12, padding:16 }}>
          <div style={{ fontSize:13, fontWeight:700, color:T.purple, marginBottom:12 }}>
            🔧 dCAPS Primers Designed
            <span style={{ fontSize:11, color:T.muted, fontWeight:400, marginLeft:8 }}>
              No natural CAPS found — primers designed with deliberate mismatch
            </span>
          </div>
          {sample.dcaps.map((d,i) => (
            <DcapsCard key={i} dcaps={d}
              isSelected={selectedCAPS?.enzyme === d.enzyme}
              onClick={() => onSelectCAPS(d)} />
          ))}
        </div>
      )}

      {/* No CAPS */}
      {!sample.hasCaps && !sample.hasDcaps && (
        <div style={{ background:`${T.gold}0C`, border:`1px solid ${T.gold}44`,
                      borderRadius:12, padding:16 }}>
          <div style={{ color:T.gold, fontSize:13 }}>
            ⚠ No CAPS or dCAPS found near the cut site for this sample.
          </div>
          <div style={{ color:T.muted, fontSize:12, marginTop:6 }}>
            Consider: longer amplicon to include more RE sites,
            or use the dCAPS method manually (dCAPS Finder 2.0 tool).
          </div>
        </div>
      )}
    </div>
  );
}

// ── CAPS Card ──────────────────────────────────────────────
function CAPSCard({ caps, type, isSelected, onClick }) {
  const color = caps.capsScore >= 75 ? T.cyan : caps.capsScore >= 55 ? T.blue : T.gold;
  return (
    <div onClick={onClick} style={{
      background: isSelected ? `${color}12` : '#060D1F',
      border: `1px solid ${isSelected ? color : T.mutedDark}`,
      borderRadius:10, padding:'12px 14px', marginBottom:10,
      cursor:'pointer', transition:'all 0.15s',
    }}>
      <div style={{ display:'flex', justifyContent:'space-between', alignItems:'center', marginBottom:8 }}>
        <div style={{ display:'flex', alignItems:'center', gap:10 }}>
          <EnzymeChip name={caps.enzyme} score={caps.capsScore} size="lg" />
          <div>
            <div style={{ fontSize:11, color:T.muted }}>
              Site at position {caps.sitePosition} · {caps.mechanism?.replace(/_/g,' ')}
            </div>
            <div style={{ fontSize:11, color:T.muted, fontFamily:'monospace', marginTop:2 }}>
              {caps.matchedSeq}
            </div>
          </div>
        </div>
        <ScoreBadge score={caps.capsScore} size="lg" />
      </div>
      <div style={{ display:'grid', gridTemplateColumns:'1fr 1fr', gap:10 }}>
        <div style={{ background:T.bg, borderRadius:8, padding:'8px 12px' }}>
          <div style={{ fontSize:9, color:T.muted, letterSpacing:1, marginBottom:4 }}>WILD-TYPE BANDS</div>
          <div style={{ fontSize:13, color:T.blue, fontWeight:700, fontFamily:'monospace' }}>
            {caps.bands?.wt?.band1} + {caps.bands?.wt?.band2} bp
          </div>
          <div style={{ fontSize:10, color:T.muted }}>2 bands (cut)</div>
        </div>
        <div style={{ background:T.bg, borderRadius:8, padding:'8px 12px' }}>
          <div style={{ fontSize:9, color:T.muted, letterSpacing:1, marginBottom:4 }}>EDITED BANDS</div>
          <div style={{ fontSize:13, color:T.red, fontWeight:700, fontFamily:'monospace' }}>
            Uncut ({caps.bands?.wt ? caps.bands.wt.band1 + caps.bands.wt.band2 : '—'} bp)
          </div>
          <div style={{ fontSize:10, color:T.muted }}>1 band (uncut)</div>
        </div>
      </div>
      <div style={{ marginTop:8, fontSize:11, color:color, fontWeight:600 }}>
        → {caps.recommendation}
      </div>
    </div>
  );
}

// ── dCAPS Card ─────────────────────────────────────────────
function DcapsCard({ dcaps, isSelected, onClick }) {
  return (
    <div onClick={onClick} style={{
      background: isSelected ? `${T.purple}12` : '#060D1F',
      border: `1px solid ${isSelected ? T.purple : T.mutedDark}`,
      borderRadius:10, padding:'12px 14px', marginBottom:10,
      cursor:'pointer', transition:'all 0.15s',
    }}>
      <div style={{ display:'flex', justifyContent:'space-between', marginBottom:8 }}>
        <div>
          <EnzymeChip name={dcaps.enzyme} score={dcaps.dcapsScore} size="lg" />
          <span style={{ marginLeft:8, fontSize:11, color:T.purple }}>dCAPS</span>
        </div>
        <ScoreBadge score={dcaps.dcapsScore} size="lg" />
      </div>

      <div style={{ background:T.bg, borderRadius:8, padding:'10px 12px', marginBottom:10 }}>
        <div style={{ fontSize:9, color:T.muted, letterSpacing:1, marginBottom:6 }}>PRIMER SEQUENCE (order this)</div>
        <code style={{ fontSize:13, color:T.cyan, letterSpacing:'1.5px', fontFamily:'monospace',
          wordBreak:'break-all', lineHeight:1.8 }}>
          {dcaps.primerSeq}
        </code>
        <div style={{ display:'flex', gap:16, marginTop:8 }}>
          {[
            { label:'Length', val:`${dcaps.primerLen} bp` },
            { label:'Tm',     val:`${dcaps.tm}°C` },
            { label:'GC%',    val:`${dcaps.gc}%` },
            { label:'Mismatches', val:`${dcaps.mismatches}` },
          ].map(({label, val}) => (
            <div key={label}>
              <div style={{ fontSize:9, color:T.muted, letterSpacing:1 }}>{label}</div>
              <div style={{ fontSize:12, color:T.text, fontWeight:700 }}>{val}</div>
            </div>
          ))}
        </div>
      </div>

      <div style={{ fontSize:11, color:T.muted, marginBottom:6 }}>
        {dcaps.detail}
      </div>
      <div style={{ fontSize:11, color:T.gold }}>
        ⚠ {dcaps.orderingNote}
      </div>
      <div style={{ marginTop:6, fontSize:11, color:T.purple, fontWeight:600 }}>
        → {dcaps.recommendation}
      </div>
    </div>
  );
}
