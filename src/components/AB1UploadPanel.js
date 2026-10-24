// ============================================================
// FILE: src/components/AB1UploadPanel.js
// VERSION: 2.0 — User-Controlled Trimming (Phase 2)
//
// Four trimming modes:
//   • Standard  — Phred ≥ 20, window 10 bp (lab default)
//   • Quality   — user sets Phred threshold (Q10–Q40) + window
//   • Manual    — drag exact start/end sliders per sample
//   • No Trim   — raw full sequence
//
// Quality bar visualisation with trim-boundary markers and
// cut-site proximity warning built in.
// ============================================================

import React, { useState, useRef, useCallback, useMemo } from 'react';
import { readMultipleAB1Files, ab1ToSample, applyTrim } from '../engine/ab1Parser';

const T = {
  bg:'#060D1F', card:'#0D1B35', border:'#1a3060',
  cyan:'#00D4AA', blue:'#4D9FFF', purple:'#BF5AF2',
  gold:'#F0C040', red:'#FF4D6D',
  text:'#E8F0FE', muted:'#8899BB', dim:'#3A4A6A',
};

const TRIM_MODES = [
  { id:'standard', label:'Standard', icon:'⚙️', desc:'Phred ≥ 20, window 10 bp', color:T.cyan },
  { id:'quality',  label:'Quality',  icon:'🎚', desc:'Set your own Phred threshold', color:T.blue },
  { id:'manual',   label:'Manual',   icon:'✂️', desc:'Drag start/end sliders', color:T.purple },
  { id:'none',     label:'No Trim',  icon:'🔓', desc:'Use full raw sequence', color:T.gold },
];

export default function AB1UploadPanel({ onSamplesLoaded, generation, ploidy, cutSite = null }) {
  const [isDragging,   setDragging]   = useState(false);
  const [isProcessing, setProcessing] = useState(false);
  const [parsedFiles,  setParsed]     = useState([]);
  const [error,        setError]      = useState(null);
  const [globalMode,   setGlobalMode] = useState('standard');
  const [globalOpts,   setGlobalOpts] = useState({ phredThreshold:20, window:10 });
  const [sampleModes,  setSampleModes] = useState({});
  const fileInputRef = useRef(null);

  const emitSamples = useCallback((results, modeMap, defMode, defOpts) => {
    const samples = results
      .filter(r => r.success)
      .map(r => {
        const ov   = modeMap[r.filename] || {};
        const mode = ov.mode    || defMode;
        const opts = ov.options || defOpts;
        return ab1ToSample(r, generation || 'T0', ploidy || 'diploid', mode, opts);
      })
      .filter(Boolean);
    if (samples.length > 0 && onSamplesLoaded) onSamplesLoaded(samples, results);
  }, [generation, ploidy, onSamplesLoaded]); // eslint-disable-line react-hooks/exhaustive-deps

  const processFiles = useCallback(async (files) => {
    if (!files || files.length === 0) return;
    const ab1Files = Array.from(files).filter(f => f.name.toLowerCase().endsWith('.ab1'));
    if (ab1Files.length === 0) { setError('No .ab1 files found.'); return; }
    setError(null); setProcessing(true);
    try {
      const results = await readMultipleAB1Files(ab1Files);
      setParsed(results);
      emitSamples(results, sampleModes, globalMode, globalOpts);
    } catch (err) { setError(`Error: ${err.message}`); }
    finally { setProcessing(false); }
  }, [generation, ploidy, sampleModes, globalMode, globalOpts, emitSamples]); // eslint-disable-line react-hooks/exhaustive-deps

  const handleGlobalMode = (mode) => {
    setGlobalMode(mode);
    if (parsedFiles.length > 0) emitSamples(parsedFiles, sampleModes, mode, globalOpts);
  };
  const handleGlobalOpts = (opts) => {
    const next = { ...globalOpts, ...opts };
    setGlobalOpts(next);
    if (parsedFiles.length > 0) emitSamples(parsedFiles, sampleModes, globalMode, next);
  };
  const handleSampleMode = (filename, mode, opts) => {
    const next = { ...sampleModes, [filename]: { mode, options: opts || globalOpts } };
    setSampleModes(next);
    if (parsedFiles.length > 0) emitSamples(parsedFiles, next, globalMode, globalOpts);
  };

  const onDrop      = useCallback((e) => { e.preventDefault(); setDragging(false); processFiles(e.dataTransfer.files); }, [processFiles]);
  const onDragOver  = (e) => { e.preventDefault(); setDragging(true); };
  const onDragLeave = () => setDragging(false);
  const onFileInput = (e) => processFiles(e.target.files);

  const ok   = parsedFiles.filter(r => r.success).length;
  const fail = parsedFiles.filter(r => !r.success).length;

  return (
    <div style={{ fontFamily:'Calibri, sans-serif', display:'flex', flexDirection:'column', gap:12 }}>

      {/* Drop zone */}
      <div onDrop={onDrop} onDragOver={onDragOver} onDragLeave={onDragLeave}
        onClick={() => fileInputRef.current?.click()}
        style={{ border:`2px dashed ${isDragging?T.cyan:T.border}`, borderRadius:12, padding:'26px 20px',
          textAlign:'center', cursor:'pointer', background:isDragging?'rgba(0,212,170,0.05)':T.card, transition:'all 0.2s' }}>
        <input ref={fileInputRef} type="file" multiple accept=".ab1" onChange={onFileInput} style={{ display:'none' }} />
        {isProcessing ? (
          <div>
            <div style={{ fontSize:26, marginBottom:8 }}>⚙️</div>
            <div style={{ color:T.cyan, fontSize:13, fontWeight:'bold' }}>Parsing .ab1 files...</div>
          </div>
        ) : (
          <div>
            <div style={{ fontSize:30, marginBottom:8 }}>{isDragging?'📂':'🧬'}</div>
            <div style={{ color:T.text, fontSize:13, fontWeight:'bold', marginBottom:5 }}>
              {isDragging?'Drop .ab1 files here':'Upload Sanger Sequencing Files'}
            </div>
            <div style={{ color:T.muted, fontSize:11, marginBottom:10 }}>Drag & drop .ab1 files or click to browse</div>
            <div style={{ display:'flex', justifyContent:'center', gap:8, flexWrap:'wrap' }}>
              {['.ab1 binary','Batch upload','Real traces + quality'].map(p=>(
                <span key={p} style={{ background:'rgba(0,212,170,0.08)', border:'1px solid rgba(0,212,170,0.2)', borderRadius:20, padding:'2px 9px', fontSize:10, color:T.cyan }}>{p}</span>
              ))}
            </div>
          </div>
        )}
      </div>

      {error && (
        <div style={{ background:'rgba(255,77,109,0.1)', border:'1px solid rgba(255,77,109,0.4)', borderRadius:8, padding:'10px 14px', color:T.red, fontSize:12 }}>⚠ {error}</div>
      )}

      {/* Global trim controls */}
      {parsedFiles.length > 0 && (
        <TrimControls
          label="Default trimming — applies to all samples"
          mode={globalMode} options={globalOpts}
          onModeChange={handleGlobalMode} onOptsChange={handleGlobalOpts}
        />
      )}

      {/* Results */}
      {parsedFiles.length > 0 && (
        <div>
          <div style={{ display:'flex', gap:10, marginBottom:10 }}>
            {ok   > 0 && <div style={{ background:'rgba(0,212,170,0.1)',  border:'1px solid rgba(0,212,170,0.4)',  borderRadius:8, padding:'6px 14px', color:T.cyan, fontSize:12, flex:1, textAlign:'center' }}>✓ {ok} file{ok!==1?'s':''} parsed</div>}
            {fail > 0 && <div style={{ background:'rgba(255,77,109,0.1)', border:'1px solid rgba(255,77,109,0.4)', borderRadius:8, padding:'6px 14px', color:T.red,  fontSize:12, flex:1, textAlign:'center' }}>✗ {fail} failed</div>}
          </div>
          <div style={{ display:'flex', flexDirection:'column', gap:9 }}>
            {parsedFiles.map((f, i) => (
              <FileCard key={i} file={f}
                defaultMode={globalMode} defaultOpts={globalOpts}
                sampleOverride={sampleModes[f.filename]}
                onModeChange={(mode, opts) => handleSampleMode(f.filename, mode, opts)}
                cutSite={cutSite}
              />
            ))}
          </div>
        </div>
      )}

      {/* Help */}
      <div style={{ padding:'10px 14px', background:'rgba(77,159,255,0.06)', border:'1px solid rgba(77,159,255,0.2)', borderRadius:8 }}>
        <div style={{ color:T.blue, fontSize:11, fontWeight:'bold', marginBottom:4 }}>💡 Choosing a trim mode</div>
        <div style={{ color:T.muted, fontSize:11, lineHeight:1.7 }}>
          <b style={{ color:T.text }}>Standard Q20:</b> best for most routine Sanger reads. &nbsp;
          <b style={{ color:T.text }}>Quality Q25–30:</b> use for high-quality publication data. &nbsp;
          <b style={{ color:T.text }}>Manual:</b> use when the cut site is near read ends and auto-trim removes it. &nbsp;
          <b style={{ color:T.text }}>No Trim:</b> when your read is already clean or you want full control.
        </div>
      </div>
    </div>
  );
}

// ── Trim Controls widget ──────────────────────────────────
function TrimControls({ label, mode, options, onModeChange, onOptsChange }) {
  return (
    <div style={{ background:T.card, border:`1px solid ${T.border}`, borderRadius:10, padding:'12px 14px' }}>
      <div style={{ fontSize:10, color:T.muted, letterSpacing:1, textTransform:'uppercase', marginBottom:10 }}>✂️ {label}</div>
      <div style={{ display:'flex', gap:6, marginBottom:10, flexWrap:'wrap' }}>
        {TRIM_MODES.map(m => (
          <button key={m.id} onClick={() => onModeChange(m.id)} title={m.desc} style={{
            background: mode===m.id ? `${m.color}20` : 'transparent',
            border: `1px solid ${mode===m.id ? m.color : T.border}`,
            borderRadius:8, padding:'5px 12px', cursor:'pointer',
            color: mode===m.id ? m.color : T.muted,
            fontSize:12, fontWeight: mode===m.id ? 700 : 400, transition:'all 0.15s',
          }}>
            {m.icon} {m.label}
          </button>
        ))}
      </div>
      {mode === 'quality' && (
        <div style={{ display:'flex', gap:16, alignItems:'center', flexWrap:'wrap' }}>
          <div>
            <div style={{ fontSize:10, color:T.muted, marginBottom:3 }}>PHRED THRESHOLD</div>
            <div style={{ display:'flex', alignItems:'center', gap:8 }}>
              <input type="range" min="10" max="40" step="5" value={options.phredThreshold||20}
                onChange={e => onOptsChange({ phredThreshold:Number(e.target.value) })}
                style={{ width:120, accentColor:T.blue }} />
              <span style={{ fontFamily:'monospace', fontSize:13, color:T.blue, fontWeight:700 }}>Q{options.phredThreshold||20}</span>
            </div>
          </div>
          <div>
            <div style={{ fontSize:10, color:T.muted, marginBottom:3 }}>WINDOW</div>
            <div style={{ display:'flex', alignItems:'center', gap:8 }}>
              <input type="range" min="5" max="25" step="5" value={options.window||10}
                onChange={e => onOptsChange({ window:Number(e.target.value) })}
                style={{ width:100, accentColor:T.blue }} />
              <span style={{ fontFamily:'monospace', fontSize:13, color:T.blue, fontWeight:700 }}>{options.window||10} bp</span>
            </div>
          </div>
          <span style={{ fontSize:11, color:T.muted }}>Error rate ≤ {(Math.pow(10,-(options.phredThreshold||20)/10)*100).toFixed(2)}%</span>
        </div>
      )}
      {mode === 'manual' && (
        <div style={{ fontSize:11, color:T.purple, background:'rgba(191,90,242,0.08)', border:'1px solid rgba(191,90,242,0.25)', borderRadius:6, padding:'6px 10px' }}>
          ✂️ Expand each file card below to set exact trim positions per sample.
        </div>
      )}
      {mode === 'none' && (
        <div style={{ fontSize:11, color:T.gold, background:'rgba(240,192,64,0.08)', border:'1px solid rgba(240,192,64,0.25)', borderRadius:6, padding:'6px 10px' }}>
          🔓 Full raw sequence used — low-quality ends may reduce alignment accuracy.
        </div>
      )}
    </div>
  );
}

// ── Per-file card ─────────────────────────────────────────
function FileCard({ file, defaultMode, defaultOpts, sampleOverride, onModeChange, cutSite }) {
  const [expanded, setExpanded] = useState(false);
  const localMode = sampleOverride?.mode    || defaultMode;
  const localOpts = sampleOverride?.options || defaultOpts;

  const trimResult = useMemo(() => {
    if (!file.success) return null;
    return applyTrim(file.sequence, file.quality, localMode, localOpts, cutSite);
  }, [file, localMode, localOpts, cutSite]);

  const qg = file.qualityGrade;

  return (
    <div style={{ background:T.card, border:`1px solid ${file.success?T.border:'rgba(255,77,109,0.4)'}`, borderRadius:10, overflow:'hidden' }}>
      {/* Header */}
      <div onClick={() => file.success && setExpanded(e=>!e)}
        style={{ padding:'9px 14px', display:'flex', justifyContent:'space-between', alignItems:'center', cursor:file.success?'pointer':'default' }}>
        <div style={{ display:'flex', alignItems:'center', gap:10 }}>
          <span style={{ fontSize:15 }}>{file.success?'🧬':'❌'}</span>
          <div>
            <div style={{ fontSize:13, fontWeight:'bold', color:file.success?T.text:T.red, fontFamily:'monospace' }}>
              {file.filename||'Unknown'}
            </div>
            {file.success && trimResult && (
              <div style={{ fontSize:11, color:T.muted, marginTop:1 }}>
                {file.seqLength} bp raw → <span style={{ color:T.cyan, fontWeight:600 }}>{trimResult.trimmedSeq.length} bp trimmed</span>
                {' · '}{TRIM_MODES.find(m=>m.id===localMode)?.icon} {localMode}
                {file.runDate && ` · ${file.runDate}`}
              </div>
            )}
            {!file.success && <div style={{ fontSize:11, color:T.red, marginTop:1 }}>{file.error}</div>}
          </div>
        </div>
        <div style={{ display:'flex', alignItems:'center', gap:8 }}>
          {trimResult?.warning && <span title={trimResult.warning} style={{ fontSize:12, cursor:'help' }}>⚠️</span>}
          {file.success && qg && (
            <span style={{ background:`${qg.color}20`, border:`1px solid ${qg.color}60`, borderRadius:12, padding:'2px 9px', color:qg.color, fontSize:10, fontWeight:'bold' }}>{qg.grade}</span>
          )}
          {file.success && <span style={{ color:T.muted, fontSize:11 }}>{expanded?'▲':'▼'}</span>}
        </div>
      </div>

      {/* Expanded */}
      {expanded && file.success && trimResult && (
        <div style={{ borderTop:`1px solid ${T.border}`, padding:'12px 14px' }}>

          {/* Cut-site warning */}
          {trimResult.warning && (
            <div style={{ background:'rgba(240,192,64,0.1)', border:'1px solid rgba(240,192,64,0.35)', borderRadius:7, padding:'7px 12px', color:T.gold, fontSize:11, marginBottom:12, lineHeight:1.5 }}>
              {trimResult.warning}
            </div>
          )}

          {/* Per-sample mode */}
          <div style={{ marginBottom:12 }}>
            <div style={{ fontSize:10, color:T.muted, letterSpacing:1, textTransform:'uppercase', marginBottom:6 }}>Override trim for this sample</div>
            <div style={{ display:'flex', gap:6, flexWrap:'wrap' }}>
              {TRIM_MODES.map(m => (
                <button key={m.id} onClick={() => onModeChange(m.id, localOpts)} style={{
                  background: localMode===m.id?`${m.color}20`:'transparent',
                  border:`1px solid ${localMode===m.id?m.color:T.border}`,
                  borderRadius:7, padding:'4px 11px', cursor:'pointer',
                  color:localMode===m.id?m.color:T.muted, fontSize:11, fontWeight:localMode===m.id?700:400,
                }}>{m.icon} {m.label}</button>
              ))}
            </div>
          </div>

          {/* Quality threshold sliders */}
          {localMode === 'quality' && (
            <div style={{ display:'flex', gap:14, marginBottom:12, flexWrap:'wrap', alignItems:'center' }}>
              <div>
                <div style={{ fontSize:10, color:T.muted, marginBottom:3 }}>PHRED THRESHOLD</div>
                <div style={{ display:'flex', alignItems:'center', gap:7 }}>
                  <input type="range" min="10" max="40" step="5" value={localOpts.phredThreshold||20}
                    onChange={e => onModeChange('quality', {...localOpts, phredThreshold:Number(e.target.value)})}
                    style={{ width:110, accentColor:T.blue }} />
                  <span style={{ fontFamily:'monospace', fontSize:12, color:T.blue, fontWeight:700 }}>Q{localOpts.phredThreshold||20}</span>
                </div>
              </div>
              <div>
                <div style={{ fontSize:10, color:T.muted, marginBottom:3 }}>WINDOW</div>
                <div style={{ display:'flex', alignItems:'center', gap:7 }}>
                  <input type="range" min="5" max="25" step="5" value={localOpts.window||10}
                    onChange={e => onModeChange('quality', {...localOpts, window:Number(e.target.value)})}
                    style={{ width:90, accentColor:T.blue }} />
                  <span style={{ fontFamily:'monospace', fontSize:12, color:T.blue, fontWeight:700 }}>{localOpts.window||10} bp</span>
                </div>
              </div>
            </div>
          )}

          {/* Manual sliders */}
          {localMode === 'manual' && (
            <div style={{ marginBottom:12 }}>
              <div style={{ display:'flex', gap:14, flexWrap:'wrap', alignItems:'center' }}>
                <div>
                  <div style={{ fontSize:10, color:T.muted, marginBottom:3 }}>TRIM START (bp)</div>
                  <div style={{ display:'flex', alignItems:'center', gap:7 }}>
                    <input type="range" min="0" max={Math.max(0, file.seqLength-10)}
                      value={localOpts.manualStart ?? trimResult.trimStart}
                      onChange={e => onModeChange('manual', {...localOpts, manualStart:Number(e.target.value)})}
                      style={{ width:130, accentColor:T.purple }} />
                    <span style={{ fontFamily:'monospace', fontSize:12, color:T.purple, fontWeight:700, minWidth:32 }}>
                      {localOpts.manualStart ?? trimResult.trimStart}
                    </span>
                  </div>
                </div>
                <div>
                  <div style={{ fontSize:10, color:T.muted, marginBottom:3 }}>TRIM END (bp)</div>
                  <div style={{ display:'flex', alignItems:'center', gap:7 }}>
                    <input type="range"
                      min={Math.min(file.seqLength, (localOpts.manualStart ?? trimResult.trimStart) + 10)}
                      max={file.seqLength}
                      value={localOpts.manualEnd ?? trimResult.trimEnd}
                      onChange={e => onModeChange('manual', {...localOpts, manualEnd:Number(e.target.value)})}
                      style={{ width:130, accentColor:T.purple }} />
                    <span style={{ fontFamily:'monospace', fontSize:12, color:T.purple, fontWeight:700, minWidth:32 }}>
                      {localOpts.manualEnd ?? trimResult.trimEnd}
                    </span>
                  </div>
                </div>
              </div>
              <div style={{ fontSize:11, color:T.muted, marginTop:5 }}>
                Region: <span style={{ color:T.purple, fontFamily:'monospace' }}>{trimResult.trimStart}–{trimResult.trimEnd}</span>
                {' '}· Length: <span style={{ color:T.purple, fontFamily:'monospace' }}>{trimResult.trimmedSeq.length} bp</span>
              </div>
            </div>
          )}

          {/* Quality bar */}
          {file.quality?.length > 0 && (
            <QualityBar quality={file.quality} trimStart={trimResult.trimStart} trimEnd={trimResult.trimEnd} cutSite={cutSite} seqLength={file.seqLength} />
          )}

          {/* Stats */}
          <div style={{ display:'grid', gridTemplateColumns:'repeat(3,1fr)', gap:8, marginTop:10 }}>
            {[
              { l:'Raw Length',   v:`${file.seqLength} bp` },
              { l:'Trimmed',      v:`${trimResult.trimmedSeq.length} bp`, c:T.cyan },
              { l:'Avg Quality',  v:`Q${file.avgQuality}`, c:file.avgQuality>=30?T.cyan:file.avgQuality>=20?T.gold:T.red },
              { l:'Trim Start',   v:`pos ${trimResult.trimStart}` },
              { l:'Trim End',     v:`pos ${trimResult.trimEnd}` },
              { l:'Traces',       v:file.syntheticTraces?'Synthetic':'Real .ab1', c:file.syntheticTraces?T.gold:T.cyan },
            ].map(({ l, v, c }) => (
              <div key={l}>
                <div style={{ fontSize:9, color:T.muted, textTransform:'uppercase', letterSpacing:'0.5px', marginBottom:2 }}>{l}</div>
                <div style={{ fontSize:12, color:c||T.text, fontFamily:'monospace', fontWeight:'bold' }}>{v}</div>
              </div>
            ))}
          </div>

          {/* Sequence preview */}
          <div style={{ marginTop:10 }}>
            <div style={{ fontSize:9, color:T.muted, textTransform:'uppercase', letterSpacing:'0.5px', marginBottom:4 }}>Trimmed Sequence Preview</div>
            <div style={{ fontFamily:'monospace', fontSize:11, color:T.cyan, letterSpacing:'1px', background:'#060D1F', padding:'6px 10px', borderRadius:6, overflowX:'auto', whiteSpace:'nowrap' }}>
              {trimResult.trimmedSeq.slice(0, 80)}
              {trimResult.trimmedSeq.length > 80 && <span style={{ color:T.muted }}>...+{trimResult.trimmedSeq.length - 80} more</span>}
            </div>
          </div>
        </div>
      )}
    </div>
  );
}

// ── Quality bar ───────────────────────────────────────────
function QualityBar({ quality, trimStart, trimEnd, cutSite, seqLength }) {
  const n      = seqLength || quality.length;
  const pts    = Math.min(n, 200);
  const step   = n / pts;

  const bars = Array.from({ length:pts }, (_, i) => {
    const idx = Math.round(i * step);
    const q   = quality[idx] || 0;
    return { q, pct:(idx/n)*100, inTrim: idx>=trimStart && idx<trimEnd };
  });

  const qcol  = q => q>=30?T.cyan : q>=20?T.blue : q>=10?T.gold : T.red;
  const sP    = (trimStart/n)*100;
  const eP    = (trimEnd/n)*100;
  const cutP  = cutSite!=null ? (cutSite/n)*100 : null;

  return (
    <div style={{ marginBottom:8 }}>
      <div style={{ fontSize:9, color:T.muted, textTransform:'uppercase', letterSpacing:'0.5px', marginBottom:4 }}>
        Quality trace — <span style={{ color:T.cyan }}>■</span> kept &nbsp;<span style={{ color:T.dim }}>■</span> trimmed{cutP!=null && <> &nbsp;<span style={{ color:T.red }}>│</span> cut site</>}
      </div>
      <div style={{ position:'relative', height:44, background:'rgba(0,0,0,0.3)', borderRadius:6, overflow:'hidden' }}>
        <div style={{ position:'absolute', left:0,     top:0, width:`${sP}%`,      height:'100%', background:'rgba(0,0,0,0.45)' }} />
        <div style={{ position:'absolute', left:`${eP}%`, top:0, width:`${100-eP}%`, height:'100%', background:'rgba(0,0,0,0.45)' }} />
        {bars.map((b, i) => (
          <div key={i} style={{
            position:'absolute', left:`${b.pct}%`, bottom:0,
            width:`${100/pts}%`, height:`${Math.max(4,(b.q/60)*100)}%`,
            background: b.inTrim ? qcol(b.q) : 'rgba(100,120,150,0.3)',
            opacity: b.inTrim ? 0.85 : 0.4,
          }} />
        ))}
        <div style={{ position:'absolute', left:`${sP}%`,  top:0, width:2, height:'100%', background:'rgba(191,90,242,0.85)' }} />
        <div style={{ position:'absolute', left:`${eP}%`,  top:0, width:2, height:'100%', background:'rgba(191,90,242,0.85)' }} />
        {cutP!=null && <div style={{ position:'absolute', left:`${cutP}%`, top:0, width:2, height:'100%', background:'rgba(255,77,109,0.9)', zIndex:2 }} />}
        <div style={{ position:'absolute', left:`${sP+0.5}%`, top:2, fontSize:8, color:'rgba(191,90,242,0.9)' }}>{trimStart}</div>
        <div style={{ position:'absolute', left:`${Math.max(0,eP-4)}%`, top:2, fontSize:8, color:'rgba(191,90,242,0.9)' }}>{trimEnd}</div>
        {cutP!=null && <div style={{ position:'absolute', left:`${cutP+0.5}%`, bottom:2, fontSize:8, color:'rgba(255,77,109,0.9)' }}>cut</div>}
      </div>
      <div style={{ display:'flex', justifyContent:'space-between', marginTop:2, fontSize:9, color:T.dim }}>
        <span>0</span><span style={{ color:T.red }}>Q10</span><span style={{ color:T.gold }}>Q20</span><span style={{ color:T.blue }}>Q30</span><span style={{ color:T.cyan }}>Q40+</span><span>{n} bp</span>
      </div>
    </div>
  );
}
