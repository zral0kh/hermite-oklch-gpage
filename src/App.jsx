import React, { useState, useMemo, useEffect, useRef } from 'react';
import { X, Plus, RotateCcw, MousePointer2, GitCommit, Copy, Check, Lock, Unlock, Monitor, FileCode, Braces, Layers, Repeat, Globe, Box, GripVertical, Settings2, AlignCenter, AlignJustify, Ruler, Compass, ArrowLeftRight } from 'lucide-react';

// ==========================================
// 1. Math & Logic Library
// ==========================================

const PI = Math.PI;
const EPS = 1e-12;

const lin = (v) => (v <= 0.04045 ? v / 12.92 : Math.pow((v + 0.055) / 1.055, 2.4));
const gam = (v) => (v <= 0.0031308 ? 12.92 * v : 1.055 * Math.pow(v, 1 / 2.4) - 0.055);

const vScale = (v, s) => v.map((x) => x * s);
const vAdd = (a, b) => a.map((x, i) => x + (b?.[i] ?? 0));
const vSub = (a, b) => a.map((x, i) => x - (b?.[i] ?? 0));
const vDot = (a, b) => a.reduce((sum, x, i) => sum + x * (b?.[i] ?? 0), 0);
const vLen = (v) => Math.sqrt(v.reduce((sum, x) => sum + x*x, 0));

// Distance excluding alpha (only first 3 components used)
function vDist3(a, b) {
  let sum = 0;
  for (let i = 0; i < 3; i++) {
    const d = (a[i] || 0) - (b[i] || 0);
    sum += d * d;
  }
  return Math.sqrt(sum);
}

/* ----------------- Color Conversions ----------------- */

function rgbToOklab(vec) {
  const [r0, g0, b0, alpha] = vec;
  const r = lin(r0), g = lin(g0), b = lin(b0);
  const l = 0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b;
  const m = 0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b;
  const s = 0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b;
  const l_ = Math.cbrt(l), m_ = Math.cbrt(m), s_ = Math.cbrt(s);
  const out = [
    0.2104542553 * l_ + 0.793617785 * m_ - 0.0040720468 * s_,
    1.9779984951 * l_ - 2.428592205 * m_ + 0.4505937099 * s_,
    0.0259040371 * l_ + 0.7827717662 * m_ - 0.808675766 * s_,
  ];
  if (alpha !== undefined) out.push(alpha);
  return out;
}

function oklabToRgb(vec) {
  const [L, a, b, alpha] = vec;
  const l_ = L + 0.3963377774 * a + 0.2158037573 * b;
  const m_ = L - 0.1055613458 * a - 0.0638541728 * b;
  const s_ = L - 0.0894841775 * a - 1.291485548 * b;
  const l = l_ ** 3, m = m_ ** 3, s = s_ ** 3;
  const r = gam(4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s);
  const g = gam(-1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s);
  const b_ = gam(-0.0041960863 * l - 0.7034186147 * m + 1.707614701 * s);
  const out = [r, g, b_];
  if (alpha !== undefined) out.push(alpha);
  return out;
}

function oklabToOklch(vec) {
  const [l, a, b, alpha] = vec;
  const c = Math.hypot(a, b);
  const h = (Math.atan2(b, a) * 180) / PI;
  const res = [l, c, (h + 360) % 360];
  if (alpha !== undefined) res.push(alpha);
  return res;
}

function oklchToOklab(vec) {
  const [l, c, h, alpha] = vec;
  const rad = (h * PI) / 180;
  const out = [l, c * Math.cos(rad), c * Math.sin(rad)];
  if (alpha !== undefined) out.push(alpha);
  return out;
}

function oklchToRgb(vec) {
    return oklabToRgb(oklchToOklab(vec));
}

// HSB Helpers
function rgbToHsb(r, g, b) {
  const max = Math.max(r, g, b), min = Math.min(r, g, b);
  let h, s, v = max;
  const d = max - min;
  s = max === 0 ? 0 : d / max;
  if (max === min) h = 0; 
  else {
    switch (max) {
      case r: h = (g - b) / d + (g < b ? 6 : 0); break;
      case g: h = (b - r) / d + 2; break;
      case b: h = (r - g) / d + 4; break;
    }
    h /= 6;
  }
  return [h * 360, s, v];
}

function hsbToRgb(h, s, b) {
  const k = (n) => (n + h / 60) % 6;
  const f = (n) => b * (1 - s * Math.max(0, Math.min(k(n), 4 - k(n), 1)));
  return [f(5), f(3), f(1)];
}
function hslToHsb(h, s, l) {
  s /= 100; l /= 100;
  const b = l + s * Math.min(l, 1 - l);
  const sb = b === 0 ? 0 : 2 * (1 - l / b);
  return [h, sb, b];
}
function hsbToHsl(h, s, b) {
  const l = b * (1 - s / 2);
  let sl = 0;
  if (l > 0 && l < 1) sl = (b - l) / Math.min(l, 1 - l);
  return { h, s: sl * 100, l: l * 100 };
}

/* ===================== Logic Utils ===================== */

function parseColor(input, targetSpace = 'oklab') {
  if (Array.isArray(input)) return input.slice();
  return input; 
}

function unwrapAngles(arr) {
  if (!arr || arr.length === 0) return [];
  const out = arr.slice();
  for (let i = 1; i < out.length; i++) {
    let delta = out[i] - out[i - 1];
    while (delta <= -180) delta += 360;
    while (delta > 180) delta -= 360;
    out[i] = out[i - 1] + delta;
  }
  return out;
}
function wrapAngle(a) { return ((a % 360) + 360) % 360; }

function getDelta(pB, pA, space) {
  const d = pB.map((v, i) => v - (pA[i] || 0));
  if (space === 'oklch') {
    let dh = d[2];
    while (dh > 180) dh -= 360;
    while (dh < -180) dh += 360;
    d[2] = dh;
  }
  return d;
}

/* ===================== Solvers ===================== */

function solveTridiagonal(n, a, b, c, d) {
  const cp = new Array(n);
  const dp = new Array(n);
  const x = new Array(n);
  cp[0] = c[0] / b[0];
  dp[0] = d[0] / b[0];
  for (let i = 1; i < n; i++) {
    const denom = b[i] - a[i] * cp[i - 1];
    if (Math.abs(denom) < EPS) return d.slice();
    cp[i] = c[i] / denom;
    dp[i] = (d[i] - a[i] * dp[i - 1]) / denom;
  }
  x[n - 1] = dp[n - 1];
  for (let i = n - 2; i >= 0; i--) x[i] = dp[i] - cp[i] * x[i + 1];
  return x;
}

function solveCyclicTridiagonal(n, a, b, c, rhs) {
  if (n <= 2) return solveTridiagonal(n, a, b, c, rhs);
  const gamma = -b[0];
  const bb = b.slice();
  bb[0] -= gamma;
  bb[n - 1] -= (c[n - 1] * a[0]) / gamma;
  const y = solveTridiagonal(n, a, bb, c, rhs);
  const u = new Array(n).fill(0);
  u[0] = gamma;
  u[n - 1] = c[n - 1];
  const z = solveTridiagonal(n, a, bb, c, u);
  const factNum = y[0] + (a[0] / gamma) * y[n - 1];
  const factDen = 1 + z[0] + (a[0] / gamma) * z[n - 1];
  const fact = factNum / factDen;
  return y.map((yi, i) => yi - fact * z[i]);
}

function fitNaturalCubic(P, chords, loop, space) {
  const n = P.length;
  const dims = P[0].length;
  const M = new Array(n).fill(0).map(() => new Array(dims).fill(0));

  for (let k = 0; k < dims; k++) {
    let vals = P.map(p => p[k] ?? 0);
    if (space === 'oklch' && k === 2) vals = unwrapAngles(vals);

    const a = new Array(n).fill(0);
    const b = new Array(n).fill(0);
    const c = new Array(n).fill(0);
    const r = new Array(n).fill(0);

    if (loop) {
      for (let i = 0; i < n; i++) {
        const prev = (i - 1 + n) % n;
        const next = (i + 1) % n;
        const hPrev = chords[prev];
        const hNext = chords[i];
        a[i] = hNext;
        b[i] = 2 * (hNext + hPrev);
        c[i] = hPrev;
        let valPrev = vals[prev], valCurr = vals[i], valNext = vals[next];
        if (space === 'oklch' && k === 2) {
          if (i === 0) {
            let d = valPrev - valCurr;
            while (d > 180) d -= 360; while (d < -180) d += 360;
            valPrev = valCurr + d;
          }
          if (i === n - 1) {
            let d = valNext - valCurr;
            while (d > 180) d -= 360; while (d < -180) d += 360;
            valNext = valCurr + d;
          }
        }
        r[i] = 3 * (hNext * (valCurr - valPrev) / hPrev + hPrev * (valNext - valCurr) / hNext);
      }
      const res = solveCyclicTridiagonal(n, a, b, c, r);
      for (let i = 0; i < n; i++) M[i][k] = res[i];
    } else {
      for (let i = 1; i < n - 1; i++) {
        const hPrev = chords[i - 1];
        const hNext = chords[i];
        a[i] = hNext; b[i] = 2 * (hNext + hPrev); c[i] = hPrev;
        r[i] = 3 * (hNext * (vals[i] - vals[i - 1]) / hPrev + hPrev * (vals[i + 1] - vals[i]) / hNext);
      }
      const h0 = chords[0];
      b[0] = 2 * h0; c[0] = h0; r[0] = 3 * (vals[1] - vals[0]);
      const hLast = chords[n - 2];
      a[n - 1] = hLast; b[n - 1] = 2 * hLast; r[n - 1] = 3 * (vals[n - 1] - vals[n - 2]);
      const res = solveTridiagonal(n, a, b, c, r);
      for (let i = 0; i < n; i++) M[i][k] = res[i];
    }
  }
  return M;
}

function fitCatmullRom(P, crChords, loop, space) {
  const n = P.length;
  const M = new Array(n);
  const segments = loop ? n : n - 1;

  for (let i = 0; i < n; i++) {
    if (!loop && (i === 0 || i === n - 1)) {
      if (i === 0) M[i] = getDelta(P[1], P[0], space);
      else M[i] = getDelta(P[n - 1], P[n - 2], space);
      continue;
    }
    const prevIdx = (i - 1 + n) % n;
    const nextIdx = (i + 1) % n;
    const dt0 = crChords[prevIdx % segments];
    const dt1 = crChords[i % segments];
    const v1 = getDelta(P[i], P[prevIdx], space);
    const v2 = getDelta(P[nextIdx], P[i], space);
    const w1 = dt1 * dt1;
    const w2 = dt0 * dt0;
    const denom = dt0 * dt1 * (dt0 + dt1) || EPS;
    M[i] = v1.map((val, k) => ((w1 * val + w2 * v2[k]) / denom) * dt1);
  }
  return M;
}

function formatColor(vec, format = 'oklab', inputSpace = 'oklab') {
  let lab = inputSpace === 'oklch' ? oklchToOklab(vec) : vec;
  if (format === 'oklab') return lab;
  if (format === 'oklch') {
    if (inputSpace === 'oklch') return { l: vec[0], c: vec[1], h: vec[2], a: vec[3] };
    const lch = oklabToOklch(lab);
    return { l: lch[0], c: lch[1], h: lch[2], a: lch[3] };
  }
  if (format === 'rgb') {
    const rgb = oklabToRgb(lab);
    return { r: Math.round(rgb[0]*255), g: Math.round(rgb[1]*255), b: Math.round(rgb[2]*255), a: rgb[3] };
  }
  return lab;
}

function findSegmentBinary(u, t) {
  const n = u.length;
  if (t <= u[0]) return 0;
  if (t >= u[n - 1]) return n - 2;
  let lo = 0, hi = n - 1;
  while (lo + 1 < hi) {
    const mid = (lo + hi) >> 1;
    if (u[mid] <= t) lo = mid;
    else hi = mid;
  }
  return Math.max(0, Math.min(n - 2, lo));
}

// Main Fit Function
function fitSpline(fixpoints, { space = 'oklab', method = 'natural-cubic', loop = false, strengths = 1, distribution = 'uniform', overrideTangents = [] } = {}) {
  const P = fixpoints.map((c) => parseColor(c, space));
  const n = P.length;
  if (n < 2) throw new Error('Need at least 2 points');

  const segments = loop ? n : n - 1;
  
  // Build chords based on distribution option
  let chords = [];
  
  if (distribution === 'uniform') {
    // Uniform spacing: all segments equal
    for (let i = 0; i < segments; i++) chords.push(1);
  } else {
    // Geometric spacing: based on spatial distance
    for (let i = 0; i < segments; i++) {
      const a = P[i];
      const b = P[(i + 1) % n];
      const d = getDelta(b, a, space);
      let dist = Math.hypot(d[0] ?? 0, d[1] ?? 0, d[2] ?? 0);
      if (dist < 1e-8) dist = 1e-8;
      
      // Apply method-specific scaling
      if (method === 'centripetal-CR') {
        chords.push(Math.sqrt(dist));
      } else if (method === 'chordal-CR') {
        chords.push(dist);
      } else {
        // natural-cubic uses geometric distance as-is
        chords.push(dist);
      }
    }
  }

  // Compute Tangents
  let M;
  if (method === 'natural-cubic') {
    M = fitNaturalCubic(P, chords, loop, space);
  } else {
    M = fitCatmullRom(P, chords, loop, space);
  }

  // Apply strengths & Manual Overrides
  const getS = (i) => (Array.isArray(strengths) ? (strengths[i] ?? 1) : strengths);
  M = M.map((t, i) => vScale(t, getS(i)));

  if (overrideTangents && overrideTangents.length) {
      for(let i=0; i<n; i++) {
          if (overrideTangents[i]) M[i] = overrideTangents[i];
      }
  }

  // Build Hermite Segments
  const segs = [];
  for (let i = 0; i < segments; i++) {
    const p0 = P[i];
    const p1 = P[(i + 1) % n];
    let m0 = M[i];
    let m1 = M[(i + 1) % n];

    // Convert tangents to Hermite form
    if (method === 'natural-cubic') {
      // Natural cubic returns derivatives; scale by chord length
      const h = chords[i];
      m0 = vScale(m0, h);
      m1 = vScale(m1, h);
    }
    // CR tangents are already in the correct units (no rescaling needed)

    let targetP1 = p1;
    if (space === 'oklch') {
      let dh = p1[2] - p0[2];
      while (dh > 180) dh -= 360;
      while (dh < -180) dh += 360;
      targetP1 = [...p1];
      targetP1[2] = p0[2] + dh;
    }

    const dims = Math.max(p0.length, targetP1.length);
    const a = new Array(dims), b = new Array(dims), c = new Array(dims), d = new Array(dims);
    for (let k = 0; k < dims; k++) {
      const P0k = p0[k] ?? 0;
      const P1k = targetP1[k] ?? 0;
      const m0k = m0[k] ?? 0;
      const m1k = m1[k] ?? 0;
      a[k] = P0k;
      b[k] = m0k;
      const delta = P1k - P0k;
      c[k] = 3 * delta - 2 * m0k - m1k;
      d[k] = -2 * delta + m0k + m1k;
    }
    segs.push({ a, b, c, d });
  }

  // Build cumulative U
  const u = [0];
  for (let i = 0; i < segments; i++) u.push(u[i] + chords[i]);

  return { segs, u, totalLen: u[u.length - 1], n, loop, space, method, distribution, tangents: M };
}

// Sampling Function
function splineColors(fixpoints, weights, { 
    format = 'oklab', 
    space = 'oklab', 
    method = 'natural-cubic', 
    loop = false, 
    strengths = 1, 
    distribution = 'uniform', 
    overrideTangents = [] 
} = {}) {
  const spline = fitSpline(fixpoints, { space, method, loop, strengths, distribution, overrideTangents });
  const segments = loop ? spline.n : spline.n - 1;

  return weights.map(w => {
    let ww = Math.max(0, Math.min(1, w));
    let segIndex = 0;
    let tLocal = 0;

    if (distribution === 'uniform') { 
      const scaled = ww * segments;
      segIndex = Math.floor(scaled);
      if (segIndex >= segments) segIndex = segments - 1;
      tLocal = scaled - segIndex;
    } else {
      const absT = ww * spline.totalLen;
      const i = findSegmentBinary(spline.u, absT);
      segIndex = i;
      const len = spline.u[i + 1] - spline.u[i];
      tLocal = len > 1e-12 ? (absT - spline.u[i]) / len : 0;
    }

    const seg = spline.segs[segIndex];
    const t2 = tLocal * tLocal, t3 = t2 * tLocal;
    const dims = seg.a.length;
    const out = new Array(dims);
    for (let k = 0; k < dims; k++) {
      out[k] = seg.a[k] + seg.b[k] * tLocal + seg.c[k] * t2 + seg.d[k] * t3;
    }
    if (space === 'oklch' && out[2] !== undefined) out[2] = wrapAngle(out[2]);
    return formatColor(out, format, space);
  });
}

// ==========================================
// 2. UI Components
// ==========================================

function Tooltip({ children, text, side = 'top' }) {
    const positionClasses = {
        top: 'bottom-full left-1/2 -translate-x-1/2 mb-2',
        bottom: 'top-full left-1/2 -translate-x-1/2 mt-2',
        left: 'right-full top-1/2 -translate-y-1/2 mr-2',
        right: 'left-full top-1/2 -translate-y-1/2 ml-2'
    };
    return (
        <div className="group relative flex items-center justify-center w-full h-full">
            {children}
            <div className={`
                absolute ${positionClasses[side]} 
                hidden group-hover:block z-[100]
                bg-gray-800 text-gray-200 text-[10px] font-medium
                px-2 py-1 rounded shadow-xl border border-gray-700
                whitespace-nowrap pointer-events-none
            `}>
                {text}
            </div>
        </div>
    );
}

// ==========================================
// Fixed Input Components
// ==========================================

function NumberInput({ value, onChange, min, max, label, step = 1 }) {
    // Local state to allow typing intermediate values (e.g. "-", "0.", empty string)
    const [strVal, setStrVal] = useState(String(value));
    const inputRef = useRef(null);

    // Sync with parent value ONLY if we are not currently editing (focused)
    useEffect(() => {
        if (document.activeElement !== inputRef.current) {
            // formatting: if step is decimal, show decimals, else integer
            const isInt = step >= 1;
            // Prevent showing excessive decimals if not needed, but keep precision if it exists
            const formatted = isInt ? Math.round(value) : Number(value).toFixed(3).replace(/\.?0+$/, '');
            setStrVal(String(formatted));
        }
    }, [value, step]);

    const commitChange = () => {
        let num = parseFloat(strVal);
        if (isNaN(num)) {
            // Revert to last valid prop value if input is garbage
            setStrVal(String(value));
            return;
        }
        
        // Clamp if min/max provided
        if (min !== undefined) num = Math.max(min, num);
        if (max !== undefined) num = Math.min(max, num);
        
        onChange(num);
    };

    const handleKeyDown = (e) => {
        if (e.key === 'Enter') {
            inputRef.current.blur(); // Triggers onBlur -> commitChange
        }
    };

    return (
        <div className="flex-1 min-w-0 group relative">
            <input 
                ref={inputRef}
                type="text" 
                value={strVal}
                onChange={(e) => setStrVal(e.target.value)}
                onBlur={commitChange}
                onKeyDown={handleKeyDown}
                className="w-full bg-gray-950 border border-gray-800 rounded px-1.5 py-0.5 text-[10px] font-mono text-gray-300 focus:border-blue-500 focus:outline-none text-right"
            />
            {label && <span className="absolute left-1 top-1/2 -translate-y-1/2 text-[8px] text-gray-600 pointer-events-none font-bold select-none">{label}</span>}
        </div>
    );
}

function HexInput({ value, onChange, label }) {
    const [strVal, setStrVal] = useState(value);
    const inputRef = useRef(null);

    useEffect(() => {
        if (document.activeElement !== inputRef.current) {
            setStrVal(value);
        }
    }, [value]);

    const commitChange = () => {
        // Validate Hex: #RRGGBB
        const match = strVal.match(/^#?([0-9A-Fa-f]{6})$/);
        if (match) {
            const formatted = `#${match[1].toUpperCase()}`;
            setStrVal(formatted);
            onChange(formatted);
        } else {
            // Revert if invalid
            setStrVal(value);
        }
    };

    const handleKeyDown = (e) => {
        if (e.key === 'Enter') inputRef.current.blur();
    };

    return (
      <div className="flex-1 min-w-0 group relative">
          <input
              ref={inputRef}
              type="text"
              value={strVal}
              onChange={(e) => setStrVal(e.target.value)}
              onBlur={commitChange}
              onKeyDown={handleKeyDown}
              className="w-full bg-gray-950 border border-gray-800 rounded px-1.5 py-0.5 text-[10px] font-mono text-gray-300 focus:border-blue-500 focus:outline-none text-right uppercase"
          />
        {label && <span className="absolute left-1 top-1/2 -translate-y-1/2 text-[8px] text-gray-600 pointer-events-none font-bold select-none">{label}</span>}
      </div>
    );
}

function ColorInputs({ color, onChange }) {
    const [h, s, b, a = 1] = color;
    const [r, g, b_rgb] = hsbToRgb(h, s, b);
    
    // Clamp RGB to 0-255 for display
    const rgb = { 
        r: Math.round(r * 255), 
        g: Math.round(g * 255), 
        b: Math.round(b_rgb * 255) 
    };
    
    // Hex calculation
    const hex = `#${((1 << 24) + (rgb.r << 16) + (rgb.g << 8) + rgb.b).toString(16).slice(1).toUpperCase()}`;
    
    const lab = rgbToOklab([r, g, b_rgb, a]);
    const oklchArr = oklabToOklch(lab);
    const oklch = { l: oklchArr[0], c: oklchArr[1], h: oklchArr[2] };
    const [copied, setCopied] = useState(null);

    const updateFromRgb = (k, v) => {
        const nextRgb = { ...rgb, [k]: v };
        onChange([...rgbToHsb(nextRgb.r/255, nextRgb.g/255, nextRgb.b/255), a]);
    };
    
    const updateFromHex = (v) => {
        const clean = v.replace('#', '');
        if(clean.length === 6) {
            const r = parseInt(clean.slice(0, 2), 16) / 255;
            const g = parseInt(clean.slice(2, 4), 16) / 255;
            const b = parseInt(clean.slice(4, 6), 16) / 255;
            onChange([...rgbToHsb(r, g, b), a]);
        }
    };

    const updateFromHsl = (k, v) => onChange([...hslToHsb({ ...hsl, [k]: v }.h, { ...hsl, [k]: v }.s, { ...hsl, [k]: v }.l), a]);
    const updateFromOklch = (k, v) => {
        const n = { ...oklch, [k]: v };
        const labNext = oklchToOklab([n.l, n.c, n.h]);
        const rgbNext = oklabToRgb(labNext);
        onChange([...rgbToHsb(rgbNext[0], rgbNext[1], rgbNext[2]), a]);
    };
    const updateAlpha = (v) => onChange([h, s, b, v]);

    const copy = (txt, key) => {
        navigator.clipboard.writeText(txt);
        setCopied(key);
        setTimeout(() => setCopied(null), 1500);
    };

    const hsl = hsbToHsl(h, s, b);

    return (
        <div className="flex flex-col gap-0.5 mt-2 p-1.5 bg-gray-950/50 rounded border border-gray-800" onDragStart={(e) => { e.preventDefault(); e.stopPropagation(); }} draggable>
            <div className="grid grid-cols-[24px_1fr_1fr_1fr_14px] gap-1 items-center">
                <span className="text-[9px] font-bold text-gray-600">RGB</span>
                <NumberInput value={rgb.r} min={0} max={255} onChange={(v) => updateFromRgb('r', v)} />
                <NumberInput value={rgb.g} min={0} max={255} onChange={(v) => updateFromRgb('g', v)} />
                <NumberInput value={rgb.b} min={0} max={255} onChange={(v) => updateFromRgb('b', v)} />
                <button onClick={() => copy(`rgba(${rgb.r}, ${rgb.g}, ${rgb.b}, ${a.toFixed(2)})`, 'rgb')} className="text-gray-500 hover:text-white">{copied === 'rgb' ? <Check size={10} className="text-green-500"/> : <Copy size={10}/>}</button>
            </div>
             <div className="grid grid-cols-[24px_1fr_1fr_1fr_14px] gap-1 items-center">
                <span className="text-[9px] font-bold text-gray-600">LCH</span>
                <NumberInput value={oklch.l} step={0.01} min={0} max={1} onChange={(v) => updateFromOklch('l', v)} />
                <NumberInput value={oklch.c} step={0.01} min={0} max={0.4} onChange={(v) => updateFromOklch('c', v)} />
                <NumberInput value={oklch.h} min={0} max={360} onChange={(v) => updateFromOklch('h', v)} />
                <button onClick={() => copy(`oklch(${oklch.l.toFixed(3)} ${oklch.c.toFixed(3)} ${oklch.h.toFixed(1)}deg / ${a.toFixed(2)})`, 'lch')} className="text-gray-500 hover:text-white">{copied === 'lch' ? <Check size={10} className="text-green-500"/> : <Copy size={10}/>}</button>
            </div>
            <div>
                <div className="grid grid-cols-[24px_1fr_14px] gap-1 items-center">
                    <span className="text-[9px] font-bold text-gray-600">HEX</span>
                    <HexInput value={hex} onChange={updateFromHex} label="HEX" />
                    <button onClick={() => copy(hex, 'hex')} className="text-gray-500 hover:text-white">{copied === 'hex' ? <Check size={10} className="text-green-500"/> : <Copy size={10}/>}</button>
                </div>
            </div>
            <div className="grid grid-cols-[24px_1fr_14px] gap-1 items-center mt-0.5 border-t border-gray-800 pt-1">
                <span className="text-[9px] font-bold text-gray-500">ALP</span>
                <NumberInput value={a} step={0.01} min={0} max={1} onChange={updateAlpha} />
                <div />
            </div>
        </div>
    );
}

function SplineEditor({ 
  points, 
  tangents, 
  tangentsState,  // ADD THIS
  setTangents, 
  activeIndex, 
  setActiveIndex, 
  viewIndex, 
  setViewIndex, 
  mode, 
  setMode,
  isRing,
  mathMode,
  onReset,
  hasEdits,
  width = 320, 
  height = 320,
  pathPoints 
}) {
  const svgRef = useRef(null);
  const [dragState, setDragState] = useState(null); 

  // Compute Basis for 3D -> 2D projection
  const basis = useMemo(() => {
    const toCart = (pt) => mathMode === 'oklch' ? oklchToOklab(pt) : pt;
    const toVec = (pt) => {
        const c = toCart(pt);
        return { x: c[1], y: c[2], z: c[0] };
    };

    if (points.length < 2) return { center: toVec(points[0]), u: {x:1,y:0,z:0}, v: {x:0,y:1,z:0} };
    
    // Basis computation index
    let basisIndex = viewIndex ?? 0;
    //if ring we have no end index, so we can use ends bc they properly defined via neighbors
    // if rotation we use the tangent anyway so it does not matter
    if (!isRing && mode !== 'rotation') {
        if (basisIndex === 0) basisIndex = 1;
        else if (basisIndex === points.length - 1) basisIndex = points.length - 2;
    }
    basisIndex = Math.max(0, Math.min(points.length - 1, basisIndex));

    const prevIdx = isRing ? (basisIndex - 1 + points.length) % points.length : Math.max(0, basisIndex - 1);
    const nextIdx = isRing ? (basisIndex + 1) % points.length : Math.min(points.length - 1, basisIndex + 1);
    const prev = toVec(points[prevIdx]);
    const next = toVec(points[nextIdx]);
    const basisCurr = toVec(points[basisIndex]);

    // Center on triangle centroid
    const centroid = {
        x: (prev.x + basisCurr.x + next.x) / 3,
        y: (prev.y + basisCurr.y + next.y) / 3,
        z: (prev.z + basisCurr.z + next.z) / 3
    };

    const sub = (a, b) => ({ x: a.x - b.x, y: a.y - b.y, z: a.z - b.z });
    const norm = (v) => { const l = Math.hypot(v.x, v.y, v.z); return l===0?v:{x:v.x/l, y:v.y/l, z:v.z/l}; };
    const cross = (a,b) => ({ x: a.y*b.z - a.z*b.y, y: a.z*b.x - a.x*b.z, z: a.x*b.y - a.y*b.x });

    if (mode === 'rotation' && tangents[basisIndex]) {
         let tanVec;
         if (mathMode === 'oklch') {
             const pt = points[basisIndex];
             const tipPolar = pt.map((v, i) => v + tangents[basisIndex][i]);
             const tipCart = oklchToOklab(tipPolar);
             const pCart = oklchToOklab(pt);
             const dLab = [tipCart[0]-pCart[0], tipCart[1]-pCart[1], tipCart[2]-pCart[2]];
             tanVec = { x: dLab[1], y: dLab[2], z: dLab[0] };
         } else {
             const t = tangents[basisIndex];
             tanVec = { x: t[1], y: t[2], z: t[0] };
         }

         let n = norm(tanVec);
         if (Math.hypot(n.x, n.y, n.z) < 1e-4) n = { x: 1, y: 0, z: 0 };
         
         let globalUp = { x: 0, y: 0, z: 1 };
         if (Math.abs(n.z) > 0.99) globalUp = { x: 0, y: 1, z: 0 };
         
         const u = norm(cross(n, globalUp));
         const v = norm(cross(n, u));

         return { center: centroid, u, v, n };
    }

    let t = sub(next, prev);
    let u = norm(t);
    if (Math.hypot(u.x, u.y, u.z) < 1e-4) u = { x: 1, y: 0, z: 0 };
    
    const v1 = sub(basisCurr, prev);
    const v2 = sub(next, basisCurr);
    let n = cross(v1, v2);
    if (Math.hypot(n.x, n.y, n.z) < 1e-4) {
         n = cross(u, {x:1, y:0, z:1});
         if (Math.hypot(n.x, n.y, n.z) < 1e-4) n = cross(u, {x:1, y:0, z:0});
    }
    n = norm(n);
    const v = norm(cross(n, u));
    return { center: centroid, u, v, n };
  }, [points, viewIndex, isRing, mathMode, mode, tangents]);

  const activeBasis = useMemo(() => {
    // Use frozen basis during rotation mode drag
    if (dragState && dragState.frozenBasis && mode === 'rotation') {
        return dragState.frozenBasis;
    }
    return basis;
  }, [basis, dragState, mode]);

  const scale = useMemo(() => {
    let maxDist = 0.05; 
    let basisIndex = viewIndex ?? 0;
    if (!isRing) {
        if (basisIndex === 0) basisIndex = 1;
        else if (basisIndex === points.length - 1) basisIndex = points.length - 2;
    }
    const toCart = (pt) => mathMode === 'oklch' ? oklchToOklab(pt) : pt;
    const toVec = (pt) => { const c = toCart(pt); return { x: c[1], y: c[2], z: c[0] }; };

    const pC = toVec(points[basisIndex]);
    const pP = toVec(points[isRing ? (basisIndex-1+points.length)%points.length : Math.max(0, basisIndex-1)]);
    const pN = toVec(points[isRing ? (basisIndex+1)%points.length : Math.min(points.length-1, basisIndex+1)]);
    
    const d3 = (a, b) => Math.hypot(a.x-b.x, a.y-b.y, a.z-b.z);
    maxDist = Math.max(d3(pC, pP), d3(pC, pN));
    if(maxDist === 0) maxDist = 0.1;
    const viewportRadius = Math.min(width, height);
    return (viewportRadius * 0.6) / maxDist;
  }, [points, width, height, isRing, mathMode, viewIndex]);

  const offsetX = width/2;
  const offsetY = height/2;

  const toScreen = (pt) => {
    const center = activeBasis.center; 
    const cart = mathMode === 'oklch' ? oklchToOklab(pt) : pt;
    const rx = cart[1] - center.x; 
    const ry = cart[2] - center.y; 
    const rz = cart[0] - center.z; 
    
    const x = rx * activeBasis.u.x + ry * activeBasis.u.y + rz * activeBasis.u.z;  
    const y = rx * activeBasis.v.x + ry * activeBasis.v.y + rz * activeBasis.v.z; 
    return { x: x * scale + offsetX, y: offsetY - y * scale };
  };

  const pathData = useMemo(() => {
      if (!pathPoints || pathPoints.length === 0) return "";
      return pathPoints.map((pt, i) => {
          const s = toScreen(pt);
          return `${i===0?'M':'L'} ${s.x} ${s.y}`;
      }).join(' ');
  }, [pathPoints, activeBasis, scale, mathMode]);

  const handlePointerDown = (e, index) => {
    e.stopPropagation();
    e.target.setPointerCapture(e.pointerId);
    setActiveIndex(index);
    setViewIndex && setViewIndex(index);
    
    const svgRect = svgRef.current.getBoundingClientRect();
    const tan = tangents[index] || [0,0,0,0];
    
    setDragState({ 
        index, 
        pointerId: e.pointerId, 
        startX: e.clientX - svgRect.left, 
        startY: e.clientY - svgRect.top,
        initialTangent: tan,
        handleVisVec: getVisualHandleVector(points[index], tan),
        frozenBasis: mode === 'rotation' ? basis : null  // Freeze camera in rotation mode
    });
  };

  const getVisualHandleVector = (pt, tan) => {
      const pCart = mathMode === 'oklch' ? oklchToOklab(pt) : pt;
      let hCart;
      if (mathMode === 'oklch') {
          const hPolar = pt.map((v, i) => v + tan[i]/3);
          hCart = oklchToOklab(hPolar);
      } else {
          hCart = pt.map((v, i) => v + tan[i]/3);
      }
      return { 
          x: hCart[1] - pCart[1], 
          y: hCart[2] - pCart[2], 
          z: hCart[0] - pCart[0] 
      };
  };

  const handlePointerMove = (e) => {
    if (!dragState) return;
    const svgRect = svgRef.current.getBoundingClientRect();
    const dxScreen = (e.clientX - svgRect.left - dragState.startX) / scale;
    const dyScreen = -((e.clientY - svgRect.top) - dragState.startY) / scale; 

    const dL = dxScreen * activeBasis.u.z + dyScreen * activeBasis.v.z; 
    const da = dxScreen * activeBasis.u.x + dyScreen * activeBasis.v.x; 
    const db = dxScreen * activeBasis.u.y + dyScreen * activeBasis.v.y; 
    
    let newTangent;

    if (mode === 'length') {
        const visHandle = dragState.handleVisVec;
        const dotHandle = visHandle.x*visHandle.x + visHandle.y*visHandle.y + visHandle.z*visHandle.z;
        if (dotHandle < 1e-8) {
            newTangent = dragState.initialTangent; 
        } else {
            const dotProduct = da*visHandle.x + db*visHandle.y + dL*visHandle.z;
            const multiplier = 1 + (dotProduct / dotHandle);
            newTangent = dragState.initialTangent.map(v => v * multiplier);
        }
    } else if (mode === 'rotation') {
        let rawTan;
        if (mathMode === 'oklch') {
            const pt = points[dragState.index];
            const tipPolar = pt.map((v, i) => v + dragState.initialTangent[i]/3);
            const tipCart = oklchToOklab(tipPolar);
            const newTipCart = [tipCart[0] + dL, tipCart[1] + da, tipCart[2] + db, pt[3]];
            const newTipPolar = oklabToOklch(newTipCart);
            let dh = newTipPolar[2] - pt[2];
            while (dh > 180) dh -= 360; while (dh < -180) dh += 360;
            rawTan = [(newTipPolar[0] - pt[0]) * 3, (newTipPolar[1] - pt[1]) * 3, dh * 3, 0];
        } else {
            rawTan = dragState.initialTangent.map((v, i) => {
                if(i===0) return v + dL * 3;
                if(i===1) return v + da * 3;
                if(i===2) return v + db * 3;
                return v;
            });
        }
        const initialLen = vDist3(dragState.initialTangent, [0,0,0,0]);
        const rawLen = vDist3(rawTan, [0,0,0,0]); 
        const scaleFactor = rawLen < 1e-8 ? 0 : initialLen / rawLen;
        newTangent = rawTan.map(v => v * scaleFactor);

    } else {
        if (mathMode === 'oklch') {
            const pt = points[dragState.index];
            const tipPolar = pt.map((v, i) => v + dragState.initialTangent[i]/3);
            const tipCart = oklchToOklab(tipPolar);
            const newTipCart = [tipCart[0] + dL, tipCart[1] + da, tipCart[2] + db, pt[3]];
            const newTipPolar = oklabToOklch(newTipCart);
            let dh = newTipPolar[2] - pt[2];
            while (dh > 180) dh -= 360; while (dh < -180) dh += 360;
            newTangent = [(newTipPolar[0] - pt[0]) * 3, (newTipPolar[1] - pt[1]) * 3, dh * 3, 0];
        } else {
            const deltaVec = [dL, da, db, 0];
            newTangent = dragState.initialTangent.map((v, i) => v + deltaVec[i] * 3); 
        }
    }
    
    if (newTangent) newTangent[3] = dragState.initialTangent[3];
    
    // Use sparse tangentsState
    const newTangents = tangentsState.slice();
    newTangents[dragState.index] = newTangent;
    setTangents(newTangents);
  };

  const invertTangent = () => {
        if (activeIndex === null || activeIndex === undefined) return;
        const tan = tangents[activeIndex] || [0,0,0,0];
        const newTangent = [ -tan[0], -tan[1], -tan[2], tan[3] ];
        const newTangents = tangentsState.slice();
        newTangents[activeIndex] = newTangent;
        setTangents(newTangents);
  }

  return (
    <div className="bg-gray-900 rounded-lg shadow-inner inline-block select-none relative overflow-hidden border border-gray-800 group/editor w-full h-full flex items-center justify-center">
       <div className="absolute top-2 left-3 text-[10px] text-gray-500 pointer-events-none z-10 font-mono">
         <div className="font-bold text-gray-300">{mathMode.toUpperCase()} SPACE</div>
         <div>Aligning to Stop #{ (viewIndex ?? 0) + 1 }</div>
       </div>
       
       <svg ref={svgRef} width={width} height={height} className="block cursor-default mx-auto" onPointerMove={handlePointerMove} onPointerUp={(e) => {e.target.releasePointerCapture(e.pointerId); setDragState(null);}}>
         <defs>
            <pattern id="grid" width="20" height="20" patternUnits="userSpaceOnUse"><path d="M 20 0 L 0 0 0 20" fill="none" stroke="#1f2937" strokeWidth="1"/></pattern>
         </defs>
         <rect width="100%" height="100%" fill="url(#grid)" />
         <path d={pathData} stroke="white" strokeWidth="2" fill="none" strokeLinecap="round" />
         {points.map((pt, i) => {
            const screenPt = toScreen(pt);
            const isActive = i === activeIndex;
            const isFocus = i === viewIndex;
            const tan = tangents[i];
            if (!tan) return null;
            
            let handlePt;
            if (mathMode === 'oklch') handlePt = pt.map((v, k) => v + tan[k] / 3); 
            else handlePt = pt.map((v, k) => v + tan[k] / 3); 

            const screenHandle = toScreen(handlePt);
            
            // Clamp visual handle to 40% of viewport
            const maxHandleDist = Math.min(width, height) * 0.4;
            const dx = screenHandle.x - screenPt.x;
            const dy = screenHandle.y - screenPt.y;
            const dist = Math.hypot(dx, dy);
            
            let clampedHandle = screenHandle;
            if (dist > maxHandleDist) {
                const scale = maxHandleDist / dist;
                clampedHandle = {
                    x: screenPt.x + dx * scale,
                    y: screenPt.y + dy * scale
                };
            }
            
            const isDragging = dragState && dragState.index === i;
            return (
              <g key={i} opacity={isFocus ? 1 : (isActive ? 0.8 : 0.4)}>
                <line x1={screenPt.x} y1={screenPt.y} x2={clampedHandle.x} y2={clampedHandle.y} stroke="#4b5563" pointerEvents="none" />
                <circle cx={screenPt.x} cy={screenPt.y} r={4} fill={isActive ? "#fff" : "#9ca3af"} 
                    onPointerDown={(e) => { e.stopPropagation(); setActiveIndex(i); setViewIndex && setViewIndex(i); }}
                    className="cursor-pointer"
                />
                <g className={mode === 'free' ? 'cursor-move' : 'cursor-pointer'} onPointerDown={(e) => handlePointerDown(e, i)} opacity= {isActive ? 1 : 0.4}>
                    <circle cx={clampedHandle.x} cy={clampedHandle.y} r={12} fill="transparent" />
                    <circle cx={clampedHandle.x} cy={clampedHandle.y} r={5} fill={isDragging ? "#ffffff" : (mode === 'free' ? "#ef4444" : mode === 'rotation' ? "#10b981" : "#3b82f6")} stroke="black" strokeWidth={1}/>
                </g>
              </g>
            );
        })}
       </svg>
       
      <div className="absolute top-2 right-2 z-20 flex flex-col gap-1 items-end">
        <div className="flex gap-1">
            <div className=" bg-gray-950 rounded border border-gray-800 p-0.5">
                <Tooltip className="px-1" text={`Invert Tangent: Inverts Tangent Direction`} side="bottom">
                    <button onClick={() => invertTangent()} className={'p-1.5 rounded text-gray-400 hover:text-gray-200 active:text-white active:bg-purple-400'}><ArrowLeftRight size={14} /></button>
                </Tooltip>
            </div>
            <div className="flex gap-1 bg-gray-950 rounded border border-gray-800 p-0.5">
                <Tooltip text={`Free Mode: Adjust length and direction`} side="bottom">
                    <button onClick={() => setMode('free')} className={`p-1.5 rounded ${mode === 'free' ? 'bg-red-600 text-white' : 'text-gray-400 hover:text-gray-200'}`}><Unlock size={14} /></button>
                </Tooltip>
                <Tooltip text="Length Mode: Adjust tension only (drag along line)" side="bottom">
                    <button onClick={() => setMode('length')} className={`p-1.5 rounded ${mode === 'length' ? 'bg-blue-600 text-white' : 'text-gray-400 hover:text-gray-200'}`}><Ruler size={14} /></button>
                </Tooltip>
                <Tooltip text="Rotation Mode: Adjust angle only (preserve length)" side="bottom">
                    <button onClick={() => setMode('rotation')} className={`p-1.5 rounded ${mode === 'rotation' ? 'bg-green-600 text-white' : 'text-gray-400 hover:text-gray-200'}`}>
                        <Compass size={14} />
                    </button>
                </Tooltip>
            </div>
        </div>
        {hasEdits && (
              <div className="bg-gray-950 rounded border border-gray-800 p-0.5">
                <div className="w-[1px] bg-gray-800 mx-0.5" />
                <Tooltip text="Reset Spatial Tangents" side="bottom">
                    <button onClick={onReset} className="p-1.5 rounded text-red-400 hover:text-red-200 hover:bg-gray-800"><RotateCcw size={14} /></button>
                </Tooltip>
              </div>
          )}
      </div>
    </div>
  );
}

function AlphaEditor({ 
  points, 
  tangents, 
  tangentsState,  // ADD THIS
  setTangents, 
  splineData,
  isRing,  // ADD THIS
  activeIndex, 
  setActiveIndex, 
  onReset,
  hasEdits,
  width = 320, 
  height = 100 
}) {
    const svgRef = useRef(null);
    const [dragState, setDragState] = useState(null);

    const uValues = useMemo(() => {
        if (splineData && splineData.u) return splineData.u;
        return points.map((_, i) => i / (points.length - 1));
    }, [points, splineData]);
    
    const totalLen = uValues[uValues.length - 1] || 1;

    const { extendedWidth, xOffset } = useMemo(() => {
      if (!isRing || points.length < 2) return { extendedWidth: width, xOffset: 0 };

      const lastSegLen = uValues[uValues.length - 1] - uValues[uValues.length - 2];
      const leftExt = (lastSegLen / totalLen) * width * 0.5;
      return { 
          extendedWidth: width, 
          xOffset: leftExt 
      };
    }, [isRing, points.length, uValues, totalLen, width]);

    // Auto-scale height
    const { actualHeight, yMin, yMax } = useMemo(() => {
      let min = Infinity;
      let max = -Infinity;
      
      points.forEach((pt, i) => {
          const alpha = pt[3] ?? 1;
          min = Math.min(min, alpha);
          max = Math.max(max, alpha);
          
          const tan = tangents[i];
          if (tan) {
              let segLen = 1;
              if (i < points.length - 1) segLen = uValues[i+1] - uValues[i];
              else if (i > 0) segLen = uValues[i] - uValues[i-1];
              
              const handleAlpha = alpha + (tan[3] * segLen) / 3;
              min = Math.min(min, handleAlpha);
              max = Math.max(max, handleAlpha);
          }
      });
      
      const range = max - min;
      const padding = range * 0.1;
      min -= padding;
      max += padding;
      
      min = Math.min(min, 0);
      max = Math.max(max, 1);
      
      const minHeight = 100;
      const calculatedHeight = Math.max(minHeight, height);
      
      return { actualHeight: calculatedHeight, yMin: min, yMax: max };
    }, [points, tangents, uValues, height]);

    const toScreen = (t, alpha) => {
      const normX = t / totalLen;
      const normY = (alpha - yMin) / (yMax - yMin);
      const padding = 20;
      const availW = width - padding * 2;
      const availH = actualHeight - padding * 2;
      return {
          x: padding + xOffset + normX * availW,
          y: actualHeight - padding - normY * availH 
      };
    };

    const pathData = useMemo(() => {
      if (!splineData) return "";
      const pathOps = [];
      const numSegments = isRing ? points.length - 1 : splineData.segs.length; // Don't draw wraparound in main loop
      
      for (let i = 0; i < numSegments; i++) {
        const p0 = points[i];
        const p1 = points[i + 1];
        let t0 = tangents[i];
        let t1 = tangents[i + 1];
        const u0 = uValues[i];
        const u1 = uValues[i + 1];
        const segLen = u1 - u0;
        
        // Apply the same rescaling as in fitSpline for non-natural-cubic methods
        if (splineData && splineData.method !== 'natural-cubic') {
            const dtCurrent = segLen;
            const dtNext = i < numSegments - 1 ? (uValues[i + 2] - uValues[i + 1]) : segLen;
            const scale = dtNext > 1e-9 ? (dtCurrent / dtNext) : 1;
            t1 = t1.map(v => v * scale);
        }
        
        const y0 = p0[3] ?? 1;
        const y1 = p1[3] ?? 1;
        const cp0y = y0 + (t0[3] * segLen) / 3;
        const cp1y = y1 - (t1[3] * segLen) / 3;
        
        const start = toScreen(u0, y0);
        const end = toScreen(u1, y1);
        const cp1 = toScreen(u0 + segLen / 3, cp0y);
        const cp2 = toScreen(u1 - segLen / 3, cp1y);
        
        if (i === 0) pathOps.push(`M ${start.x} ${start.y}`);
        pathOps.push(`C ${cp1.x} ${cp1.y}, ${cp2.x} ${cp2.y}, ${end.x} ${end.y}`);
      }
      
      // Add wraparound visual for loop mode
      if (isRing && points.length > 1) {
        const lastIdx = points.length - 1;
        const p0 = points[lastIdx];
        const p1 = points[0];
        let t0 = tangents[lastIdx];
        let t1 = tangents[0];
        
        const lastSegLen = totalLen - uValues[lastIdx];
        
        // Apply the same rescaling as in fitSpline for non-natural-cubic methods
        if (splineData && splineData.method !== 'natural-cubic') {
            const dtCurrent = lastSegLen; // chord[n-1]
            const dtNext = uValues[1] - uValues[0]; // chord[0]
            const scale = dtNext > 1e-9 ? (dtCurrent / dtNext) : 1;
            t1 = t1.map(v => v * scale);
        }
        
        const y0 = p0[3] ?? 1;
        const y1 = p1[3] ?? 1;
        const cp0y = y0 + (t0[3] * lastSegLen) / 3;
        const cp1y = y1 - (t1[3] * lastSegLen) / 3;
        
        // Right half: last point to middle of wraparound
        const uMid = totalLen - lastSegLen * 0.5;
        const t = 0.5;
        const yMid = (1-t)**3 * y0 + 3*(1-t)**2*t * cp0y + 3*(1-t)*t**2 * cp1y + t**3 * y1;

        const uStart_r = totalLen - lastSegLen;
        const uEnd = uMid;
        const cp1U = uStart_r + lastSegLen / 3;
        const cp2U = uStart_r + lastSegLen / 3 + lastSegLen / 6;

        const startRight = toScreen(uStart_r, y0);
        const endRight = toScreen(uEnd, yMid);
        const cp1Right = toScreen(cp1U, cp0y);
        const cp2Right = toScreen(cp2U, (cp0y + cp1y) * 0.5);

        pathOps.push(`M ${startRight.x} ${startRight.y}`);
        pathOps.push(`C ${cp1Right.x} ${cp1Right.y}, ${cp2Right.x} ${cp2Right.y}, ${endRight.x} ${endRight.y}`);
        
        // Left half: middle of wraparound to first point
        const uStart = -lastSegLen * 0.5;
        
        const startLeft = toScreen(uStart, yMid);
        const endLeft = toScreen(0, y1);
        const cp1Left = toScreen(uStart + lastSegLen / 3, cp1y);
        const cp2Left = toScreen(-lastSegLen / 6, y1 + (cp1y - y1) * 0.5);
        
        pathOps.push(`M ${startLeft.x} ${startLeft.y}`);
        pathOps.push(`C ${cp1Left.x} ${cp1Left.y}, ${cp2Left.x} ${cp2Left.y}, ${endLeft.x} ${endLeft.y}`);
      }
      
      return pathOps.join(" ");
    }, [points, tangents, uValues, totalLen, width, actualHeight, splineData, yMin, yMax, isRing]);

    const handlePointerDown = (e, index) => {
        e.stopPropagation();
        e.target.setPointerCapture(e.pointerId);
        setActiveIndex(index);
        const svgRect = svgRef.current.getBoundingClientRect();
        setDragState({
            index,
            pointerId: e.pointerId,
            startY: e.clientY - svgRect.top,
            initialAlphaTan: tangents[index][3]
        });
    };

    const handlePointerMove = (e) => {
      if (!dragState) return;
      const svgRect = svgRef.current.getBoundingClientRect();
      const dyScreen = (e.clientY - svgRect.top) - dragState.startY;
      const padding = 20;
      const availH = actualHeight - padding * 2;
      const dAlpha = -(dyScreen / availH) * (yMax - yMin);

      let segLen = 1;
      if (points.length > 1) {
          const idx = dragState.index;
          if (idx < points.length - 1) segLen = uValues[idx+1] - uValues[idx];
          else segLen = uValues[idx] - uValues[idx-1];
      }
      
      const newTanAlpha = dragState.initialAlphaTan + (dAlpha * 3) / segLen;
      
      // Use sparse tangentsState
      const newTangents = tangentsState.slice();
      const currentTangent = tangents[dragState.index];
      newTangents[dragState.index] = [
          currentTangent[0],
          currentTangent[1],
          currentTangent[2],
          newTanAlpha
      ];
      setTangents(newTangents);
    };

    return (
        <div className="bg-gray-900 rounded-lg shadow-inner relative overflow-hidden border border-gray-800 w-full h-full flex items-center justify-center group/alpha">
            <div className="absolute top-2 left-3 text-[10px] text-gray-500 font-mono font-bold pointer-events-none z-10">
                ALPHA CURVE
            </div>
             <div className="absolute top-3 right-3 z-20 flex bg-gray-950 rounded border border-gray-800 p-0.5">
               {hasEdits && (
                    <Tooltip text="Reset Alpha Tangents" side="bottom">
                        <button onClick={onReset} className="p-1.5 rounded text-red-400 hover:text-red-200 hover:bg-gray-800"><RotateCcw size={14} /></button>
                    </Tooltip>
               )}
            </div>
            
            <svg 
              ref={svgRef} 
              width={extendedWidth} 
              height={actualHeight}  
              className="block cursor-default"
              onPointerMove={handlePointerMove}
              onPointerUp={(e) => { e.target.releasePointerCapture(e.pointerId); setDragState(null); }}
            >
                 <defs>
                    <pattern id="gridAlpha" width="20" height="20" patternUnits="userSpaceOnUse"><path d="M 20 0 L 0 0 0 20" fill="none" stroke="#1f2937" strokeWidth="1"/></pattern>
                 </defs>
                 <rect width="100%" height="100%" fill="url(#gridAlpha)" />
                 <line x1="0" y1={toScreen(0,0).y} x2={width} y2={toScreen(0,0).y} stroke="#374151" strokeDasharray="2 2" />
                 <line x1="0" y1={toScreen(0,1).y} x2={width} y2={toScreen(0,1).y} stroke="#374151" strokeDasharray="2 2" />
                 <path d={pathData} stroke="#60a5fa" strokeWidth="2" fill="none" />
                 {points.map((pt, i) => {
                    const u = uValues[i];
                    const alpha = pt[3] ?? 1;  // Direct alpha value
                    const tan = tangents[i];
                    if(!tan) return null;
                    
                    const screenPt = toScreen(u, alpha);
                    let segLen = 1;
                    if (i < points.length - 1) segLen = uValues[i+1] - uValues[i];
                    else if (i > 0) segLen = uValues[i] - uValues[i-1];

                    const handleU = u + segLen / 3;
                    const handleAlpha = alpha + (tan[3] * segLen) / 3;
                    const screenHandle = toScreen(handleU, handleAlpha);
                    const isActive = i === activeIndex;
                    const isDragging = dragState && dragState.index === i;

                    return (
                      <g key={i}>
                          <line x1={screenPt.x} y1={screenPt.y} x2={screenPt.x} y2={screenHandle.y} stroke="#4b5563" strokeWidth="1" strokeDasharray="2 2" />
                          <circle cx={screenPt.x} cy={screenPt.y} r={3} fill={isActive ? "#fff" : "#6b7280"} />
                          <line x1={screenPt.x} y1={screenPt.y} x2={screenPt.x} y2={screenHandle.y} stroke="#9ca3af" />
                          <circle 
                            cx={screenPt.x} 
                            cy={screenHandle.y} 
                            r={isDragging ? 5 : 4} 
                            fill={isDragging ? "#fff" : "#3b82f6"} 
                            className="cursor-ns-resize"
                            onPointerDown={(e) => handlePointerDown(e, i)}
                          />
                      </g>
                     );
                 })}
            </svg>
        </div>
    );
}

function ColorPicker({ color, onChange, isActive, onActivate, onRemove, canRemove, onDragStart, onDragEnter, onDrop }) {
    const [h, s, b, a = 1] = color;
    const hueRef = useRef(null);
    const sbRef = useRef(null);
    const alphaRef = useRef(null);
  
    const handleHue = (e) => {
      const rect = hueRef.current.getBoundingClientRect();
      const y = Math.max(0, Math.min(1, (e.clientY - rect.top) / rect.height));
      onChange([y * 360, s, b, a]);
    };
    const handleSB = (e) => {
      const rect = sbRef.current.getBoundingClientRect();
      const x = Math.max(0, Math.min(1, (e.clientX - rect.left) / rect.width));
      const y = Math.max(0, Math.min(1, (e.clientY - rect.top) / rect.height));
      onChange([h, x, 1 - y, a]);
    };
    const handleAlpha = (e) => {
      const rect = alphaRef.current.getBoundingClientRect();
      const y = Math.max(0, Math.min(1, (e.clientY - rect.top) / rect.height));
      onChange([h, s, b, 1 - y]); 
    };
    
    return (
      <div 
          className={`flex flex-col gap-2 transition-all p-1.5 rounded-lg border ${isActive ? 'bg-gray-800/50 border-blue-500/30' : 'hover:bg-gray-800/30 border-gray-800'}`} 
          onClick={onActivate}
          draggable
          onDragStart={onDragStart}
          onDragEnter={onDragEnter}
          onDragEnd={onDrop}
          onDragOver={(e) => e.preventDefault()}
      >
          <div className="flex justify-between items-center px-1 mb-1">
               <div className="cursor-grab text-gray-600 hover:text-gray-400 active:cursor-grabbing">
                  <GripVertical size={12} />
               </div>
               {canRemove && (
                   <button onClick={(e) => { e.stopPropagation(); onRemove(); }} className="text-gray-600 hover:text-red-400 p-0.5">
                      <X size={12} />
                   </button>
               )}
          </div>
          
          <div className="flex gap-2 h-36" onDragStart={(e) => { e.preventDefault(); e.stopPropagation(); }} draggable>
              <div ref={hueRef} className="relative w-4 h-full rounded-full overflow-hidden cursor-pointer shadow-sm touch-none ring-1 ring-white/10 flex-none" style={{ background: 'linear-gradient(to bottom, #f00 0%, #ff0 17%, #0f0 33%, #0ff 50%, #00f 67%, #f0f 83%, #f00 100%)' }}
                  onPointerDown={(e) => { e.target.setPointerCapture(e.pointerId); onActivate(); handleHue(e); }}
                  onPointerMove={(e) => e.buttons === 1 && handleHue(e)}
                  onPointerUp={(e) => e.target.releasePointerCapture(e.pointerId)}>
                  <div className="absolute left-0 right-0 h-1 border border-white bg-black/20 rounded-full" style={{ top: `${h / 360 * 100}%`, transform: 'translateY(-50%)' }} />
              </div>
              <div ref={sbRef} className="relative flex-1 h-full rounded-md overflow-hidden cursor-crosshair shadow-sm touch-none ring-1 ring-white/10" style={{ background: `linear-gradient(to top, black, transparent), linear-gradient(to right, white, hsl(${h}, 100%, 50%))` }}
                  onPointerDown={(e) => { e.target.setPointerCapture(e.pointerId); onActivate(); handleSB(e); }}
                  onPointerMove={(e) => e.buttons === 1 && handleSB(e)}
                  onPointerUp={(e) => e.target.releasePointerCapture(e.pointerId)}>
                  <div className="absolute w-3 h-3 border-2 border-white rounded-full shadow-md" style={{ left: `${s * 100}%`, top: `${(1 - b) * 100}%`, transform: 'translate(-50%, -50%)' }} />
              </div>
              <div ref={alphaRef} className="relative w-4 h-full rounded-full overflow-hidden cursor-pointer shadow-sm touch-none ring-1 ring-white/10 flex-none bg-[url('data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSI4IiBoZWlnaHQ9IjgiPjxwYXRoIGQ9Ik0wIDBoNHY0SDB6bTQgNGg0djRINHoiIGZpbGw9IiMzMzMiLz48L3N2Zz4=')] bg-repeat"
                  onPointerDown={(e) => { e.target.setPointerCapture(e.pointerId); onActivate(); handleAlpha(e); }}
                  onPointerMove={(e) => e.buttons === 1 && handleAlpha(e)}
                  onPointerUp={(e) => e.target.releasePointerCapture(e.pointerId)}>
                  <div className="absolute inset-0" style={{ background: `linear-gradient(to top, transparent, hsl(${h}, 100%, 50%))` }} />
                  <div className="absolute left-0 right-0 h-1 border border-white bg-black/50 rounded-full" style={{ top: `${(1 - a) * 100}%`, transform: 'translateY(-50%)' }} />
              </div>
          </div>
          <ColorInputs color={color} onChange={onChange} />
      </div>
    );
  }


function CodeBlock({ label, value }) {
    const [copied, setCopied] = useState(false);
    const handleCopy = () => {
        navigator.clipboard.writeText(value);
        setCopied(true);
        setTimeout(() => setCopied(false), 2000);
    };
    return (
        <div className="bg-gray-900 rounded-lg border border-gray-800 overflow-hidden flex flex-col h-full">
            <div className="flex justify-between items-center px-3 py-2 bg-gray-800/50 border-b border-gray-800 whitespace-nowrap">
                <span className="text-[10px] font-bold uppercase tracking-widest text-gray-400 flex items-center gap-2 flex-1 overflow-hidden text-ellipsis">
                    {label.includes('Strength') ? <Braces size={12}/> : <FileCode size={12}/>} 
                    {label}
                </span>
                <Tooltip text="Copy content" side="left">
                    <button onClick={handleCopy} className="text-gray-400 hover:text-white transition-colors flex-none ml-2">
                        {copied ? <Check size={14} className="text-green-400"/> : <Copy size={14}/>}
                    </button>
                </Tooltip>
            </div>
            <textarea readOnly value={value} className="w-full h-24 bg-transparent p-3 text-[10px] font-mono text-gray-400 resize-none focus:outline-none focus:text-gray-200" onClick={(e) => e.target.select()}/>
        </div>
    );
}

// ==========================================
// 3. Main Application
// ==========================================

export default function App() {
  const [colors, setColors] = useState([ [0, 1, 1, 1], [120, 1, 1, 1], [240, 1, 1, 1] ]);
  const [activePointIndex, setActivePointIndex] = useState(1);
  const [viewPointIndex, setViewPointIndex] = useState(1);
  const [position, setPosition] = useState(0.5); 
  const [steps, setSteps] = useState(40); 
  const [editorMode, setEditorMode] = useState('length')
  const [copied, setCopied] = useState(false);
  const [isRing, setIsRing] = useState(false);
  const [mathMode, setMathMode] = useState('oklab'); 
  const [splineMethod, setSplineMethod] = useState('natural-cubic');
  const [distribution, setDistribution] = useState('uniform'); 
  const [tangentsState, setTangentsState] = useState([]); 

  const dragItem = useRef(null);
  const dragOverItem = useRef(null);
  const handleDragStart = (e, index) => { dragItem.current = index; };
  const handleDragEnter = (e, index) => { dragOverItem.current = index; };
  const handleDrop = (e) => {
      const copyListItems = [...colors];
      const dragItemContent = copyListItems[dragItem.current];
      copyListItems.splice(dragItem.current, 1);
      copyListItems.splice(dragOverItem.current, 0, dragItemContent);
      dragItem.current = null;
      dragOverItem.current = null;
      setColors(copyListItems);
      setActivePointIndex(0);
      setViewPointIndex(0);
      setTangentsState([]);
  };

  const minSteps = Math.max(2, colors.length + (isRing ? 1 : 0));
  useEffect(() => { if(steps < minSteps) setSteps(minSteps); }, [minSteps, steps]);

  // Tangent Reset Logic: Ensures manual edits are discarded when the fundamental geometry changes
  useEffect(() => {
     setTangentsState([]);
  }, [colors.length, mathMode, isRing, splineMethod, distribution]);


  // 1. Prepare Points (Correct Space)
  const fixpoints = useMemo(() => {
      return colors.map(c => {
           const rgb = hsbToRgb(c[0], c[1], c[2]); 
           const lab = rgbToOklab([...rgb, c[3] ?? 1]);
           if (mathMode === 'oklch') return oklabToOklch(lab);
           return lab; 
      });
  }, [colors, mathMode]);

  // 2. Generate Spline Data
  const splineData = useMemo(() => {
      try {
          return fitSpline(fixpoints, { 
              space: mathMode, 
              method: splineMethod, 
              loop: isRing, 
              strengths: 1,
              overrideTangents: tangentsState 
          });
      } catch (e) {
          console.error(e);
          return null;
      }
  }, [fixpoints, mathMode, splineMethod, isRing, tangentsState]);

  // 3. Generate Colors
  const gradientStops = useMemo(() => {
      if (!splineData) return [];
      const weights = Array.from({length: steps + 1}, (_, i) => i / steps);
      const cols = splineColors(fixpoints, weights, {
          space: mathMode,
          method: splineMethod,
          loop: isRing,
          strengths: 1,
          distribution,
          overrideTangents: tangentsState,
          format: 'oklch'
      });
      return cols.map((c, i) => {
          const w = weights[i] * 100;
          return `oklch(${c.l.toFixed(4)} ${c.c.toFixed(4)} ${c.h.toFixed(2)}deg / ${c.a?.toFixed(3) ?? 1}) ${w.toFixed(1)}%`;
      });
  }, [fixpoints, steps, mathMode, splineMethod, isRing, distribution, splineData, tangentsState]);

  const polyGradient = `linear-gradient(to right, ${gradientStops.join(', ')})`;

  const { polyVal, polyRgb } = useMemo(() => {
      if(!splineData) return { polyVal: {l:0,c:0,h:0,a:1}, polyRgb: {r:0,g:0,b:0} };
      const col = splineColors(fixpoints, [position], {
          space: mathMode,
          method: splineMethod,
          loop: isRing,
          distribution,
          overrideTangents: tangentsState,
          format: 'oklch' 
      })[0];
      
      const alpha = Math.max(Math.min(col.a ?? 1, 1), 0);
      const rgb = oklchToRgb([col.l, col.c, col.h, alpha]).map(v => Math.max(0, Math.min(255, Math.round(v * 255))));
      
      return { 
          polyVal: { l: col.l, c: col.c, h: col.h, a: alpha }, 
          polyRgb: { r: rgb[0], g: rgb[1], b: rgb[2] }
      };
  }, [fixpoints, position, mathMode, splineMethod, isRing, distribution, splineData, tangentsState]);

  // 5. Linear Reference (Manual Interp)
  const { linearGradient, linearVal, linearRgb } = useMemo(() => {
      const stops = colors.map((c, i, arr) => {
           const rgb = hsbToRgb(c[0], c[1], c[2]);
           const lab = rgbToOklab([...rgb, c[3]??1]);
           const lch = oklabToOklch(lab);
           return `oklch(${lch[0].toFixed(4)} ${lch[1].toFixed(4)} ${lch[2].toFixed(2)}deg / ${lch[3]?.toFixed(3)})`;
      });
      if(isRing) stops.push(stops[0]); 
      const distributed = stops.map((s, i) => `${s} ${(i/(stops.length-1)*100).toFixed(1)}%`);
      const linearGradient = `linear-gradient(to right, ${distributed.join(', ')})`;

      const segments = isRing ? colors.length : colors.length - 1;
      const scaledP = position * segments;
      let idx0 = Math.floor(scaledP);
      if (idx0 >= segments) idx0 = segments - 1;
      const t = scaledP - idx0;
      const idx1 = (idx0 + 1) % colors.length;

      const c0 = fixpoints[idx0];
      const c1 = fixpoints[idx1 % fixpoints.length];

      let finalArr;
      if (mathMode === 'oklch') {
          // c0/c1 are OKLCH arrays
          let dh = c1[2] - c0[2];
          while (dh > 180) dh -= 360; while (dh < -180) dh += 360;
          finalArr = [
              c0[0] + (c1[0] - c0[0]) * t,
              c0[1] + (c1[1] - c0[1]) * t,
              c0[2] + dh * t,
              (c0[3]??1) + ((c1[3]??1) - (c0[3]??1)) * t
          ];
          finalArr[2] = (finalArr[2] % 360 + 360) % 360;
      } else {
          // c0/c1 are OKLAB arrays
          finalArr = [
              c0[0] + (c1[0] - c0[0]) * t,
              c0[1] + (c1[1] - c0[1]) * t,
              c0[2] + (c1[2] - c0[2]) * t,
              (c0[3]??1) + ((c1[3]??1) - (c0[3]??1)) * t
          ];
      }

      const linLCH = mathMode === 'oklch' 
          ? { l: finalArr[0], c: finalArr[1], h: finalArr[2], a: finalArr[3] }
          : formatColor(finalArr, 'oklch', 'oklab');
      
      const linRGBArr = oklchToRgb([linLCH.l, linLCH.c, linLCH.h, linLCH.a]);
      
      const clamp = (v) => Math.max(0, Math.min(255, Math.round(v * 255)));
      const linRGB = { r: clamp(linRGBArr[0]), g: clamp(linRGBArr[1]), b: clamp(linRGBArr[2]) };

      return { linearGradient, linearVal: linLCH, linearRgb: linRGB };
  }, [colors, isRing, fixpoints, position, mathMode]);

  // 6. Tangent Strength (Directional Logic)
  const tangentStrengthList = useMemo(() => {
      const defaultData = fitSpline(fixpoints, { space: mathMode, method: splineMethod, loop: isRing, strengths: 1, overrideTangents: [] });
      const currentTangents = splineData ? splineData.tangents : defaultData.tangents;
      
      const strengths = currentTangents.map((t, i) => {
          const defTan = defaultData.tangents[i];
          const defLen = vLen(defTan);
          const curLen = vLen(t);
          
          if (defLen < 1e-6) return 1.0;
          
          // Calculate sign based on dot product direction
          const dot = vDot(t, defTan);
          const sign = dot < 0 ? -1 : 1;
          
          return parseFloat((sign * (curLen / defLen)).toFixed(2));
      });
      return JSON.stringify(strengths, null, 0).replace(/,/g, ', ');
  }, [fixpoints, mathMode, splineMethod, isRing, splineData]);

  // 7. Editor Data
  const editorPoints = useMemo(() => fixpoints, [fixpoints]);

  const editorPathPoints = useMemo(() => {
      const res = 100;
      const weights = Array.from({length: res + 1}, (_, i) => i / res);
      const raw = splineColors(fixpoints, weights, {
          space: mathMode,
          method: splineMethod,
          loop: isRing,
          distribution,
          overrideTangents: tangentsState,
          format: mathMode 
      });
      if (mathMode === 'oklch') {
          return raw.map(o => [o.l, o.c, o.h, o.a ?? 1]);
      }
      return raw;
  }, [fixpoints, mathMode, splineMethod, isRing, distribution, tangentsState]);

  const displayTangents = useMemo(() => {
     if(!splineData) return [];
     return splineData.tangents;
  }, [splineData]);

  const defaultSplineData = useMemo(() => {
      try {
          return fitSpline(fixpoints, { 
              space: mathMode, 
              method: splineMethod, 
              loop: isRing, 
              strengths: 1, 
              overrideTangents: [] // Force empty overrides
          });
      } catch (e) { return null; }
  }, [fixpoints, mathMode, splineMethod, isRing]);

  // 2. Reset Handlers
  const handleResetSpatial = () => {
      if (!defaultSplineData) return;
      // Revert indices 0,1,2 to default, keep 3 (Alpha) from current state
      const defaults = defaultSplineData.tangents;
      const mixed = defaults.map((def, i) => {
          const current = tangentsState[i] || def;
          return [def[0], def[1], def[2], current[3] ?? def[3]];
      });
      setTangentsState(mixed);
  };

  const handleResetAlpha = () => {
      if (!defaultSplineData) return;
      // Keep indices 0,1,2 from current state, revert 3 to default
      const defaults = defaultSplineData.tangents;
      const mixed = defaults.map((def, i) => {
          const current = tangentsState[i] || def;
          return [current[0], current[1], current[2], def[3]];
      });
      setTangentsState(mixed);
  };

  // Check if we have edits to show buttons
  const hasSpatialEdits = useMemo(() => {
      if (!defaultSplineData || !tangentsState.length) return false;
      return tangentsState.some((t, i) => {
         if (!t) return false;
         const d = defaultSplineData.tangents[i];
         // Check dist for first 3 dims
         return vDist3(t, d) > 1e-4;
      });
  }, [tangentsState, defaultSplineData]);

  const hasAlphaEdits = useMemo(() => {
      if (!defaultSplineData || !tangentsState.length) return false;
      return tangentsState.some((t, i) => {
         if (!t) return false;
         const d = defaultSplineData.tangents[i];
         return Math.abs(t[3] - d[3]) > 1e-4;
      });
  }, [tangentsState, defaultSplineData]);

  const copyToClipboard = () => {
    navigator.clipboard.writeText(polyGradient);
    setCopied(true);
    setTimeout(() => setCopied(false), 2000);
  };

  return (
    <div className="min-h-screen bg-gray-950 text-gray-200 py-8 px-4 font-sans flex justify-center">
      <div className="w-full max-w-6xl flex flex-col gap-6">
        
        <header className="flex flex-col md:flex-row justify-between items-start md:items-end border-b border-gray-800 pb-4 gap-4">
          <div>
            <h1 className="text-2xl font-bold text-white mb-1 tracking-tight flex items-center gap-2">
               <GitCommit className="text-blue-500" />
               {mathMode ==="oklab" ? "OKLAB" : "OKLCH"} Spline Gradient
            </h1>
            <p className="text-sm text-gray-400">Generates smooth interpolations using various spline methods.</p>
          </div>
          <div className="flex gap-2 items-center flex-wrap md:flex-nowrap">
            
            <div className="flex bg-gray-800 rounded p-1 gap-1">
                <Tooltip text="Cartesian Space (Physically Accurate)" side="bottom">
                    <button onClick={() => setMathMode('oklab')} className={`p-1.5 rounded flex items-center gap-2 text-xs font-bold ${mathMode === 'oklab' ? 'bg-blue-600 text-white' : 'text-gray-400 hover:text-gray-200'}`}>
                        <Box size={14}/> OKLAB
                    </button>
                </Tooltip>
                <Tooltip text="Polar Space (Preserves Saturation, Rainbows)" side="bottom">
                    <button onClick={() => setMathMode('oklch')} className={`p-1.5 rounded flex items-center gap-2 text-xs font-bold ${mathMode === 'oklch' ? 'bg-pink-600 text-white' : 'text-gray-400 hover:text-gray-200'}`}>
                        <Globe size={14}/> OKLCH
                    </button>
                </Tooltip>
            </div>

             <div className="relative group z-50">
                <div className="flex items-center gap-2 bg-gray-800 text-xs font-medium px-3 py-1.5 rounded text-gray-300 border border-gray-700">
                    {distribution === 'uniform' ? <AlignJustify size={14}/> : <AlignCenter size={14}/>}
                    <select 
                        value={distribution} 
                        onChange={(e) => setDistribution(e.target.value)}
                        className="bg-gray-800 text-gray-200 appearance-none outline-none cursor-pointer pr-4"
                    >
                        <option value="uniform">Uniform</option>
                        <option value="geometric">Geometric</option>
                    </select>
                </div>
            </div>

             <div className="relative group z-50">
                <div className="flex items-center gap-2 bg-gray-800 text-xs font-medium px-3 py-1.5 rounded text-gray-300 border border-gray-700">
                    <Settings2 size={14}/>
                    <select 
                        value={splineMethod} 
                        onChange={(e) => setSplineMethod(e.target.value)}
                        className="bg-gray-800 text-gray-200 appearance-none outline-none cursor-pointer pr-4"
                    >
                        <option value="chordal-CR">Chordal Catmull-Rom</option>
                        <option value="centripetal-CR">Centripetal Catmull-Rom</option>
                        <option value="natural-cubic">Natural Cubic</option>
                    </select>
                </div>
            </div>

            {tangentsState.some(t => t) && (
                <Tooltip text="Reset manual edits" side="bottom">
                    <button onClick={() => setTangentsState([])} className="flex items-center gap-2 text-xs font-medium bg-gray-800 hover:bg-gray-700 px-3 py-1.5 rounded text-blue-300 transition-colors whitespace-nowrap">
                        <RotateCcw size={14} /> Reset
                    </button>
                </Tooltip>
            )}

            <Tooltip text="Copy the Spline Gradient CSS" side="bottom">
                <button onClick={copyToClipboard} className="flex items-center gap-2 text-xs font-medium bg-blue-600 hover:bg-blue-500 px-3 py-1.5 rounded text-white transition-colors whitespace-nowrap">
                    {copied ? <Check size={14} /> : <Copy size={14} />} 
                    {copied ? 'Copied' : 'Copy CSS'}
                </button>
            </Tooltip>
          </div>
        </header>

        <div className="grid grid-cols-1 lg:grid-cols-2 gap-8">
            <div className="space-y-6">
                <div className="space-y-1">
                    <div className="flex justify-between items-center text-[10px] font-bold uppercase tracking-widest text-gray-500 px-1">
                        <span className="flex items-center gap-2">
                             Spline ({mathMode.toUpperCase()} / {splineMethod} / {distribution})
                        </span>
                        <span className="flex items-center gap-2">
                             <Tooltip text="Gradient Resolution (Steps)" side="right">
                                <input type="range" min={minSteps} max="100" value={steps} onChange={(e) => setSteps(Number(e.target.value))} className="w-16 accent-blue-500 cursor-pointer"/>
                             </Tooltip>
                             <span className="text-gray-600 font-mono normal-case">{steps} stops</span>
                        </span>
                    </div>
                    
                    <div className="relative group select-none rounded-lg shadow-lg ring-1 ring-white/10 overflow-hidden">
                        <div className="absolute inset-0 bg-[url('data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgaGVpZ2h0PSIxNiI+PHBhdGggZD0iTTAgMGg4djhIMHptOCA4aDh2OEg4eiIgZmlsbD0iIzIyMiIvPjwvc3ZnPg==')] z-0" />
                        <div className="h-16 w-full z-10 relative" style={{ background: polyGradient }} />
                        <div className="h-16 w-full z-10 relative border-t-2 border-black/100 " style={{ background: linearGradient }} />
                        <div className="absolute top-0 bottom-0 z-20 pointer-events-none flex flex-col items-center justify-center w-full" style={{ left: `${position * 100}%`, transform: 'translateX(-50%)' }}>
                            <div className="w-0.5 h-full bg-white shadow-[0_0_10px_rgba(0,0,0,0.5)]" />
                            <div className="absolute top-1/2 -translate-y-1/2 bg-white text-gray-950 text-[10px] font-bold px-1.5 py-0.5 rounded-full shadow-md whitespace-nowrap">{(position * 100).toFixed(1)}%</div>
                        </div>
                        <div className="absolute inset-0 z-30">
                             <Tooltip text="Scrub to compare values" side="top">
                                <input 
                                    type="range" 
                                    min="0" 
                                    max="1" 
                                    step="0.001" 
                                    value={position} 
                                    onChange={(e) => setPosition(parseFloat(e.target.value))} 
                                    className="absolute inset-0 w-full h-full opacity-0 cursor-ew-resize"
                                />
                             </Tooltip>
                        </div>
                    </div>
                    <div className="flex justify-between text-[10px] font-bold uppercase tracking-widest text-gray-500 px-1 pt-1"><span>Linear (Reference)</span></div>
                </div>

                <div className="grid grid-cols-2 gap-4">
                     <div className="bg-gray-900 p-3 rounded border border-gray-800">
                        <div className="text-[10px] text-blue-400 mb-1 uppercase tracking-wider font-bold">Spline Value</div>
                        <div className="flex items-start gap-3">
                             <div className="relative w-24 h-24 rounded shadow-sm ring-1 ring-white/10 overflow-hidden">
                                 <div className="absolute inset-0 bg-[url('data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSI4IiBoZWlnaHQ9IjgiPjxwYXRoIGQ9Ik0wIDBoNHY0SDB6bTQgNGg0djRINHoiIGZpbGw9IiMzMzMiLz48L3N2Zz4=')] opacity-50"/>
                                 <div className="absolute inset-0" style={{ background: `oklch(${polyVal.l} ${polyVal.c} ${polyVal.h}deg / ${polyVal.a})` }} />
                             </div>
                             <div className="space-y-2 flex-1">
                                 <div className="text-xs font-mono text-gray-300">
                                    <div className="flex justify-between"><span>L</span> <span>{polyVal.l.toFixed(3)}</span></div>
                                    <div className="flex justify-between"><span>C</span> <span>{polyVal.c.toFixed(3)}</span></div>
                                    <div className="flex justify-between"><span>H</span> <span>{polyVal.h.toFixed(1)}</span></div>
                                    <div className="flex justify-between text-gray-500"><span>A</span> <span>{polyVal.a.toFixed(3)}</span></div>
                                 </div>
                                 <div className="border-t border-gray-800 pt-1 text-[10px] font-mono text-gray-500 flex justify-between items-center">
                                    <Monitor size={10} /><span>{polyRgb.r}, {polyRgb.g}, {polyRgb.b}</span>
                                 </div>
                             </div>
                        </div>
                     </div>
                     <div className="bg-gray-900 p-3 rounded border border-gray-800">
                        <div className="text-[10px] text-gray-500 mb-1 uppercase tracking-wider font-bold">Linear Value</div>
                        <div className="flex items-start gap-3">
                             <div className="relative w-24 h-24 rounded shadow-sm ring-1 ring-white/10 overflow-hidden">
                                 <div className="absolute inset-0 bg-[url('data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSI4IiBoZWlnaHQ9IjgiPjxwYXRoIGQ9Ik0wIDBoNHY0SDB6bTQgNGg0djRINHoiIGZpbGw9IiMzMzMiLz48L3N2Zz4=')] opacity-50"/>
                                 <div className="absolute inset-0" style={{ background: `oklch(${linearVal.l.toFixed(3)} ${linearVal.c.toFixed(3)} ${linearVal.h.toFixed(1)}deg / ${linearVal.a.toFixed(3)})` }} />
                             </div>
                             <div className="space-y-2 flex-1">
                                 <div className="text-xs font-mono text-gray-300">
                                    <div className="flex justify-between"><span>L</span> <span>{linearVal.l.toFixed(3)}</span></div>
                                    <div className="flex justify-between"><span>C</span> <span>{linearVal.c.toFixed(3)}</span></div>
                                    <div className="flex justify-between"><span>H</span> <span>{linearVal.h.toFixed(1)}</span></div>
                                    <div className="flex justify-between text-gray-500"><span>A</span> <span>{linearVal.a.toFixed(3)}</span></div>
                                 </div>
                                 <div className="border-t border-gray-800 pt-1 text-[10px] font-mono text-gray-500 flex justify-between items-center">
                                    <Monitor size={10} /><span>{linearRgb.r}, {linearRgb.g}, {linearRgb.b}</span>
                                 </div>
                             </div>
                        </div>
                     </div>
                </div>
            </div>

            <div className="space-y-2">
                 <label className="text-[10px] font-bold uppercase tracking-widest text-gray-500 flex justify-between">
                      <span className="flex items-center gap-2"><Layers size={12}/> Space Visualizer</span>
                      <span className="text-gray-600 flex items-center gap-1"><MousePointer2 size={10}/> Click stops to re-align view</span>
                 </label>
              <div className="w-full flex flex-col gap-3">
                {/* Spatial Editor Container */}
                <div className="w-full h-fit bg-gray-900/50 rounded-xl border border-gray-800 flex items-center justify-center p-2 gap-2 grid grid-cols-1">
                    <SplineEditor 
                      points={editorPoints} 
                      tangents={displayTangents}
                      tangentsState={tangentsState}  // ADD THIS
                      setTangents={setTangentsState} 
                      activeIndex={activePointIndex}
                      setActiveIndex={setActivePointIndex}
                      viewIndex={viewPointIndex} 
                      setViewIndex={setViewPointIndex}
                      mode={editorMode}
                      setMode={setEditorMode}
                      isRing={isRing}
                      mathMode={mathMode}
                      onReset={handleResetSpatial}
                      hasEdits={hasSpatialEdits}
                      width={240} 
                      height={180} 
                      pathPoints={editorPathPoints}
                  />

                    {/* Alpha Editor Container */}
                
                  <AlphaEditor
                    points={editorPoints}
                    tangents={displayTangents}
                    tangentsState={tangentsState}  // ADD THIS
                    setTangents={setTangentsState}
                    splineData={splineData} 
                    isRing={isRing}  // ADD THIS
                    activeIndex={activePointIndex}
                    setActiveIndex={setActivePointIndex}
                    onReset={handleResetAlpha}
                    hasEdits={hasAlphaEdits}
                    width={240} 
                    height={120}
                />
                </div>
              </div>
                 <p className="text-[10px] text-gray-600 text-center">
                    Visualizes the curve in <strong>{mathMode === 'oklch' ? 'CARTESIAN OKLAB' : 'OKLAB'}</strong> space.
                 </p>
            </div>
        </div>

        <div className="border-t border-gray-800 pt-6">
            <div className="flex justify-between items-center mb-4">
                <h2 className="text-sm font-semibold uppercase tracking-wider text-gray-400">Color Stops</h2>
                <div className="flex gap-2">
                    <button onClick={() => setIsRing(!isRing)} className={`flex items-center gap-2 text-xs font-medium px-3 py-1.5 rounded transition-colors border ${isRing ? 'bg-purple-900/50 border-purple-500 text-purple-200' : 'bg-gray-800 border-gray-700 text-gray-400 hover:text-gray-200'}`}>
                        <Repeat size={14} /> Loop
                    </button>
                    <Tooltip text="Add a new color stop" side="left">
                        <button onClick={() => { const newIdx = colors.length; setColors([...colors, [Math.random() * 360, 1, 1, 1]]); setActivePointIndex(newIdx); setViewPointIndex(newIdx); }} className="bg-blue-600 hover:bg-blue-500 text-white p-2 rounded shadow-lg shadow-blue-900/20 transition-all active:scale-95">
                            <Plus size={18} />
                        </button>
                    </Tooltip>
                </div>
            </div>
            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 xl:grid-cols-5 gap-4">
                {colors.map((color, i) => (
                    <ColorPicker 
                        key={i}
                        color={color} 
                        onChange={(c) => { const newColors = [...colors]; newColors[i] = c; setColors(newColors); }} 
                        isActive={activePointIndex === i} 
                        onActivate={() => { setActivePointIndex(i); setViewPointIndex(i); }} 
                        canRemove={colors.length > 2}
                        onRemove={() => { const newColors = colors.filter((_, idx) => idx !== i); setColors(newColors); setActivePointIndex(0); setViewPointIndex(0); }}
                        onDragStart={(e) => handleDragStart(e, i)}
                        onDragEnter={(e) => handleDragEnter(e, i)}
                        onDrop={handleDrop}
                    />
                ))}
            </div>
        </div>

        <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <CodeBlock label="Spline Gradient (CSS)" value={polyGradient} />
            <CodeBlock label="Tangent Strengths (JSON)" value={tangentStrengthList} />
            <CodeBlock label="Linear Gradient (Reference)" value={linearGradient} />
        </div>

      </div>
    </div>
  );
}