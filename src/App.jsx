import React, { useState, useMemo, useEffect, useRef } from 'react';
import { X, Plus, RotateCcw, MousePointer2, GitCommit, Copy, Check, Lock, Unlock, Monitor, FileCode, Braces, Layers, Repeat, Globe, Box, GripVertical } from 'lucide-react';

// ==========================================
// 1. Math & Color Logic
// ==========================================

const fromLinear = (val) => {
  const v = val <= 0.0031308 ? 12.92 * val : 1.055 * Math.pow(val, 1/2.4) - 0.055;
  return Math.max(0, Math.min(1, v));
};
const toLinear = (val) => val <= 0.04045 ? val / 12.92 : Math.pow((val + 0.055) / 1.055, 2.4);

// OKLCH <-> OKLAB
function oklchToOklab(l, c, h) {
  const rad = h * Math.PI / 180;
  return { 
    x: c * Math.cos(rad), // a
    y: c * Math.sin(rad), // b
    z: l                  // L
  };
}

function oklabToOklch(a, b, l) {
  const c = Math.hypot(a, b);
  const h = (Math.atan2(b, a) * 180 / Math.PI + 360) % 360;
  return { l, c, h };
}

// OKLCH -> RGB
function oklchToRgb(l, c, h) {
  const hRad = h * Math.PI / 180;
  const a = c * Math.cos(hRad);
  const b = c * Math.sin(hRad);

  const l_ = l + 0.3963377774 * a + 0.2158037573 * b;
  const m_ = l - 0.1055613458 * a - 0.0638541728 * b;
  const s_ = l - 0.0894841775 * a - 1.2914855480 * b;

  const l1 = l_ * l_ * l_;
  const m1 = m_ * m_ * m_;
  const s1 = s_ * s_ * s_;

  let lr = +4.0767416621 * l1 - 3.3077115913 * m1 + 0.2309699292 * s1;
  let lg = -1.2684380046 * l1 + 2.6097574011 * m1 - 0.3413193965 * s1;
  let lb = -0.0041960863 * l1 - 0.7034186147 * m1 + 1.7076147010 * s1;

  return { 
    r: Math.round(fromLinear(lr) * 255), 
    g: Math.round(fromLinear(lg) * 255), 
    b: Math.round(fromLinear(lb) * 255) 
  };
}

// RGB <-> HSB
function rgbToHsb(r, g, b) {
  r /= 255; g /= 255; b /= 255;
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
  return { 
    r: Math.round(f(5) * 255), 
    g: Math.round(f(3) * 255), 
    b: Math.round(f(1) * 255) 
  };
}

// RGB -> OKLCH (Approximation via Oklab)
function rgbToOklch(r, g, b) {
  const lr = toLinear(r / 255);
  const lg = toLinear(g / 255);
  const lb = toLinear(b / 255);

  const l = 0.4122214708 * lr + 0.5363325363 * lg + 0.0514459929 * lb;
  const m = 0.2119034982 * lr + 0.6806995451 * lg + 0.1073969566 * lb;
  const s = 0.0883024619 * lr + 0.2817188376 * lg + 0.6299787005 * lb;

  const l_ = Math.cbrt(l), m_ = Math.cbrt(m), s_ = Math.cbrt(s);
  const L = 0.2104542553 * l_ + 0.7936177850 * m_ - 0.0040720468 * s_;
  const a = 1.9779984951 * l_ - 2.4285922050 * m_ + 0.4505937099 * s_;
  const b_ = 0.0259040371 * l_ + 0.7827717662 * m_ - 0.8086757660 * s_;

  const C = Math.sqrt(a * a + b_ * b_);
  const H = (Math.atan2(b_, a) * 180 / Math.PI + 360) % 360;
  return { l: L, c: C, h: H };
}

// HSL Helpers
function hsbToHsl(h, s, b) {
  const l = b * (1 - s / 2);
  let sl = 0;
  if (l > 0 && l < 1) sl = (b - l) / Math.min(l, 1 - l);
  return { h, s: sl * 100, l: l * 100 };
}
function hslToHsb(h, s, l) {
  s /= 100; l /= 100;
  const b = l + s * Math.min(l, 1 - l);
  const sb = b === 0 ? 0 : 2 * (1 - l / b);
  return [h, sb, b];
}

function hsbToOklch(h, s, b) {
  const rgb = hsbToRgb(h, s, b);
  return rgbToOklch(rgb.r, rgb.g, rgb.b);
}

// --- Vector Operations (4D) ---
const vSub = (a, b) => ({ x: a.x - b.x, y: a.y - b.y, z: a.z - b.z, w: (a.w ?? 1) - (b.w ?? 1) });
const vAdd = (a, b) => ({ x: a.x + b.x, y: a.y + b.y, z: a.z + b.z, w: (a.w ?? 1) + (b.w ?? 1) });
const vScale = (v, s) => ({ x: v.x * s, y: v.y * s, z: v.z * s, w: (v.w ?? 0) * s });
const vDot = (a, b) => a.x * b.x + a.y * b.y + a.z * b.z + (a.w ?? 0) * (b.w ?? 0);
// Full 4D Length (for Spline Math)
const vLen = (v) => Math.sqrt(v.x*v.x + v.y*v.y + v.z*v.z + (v.w ?? 0)*(v.w ?? 0));
// Spatial 3D Length (for Camera Scaling - Fix #5)
const vLen3 = (v) => Math.hypot(v.x, v.y, v.z);

const vNorm = (v) => {
  const l = vLen(v);
  return l === 0 ? { x: 0, y: 0, z: 0, w: 0 } : vScale(v, 1 / l);
};
// Cross Product is 3D Spatial Only
const vCross = (a, b) => ({
  x: a.y * b.z - a.z * b.y,
  y: a.z * b.x - a.x * b.z,
  z: a.x * b.y - a.y * b.x
});

// --- Spline Logic ---

function getDelta(p1, p2, mathMode) {
    let d = vSub(p1, p2);
    if (mathMode === 'oklch') {
        let dh = d.y;
        while (dh > 180) dh -= 360;
        while (dh < -180) dh += 360;
        d.y = dh;
    }
    return d;
}

function computeTangents(points, isRing, mathMode) {
  const n = points.length;
  if (n < 2 && !isRing) return points.map(() => ({ x: 0, y: 0, z: 0 }));
  
  const tangents = new Array(n);
  
  for (let i = 0; i < n; i++) {
    let idxPrev, idxNext;
    
    if (isRing) {
        idxPrev = (i - 1 + n) % n;
        idxNext = (i + 1) % n;
    } else {
        if (i === 0) { idxPrev = 0; idxNext = 1; }
        else if (i === n - 1) { idxPrev = n - 2; idxNext = n - 1; }
        else { idxPrev = i - 1; idxNext = i + 1; }
    }
    
    const prev = points[idxPrev];
    const next = points[idxNext];
    
    const diff = getDelta(next, prev, mathMode);
    let t = vScale(diff, 0.5);

    tangents[i] = t;
  }
  return tangents;
}

function hermite(p0, m0, p1, m1, t, mathMode) {
  const t2 = t * t;
  const t3 = t2 * t;
  const h00 = 2 * t3 - 3 * t2 + 1;
  const h10 = t3 - 2 * t2 + t;
  const h01 = -2 * t3 + 3 * t2;
  const h11 = t3 - t2;
  
  let targetP1 = p1;
  if (mathMode === 'oklch') {
      let dh = p1.y - p0.y;
      while (dh > 180) dh -= 360;
      while (dh < -180) dh += 360;
      targetP1 = { ...p1, y: p0.y + dh };
  }

  // Alpha (w) interpolation
  const w0 = p0.w ?? 1;
  const w1 = targetP1.w ?? 1;
  const mw0 = m0.w ?? 0;
  const mw1 = m1.w ?? 0;

  const val = {
    x: h00 * p0.x + h10 * m0.x + h01 * targetP1.x + h11 * m1.x,
    y: h00 * p0.y + h10 * m0.y + h01 * targetP1.y + h11 * m1.y,
    z: h00 * p0.z + h10 * m0.z + h01 * targetP1.z + h11 * m1.z,
    w: h00 * w0 + h10 * mw0 + h01 * w1 + h11 * mw1
  };

  if (mathMode === 'oklch') {
      val.y = (val.y % 360 + 360) % 360;
  }
  return val;
}

function toCartesian(pt, mathMode) {
    if (mathMode === 'oklch') {
        return oklchToOklab(pt.z, pt.x, pt.y);
    }
    return pt;
}

function unwrapPoint(source, target, mathMode) {
    if (mathMode !== 'oklch') return target;
    const d = target.y - source.y;
    let shift = 0;
    if (d > 180) shift = -360;
    else if (d < -180) shift = 360;
    if (shift !== 0) return { ...target, y: target.y + shift };
    return target;
}

function getLocalBasis(points, index, isRing, mathMode) {
  if (points.length < 2) return { center: toCartesian(points[0], mathMode), u: {x:1,y:0,z:0}, v: {x:0,y:1,z:0} };
  
  let basisIndex = index;
  if (!isRing) {
      if (index === 0) basisIndex = 1;
      else if (index === points.length - 1) basisIndex = points.length - 2;
  }
  basisIndex = Math.max(0, Math.min(points.length - 1, basisIndex));

  const prevIdx = isRing ? (basisIndex - 1 + points.length) % points.length : Math.max(0, basisIndex - 1);
  const nextIdx = isRing ? (basisIndex + 1) % points.length : Math.min(points.length - 1, basisIndex + 1);
  
  const currCart = toCartesian(points[basisIndex], mathMode);
  const prevCart = toCartesian(points[prevIdx], mathMode);
  const nextCart = toCartesian(points[nextIdx], mathMode);

  let tangent = vSub(nextCart, prevCart); 
  let u = vNorm(tangent);
  if (vLen3(u) < 0.0001) u = { x: 1, y: 0, z: 0 };

  const v1 = vSub(currCart, prevCart);
  const v2 = vSub(nextCart, currCart);
  let normal = vCross(v1, v2);

  if (vLen3(normal) < 0.0001) {
     const obliqueUp = vNorm({ x: 1, y: 0, z: 1 });
     normal = vCross(u, obliqueUp); 
     if (vLen3(normal) < 0.0001) normal = vCross(u, { x: 1, y: 0, z: 0 }); 
  }
  normal = vNorm(normal);
  const v = vNorm(vCross(normal, u));
  const realCenter = toCartesian(points[index], mathMode);
  return { center: realCenter, u: u, v: v, n: normal }; 
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

function NumberInput({ value, onChange, min, max, label }) {
    return (
        <div className="flex-1 min-w-0 group relative" onDragStart={(e) => { e.preventDefault(); e.stopPropagation(); }} draggable>
            <input 
                type="text" 
                value={value}
                onChange={(e) => {
                    let val = parseFloat(e.target.value);
                    if (isNaN(val)) val = 0; 
                    onChange(val);
                }}
                onBlur={(e) => {
                    let val = parseFloat(e.target.value);
                    if (isNaN(val)) val = min;
                    val = Math.max(min, Math.min(max, val));
                    onChange(val);
                }}
                className="w-full bg-gray-950 border border-gray-800 rounded px-1.5 py-0.5 text-[10px] font-mono text-gray-300 focus:border-blue-500 focus:outline-none text-right"
            />
            {label && <span className="absolute left-1 top-1/2 -translate-y-1/2 text-[8px] text-gray-600 pointer-events-none font-bold select-none">{label}</span>}
        </div>
    );
}

function ColorInputs({ color, onChange }) {
    const [h, s, b, a = 1] = color;
    const rgb = hsbToRgb(h, s, b);
    const hsl = hsbToHsl(h, s, b);
    const oklch = rgbToOklch(rgb.r, rgb.g, rgb.b);
    const [copied, setCopied] = useState(null);

    const updateFromRgb = (k, v) => onChange([...rgbToHsb({ ...rgb, [k]: v }.r, { ...rgb, [k]: v }.g, { ...rgb, [k]: v }.b), a]);
    const updateFromHsl = (k, v) => onChange([...hslToHsb({ ...hsl, [k]: v }.h, { ...hsl, [k]: v }.s, { ...hsl, [k]: v }.l), a]);
    const updateFromOklch = (k, v) => {
        const n = { ...oklch, [k]: v };
        const lab = oklchToOklab(n.l, n.c, n.h);
        const rgbVal = oklchToRgb(lab.z, oklabToOklch(lab.x, lab.y, lab.z).c, oklabToOklch(lab.x, lab.y, lab.z).h);
        onChange([...rgbToHsb(rgbVal.r, rgbVal.g, rgbVal.b), a]);
    };
    const updateAlpha = (v) => onChange([h, s, b, v]);

    const copy = (txt, key) => {
        navigator.clipboard.writeText(txt);
        setCopied(key);
        setTimeout(() => setCopied(null), 1500);
    };

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
                <span className="text-[9px] font-bold text-gray-600">HSL</span>
                <NumberInput value={Math.round(hsl.h)} min={0} max={360} onChange={(v) => updateFromHsl('h', v)} />
                <NumberInput value={Math.round(hsl.s)} min={0} max={100} onChange={(v) => updateFromHsl('s', v)} />
                <NumberInput value={Math.round(hsl.l)} min={0} max={100} onChange={(v) => updateFromHsl('l', v)} />
                <button onClick={() => copy(`hsla(${hsl.h.toFixed(0)}, ${hsl.s.toFixed(0)}%, ${hsl.l.toFixed(0)}%, ${a.toFixed(2)})`, 'hsl')} className="text-gray-500 hover:text-white">{copied === 'hsl' ? <Check size={10} className="text-green-500"/> : <Copy size={10}/>}</button>
            </div>
             <div className="grid grid-cols-[24px_1fr_1fr_1fr_14px] gap-1 items-center">
                <span className="text-[9px] font-bold text-gray-600">LCH</span>
                <NumberInput value={oklch.l.toFixed(3)} step={0.01} min={0} max={1} onChange={(v) => updateFromOklch('l', v)} />
                <NumberInput value={oklch.c.toFixed(3)} step={0.01} min={0} max={0.4} onChange={(v) => updateFromOklch('c', v)} />
                <NumberInput value={Math.round(oklch.h)} min={0} max={360} onChange={(v) => updateFromOklch('h', v)} />
                <button onClick={() => copy(`oklch(${oklch.l.toFixed(3)} ${oklch.c.toFixed(3)} ${oklch.h.toFixed(1)}deg / ${a.toFixed(2)})`, 'lch')} className="text-gray-500 hover:text-white">{copied === 'lch' ? <Check size={10} className="text-green-500"/> : <Copy size={10}/>}</button>
            </div>
            <div className="grid grid-cols-[24px_1fr_14px] gap-1 items-center mt-0.5 border-t border-gray-800 pt-1">
                <span className="text-[9px] font-bold text-gray-500">ALP</span>
                <NumberInput value={a.toFixed(2)} step={0.01} min={0} max={1} onChange={updateAlpha} />
                <div />
            </div>
        </div>
    );
}

function SplineEditor({ 
  points, 
  tangents, 
  setTangents, 
  activeIndex, 
  setActiveIndex, 
  viewIndex, 
  setViewIndex, 
  mode, 
  setMode,
  isRing,
  mathMode,
  width = 320, 
  height = 320 
}) {
  const svgRef = useRef(null);
  const [dragState, setDragState] = useState(null); 

  const basis = useMemo(() => {
    return getLocalBasis(points, viewIndex ?? 0, isRing, mathMode);
  }, [points, viewIndex, isRing, mathMode]);

  const scale = useMemo(() => {
    let maxDist = 0.05; 
    let basisIndex = viewIndex ?? 0;
    if (!isRing) {
        if (basisIndex === 0) basisIndex = 1;
        else if (basisIndex === points.length - 1) basisIndex = points.length - 2;
    }
    basisIndex = Math.max(0, Math.min(points.length - 1, basisIndex));
    const basisCenter = toCartesian(points[basisIndex], mathMode);

    const prevIdx = isRing ? (basisIndex - 1 + points.length) % points.length : Math.max(0, basisIndex - 1);
    const nextIdx = isRing ? (basisIndex + 1) % points.length : Math.min(points.length - 1, basisIndex + 1);
    
    const prevCart = toCartesian(points[prevIdx], mathMode);
    const nextCart = toCartesian(points[nextIdx], mathMode);

    // Fix #5: Use 3D length (spatial) for zoom scaling to ignore alpha changes
    const d1 = vLen3(vSub(prevCart, basisCenter));
    const d2 = vLen3(vSub(nextCart, basisCenter));
    
    maxDist = Math.max(d1, d2);
    if(maxDist === 0) maxDist = 0.1;

    const viewportRadius = Math.min(width, height) / 2;
    return (viewportRadius * 0.6) / maxDist;
  }, [points, width, height, isRing, mathMode, viewIndex]);

  const offsetX = width / 2;
  const offsetY = height / 2;

  const toScreen = (pt) => {
    const uPt = unwrapPoint(basis.center, pt, mathMode);
    const rel = vSub(toCartesian(uPt, mathMode), basis.center);
    const x = vDot(rel, basis.u);
    const y = vDot(rel, basis.v);
    return { x: x * scale + offsetX, y: offsetY - y * scale };
  };

  const fromScreen = (sx, sy) => {
    const x = (sx - offsetX) / scale;
    const y = (offsetY - sy) / scale;
    const wu = vScale(basis.u, x);
    const wv = vScale(basis.v, y);
    return vAdd(basis.center, vAdd(wu, wv));
  };

  const handlePointerDown = (e, index) => {
    e.stopPropagation();
    e.target.setPointerCapture(e.pointerId);
    setActiveIndex(index);
    if(setViewIndex) setViewIndex(index);
    
    const pt = points[index];
    let handlePos3D;
    if (mathMode === 'oklch') {
       const polarHandle = vAdd(pt, vScale(tangents[index], 1/3));
       handlePos3D = polarHandle; 
    } else {
       handlePos3D = vAdd(pt, vScale(tangents[index], 1/3));
    }

    const handlePosScreen = toScreen(handlePos3D);
    const svgRect = svgRef.current.getBoundingClientRect();
    const ex = e.clientX - svgRect.left;
    const ey = e.clientY - svgRect.top;
    
    setDragState({ 
        index, 
        pointerId: e.pointerId,
        initialTangent: tangents[index],
        dragOffset: { x: ex - handlePosScreen.x, y: ey - handlePosScreen.y }
    });
  };

  const handlePointerMove = (e) => {
    if (!dragState) return;
    const svgRect = svgRef.current.getBoundingClientRect();
    const ex = e.clientX - svgRect.left;
    const ey = e.clientY - svgRect.top;
    const correctedX = ex - dragState.dragOffset.x;
    const correctedY = ey - dragState.dragOffset.y;
    
    const newHandlePosCart = fromScreen(correctedX, correctedY);
    const pt = points[dragState.index];
    
    let newTangent;

    if (mathMode === 'oklch') {
        const polarHandle = oklabToOklch(newHandlePosCart.x, newHandlePosCart.y, newHandlePosCart.z);
        const pHandle = { x: polarHandle.c, y: polarHandle.h, z: polarHandle.l };
        let dh = pHandle.y - pt.y;
        while (dh > 180) dh -= 360; while (dh < -180) dh += 360;
        const unscaledT = { x: pHandle.x - pt.x, y: dh, z: pHandle.z - pt.z };
        newTangent = vScale(unscaledT, 3);
        
        if (mode === 'strength' && dragState.initialTangent) {
            const initT = dragState.initialTangent;
            const len = vLen(initT);
            if (len > 0.0001) {
                const dir = vScale(initT, 1/len);
                const currentDelta = vScale(newTangent, 1/3);
                const projLen = vDot(currentDelta, dir);
                newTangent = vScale(dir, projLen * 3);
            }
        }

    } else {
        const ptCart = toCartesian(pt, mathMode);
        let delta = vSub(newHandlePosCart, ptCart);
        if (mode === 'strength' && dragState.initialTangent) {
            const initT = dragState.initialTangent;
            const len = vLen(initT);
            if (len > 0.0001) {
                const dir = vScale(initT, 1/len);
                const projectedLen = vDot(delta, dir);
                delta = vScale(dir, projectedLen);
            }
        }
        newTangent = vScale(delta, 3);
    }

    const newTangents = [...tangents];
    newTangents[dragState.index] = newTangent;
    setTangents(newTangents);
  };

  const handlePointerUp = (e) => {
    if(dragState) { e.target.releasePointerCapture(e.pointerId); setDragState(null); }
  };

  const pathData = useMemo(() => {
    if (points.length < 2) return "";
    let d = "";
    const samples = 40; 
    const segments = isRing ? points.length : points.length - 1;
    for (let i = 0; i < segments; i++) {
      const p0 = points[i];
      const p1 = points[(i+1) % points.length];
      const m0 = tangents[i];
      const m1 = tangents[(i+1) % points.length];
      for (let s = 0; s <= samples; s++) {
        const t = s / samples;
        const ptNative = hermite(p0, m0, p1, m1, t, mathMode);
        const screenPt = toScreen(ptNative);
        if (i === 0 && s === 0) d += `M ${screenPt.x} ${screenPt.y}`;
        else d += ` L ${screenPt.x} ${screenPt.y}`;
      }
    }
    return d;
  }, [points, tangents, basis, scale, isRing, mathMode]);

  return (
    <div className="bg-gray-900 rounded-lg shadow-inner inline-block select-none relative overflow-hidden border border-gray-800 group/editor">
       <div className="absolute top-3 left-3 text-[10px] text-gray-500 pointer-events-none z-10 font-mono">
         <div className="font-bold text-gray-300">{mathMode.toUpperCase()} SPACE EDITOR</div>
         <div>Aligning to Stop #{ (viewIndex ?? 0) + 1 }</div>
         <div className="text-[9px] opacity-50 mt-1">
             {mathMode === 'oklch' ? 'X: C, Y: H, Z: L' : 'X: a, Y: b, Z: L'}
         </div>
       </div>
       <div className="absolute top-3 right-3 z-20 flex bg-gray-950 rounded border border-gray-800 p-0.5">
          <Tooltip text={`Free Mode: Adjust angle & length (${mathMode})`} side="bottom">
              <button onClick={() => setMode('free')} className={`p-1.5 rounded ${mode === 'free' ? 'bg-blue-600 text-white' : 'text-gray-400 hover:text-gray-200'}`}><Unlock size={14} /></button>
          </Tooltip>
          <Tooltip text="Strength Mode: Adjust tension only" side="bottom">
              <button onClick={() => setMode('strength')} className={`p-1.5 rounded ${mode === 'strength' ? 'bg-blue-600 text-white' : 'text-gray-400 hover:text-gray-200'}`}><Lock size={14} /></button>
          </Tooltip>
       </div>
       <svg ref={svgRef} width={width} height={height} className="block cursor-default" onPointerMove={handlePointerMove} onPointerUp={handlePointerUp}>
         <defs>
            <pattern id="grid" width="20" height="20" patternUnits="userSpaceOnUse"><path d="M 20 0 L 0 0 0 20" fill="none" stroke="#1f2937" strokeWidth="1"/></pattern>
         </defs>
         <rect width="100%" height="100%" fill="url(#grid)" />
         <line x1={offsetX - 10} y1={offsetY} x2={offsetX + 10} y2={offsetY} stroke="#374151" />
         <line x1={offsetX} y1={offsetY - 10} x2={offsetX} y2={offsetY + 10} stroke="#374151" />
         <path d={pathData} stroke="white" strokeWidth="2" fill="none" strokeLinecap="round" />
         {points.map((pt, i) => {
            if (!tangents[i]) return null;
            const isActive = i === activeIndex;
            const screenPt = toScreen(pt);
            let handlePos3D;
            if(mathMode === 'oklch') {
               handlePos3D = vAdd(pt, vScale(tangents[i], 1/3));
            } else {
               handlePos3D = vAdd(pt, vScale(tangents[i], 1/3));
            }
            const screenHandle = toScreen(handlePos3D);
            const deltaScreen = vSub(screenHandle, screenPt);
            const screenBack = { x: screenPt.x - deltaScreen.x, y: screenPt.y - deltaScreen.y };
            const isFocus = i === viewIndex;
            const opacity = isFocus ? 1 : (isActive ? 0.8 : 0.25);
            return (
              <g key={i} opacity={opacity} className="transition-opacity duration-200">
                <line x1={screenBack.x} y1={screenBack.y} x2={screenHandle.x} y2={screenHandle.y} stroke={isActive ? "#9ca3af" : "#4b5563"} strokeWidth="1" pointerEvents="none" />
                <circle cx={screenPt.x} cy={screenPt.y} r={4} fill={isActive ? "#fff" : "#9ca3af"} 
                    onPointerDown={(e) => { e.stopPropagation(); setActiveIndex(i); if(setViewIndex) setViewIndex(i); }}
                    className="cursor-pointer hover:fill-white" />
                {isActive && mode === 'strength' && (
                     <line x1={screenPt.x} y1={screenPt.y} x2={screenPt.x + (screenHandle.x - screenPt.x) * 10} y2={screenPt.y + (screenHandle.y - screenPt.y) * 10} stroke="rgba(59, 130, 246, 0.2)" strokeDasharray="2 2" pointerEvents="none" />
                )}
                <circle cx={screenHandle.x} cy={screenHandle.y} r={6} fill={mode === 'strength' ? "#ec4899" : "#3b82f6"} stroke="white" strokeWidth={isActive ? 2 : 1}
                  className={`${isActive ? 'cursor-move' : 'cursor-pointer'} hover:scale-125`}
                  onPointerDown={(e) => handlePointerDown(e, i)} />
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
        
        {/* Fix 1: Taller square (h-36) and Stop drag prop */}
        <div 
            className="flex gap-2 h-36" 
            onDragStart={(e) => { e.preventDefault(); e.stopPropagation(); }} 
            draggable
        >
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
            {/* Fix 3: Flex header wrapping issue */}
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
  const [tangentsState, setTangentsState] = useState([]);
  const [isManual, setIsManual] = useState(false);
  const [position, setPosition] = useState(0.5); 
  const [steps, setSteps] = useState(40); 
  const [editorMode, setEditorMode] = useState('strength');
  const [copied, setCopied] = useState(false);
  const [isRing, setIsRing] = useState(false);
  const [mathMode, setMathMode] = useState('oklab'); 

  // Drag & Drop State
  const dragItem = useRef(null);
  const dragOverItem = useRef(null);

  const handleDragStart = (e, index) => {
      dragItem.current = index;
  };

  const handleDragEnter = (e, index) => {
      dragOverItem.current = index;
  };

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
  };

  const minSteps = Math.max(2, colors.length + (isRing ? 1 : 0));
  
  useEffect(() => {
     if(steps < minSteps) setSteps(minSteps);
  }, [minSteps, steps]);

  const points = useMemo(() => {
    return colors.map(([h, s, b, a]) => {
      const { l, c, h: oh } = hsbToOklch(h, s, b);
      // Pass 'a' as 4th component
      const alpha = a ?? 1; 
      if (mathMode === 'oklch') {
          return { x: c, y: oh, z: l, w: alpha }; 
      } else {
          const lab = oklchToOklab(l, c, oh);
          return { ...lab, w: alpha }; 
      }
    });
  }, [colors, mathMode]);

  useEffect(() => {
    setTangentsState([]);
    setIsManual(false);
  }, [points.length, isRing, mathMode]);

  const effectiveTangents = useMemo(() => {
    if (tangentsState.length === points.length) return tangentsState;
    return computeTangents(points, isRing, mathMode);
  }, [points, tangentsState, isRing, mathMode]);

  const handleTangentUpdate = (newT) => {
    setTangentsState(newT);
    setIsManual(true);
  };

  const tangentStrengthList = useMemo(() => {
      const defaults = computeTangents(points, isRing, mathMode);
      const strengths = effectiveTangents.map((t, i) => {
          const defLen = vLen(defaults[i]);
          const curLen = vLen(t);
          if (defLen < 0.00001) return 1.0;
          return parseFloat((curLen / defLen).toFixed(3));
      });
      return JSON.stringify(strengths, null, 0).replace(/,/g, ', ');
  }, [points, effectiveTangents, isRing, mathMode]);


  const { polyGradient, linearGradient, polyVal, linearVal, polyRgb, linearRgb } = useMemo(() => {
    if (points.length < 2 && !isRing) return { polyGradient:'', linearGradient: '', polyVal: null, linearVal: null };
    
    // --- Advanced Spline Math: Non-Uniform Catmull-Rom (Chordal) ---
    const n = points.length;
    // For a ring, we have n segments (0->1, 1->2, ... n-1->0).
    // For a line, we have n-1 segments.
    const segments = isRing ? n : n - 1;
    
    // 1. Calculate Chords (Distances) in NATIVE space
    // This fixes the bug: we calculate distance based on the dimensions we are interpolating.
    const chords = [];
    for (let i = 0; i < segments; i++) {
        const pA = points[i];
        const pB = points[(i + 1) % n];
        // Use getDelta to handle Hue wrapping automatically in OKLCH
        let dist = vLen(getDelta(pB, pA, mathMode));
        if (dist < 0.0001) dist = 0.0001; // Avoid divide by zero
        chords.push(dist);
    }

    // 2. Compute Tangents (m)
    // Non-Uniform Catmull-Rom Tangents adapted for Hermite Basis
    const realTangents = new Array(n);
    
    for (let i = 0; i < n; i++) {
        // Handle Open Ends
        if (!isRing) {
            if (i === 0) {
                // Start: Linear projection to 2nd point
                realTangents[i] = getDelta(points[1], points[0], mathMode);
                continue;
            } else if (i === n - 1) {
                // End: Linear projection from 2nd to last point
                realTangents[i] = getDelta(points[n-1], points[n-2], mathMode);
                continue;
            }
        }

        // Internal / Ring Points
        const prevIdx = (i - 1 + n) % n;
        const nextIdx = (i + 1) % n;
        
        const pPrev = points[prevIdx];
        const pCurr = points[i];
        const pNext = points[nextIdx];

        // Indices for chords array
        // Chord[k] is distance between P[k] and P[k+1]
        // Incoming segment length: distance(prev, curr)
        const dt0 = chords[(prevIdx) % segments]; 
        // Outgoing segment length: distance(curr, next)
        const dt1 = chords[i % segments];         

        // Vectors between points (handling hue wrap)
        const v1 = getDelta(pCurr, pPrev, mathMode); // P_i - P_i-1
        const v2 = getDelta(pNext, pCurr, mathMode); // P_i+1 - P_i
        
        // Weighted Tangent Formula (Centripetal/Chordal)
        // Scaled to match the Hermite interval [0,1] for the *outgoing* segment
        // Multiplier derived from: T = (dt1^2 * v1 + dt0^2 * v2) / (dt0 * dt1 * (dt0 + dt1))
        // To normalize for Hermite (velocity * duration), we multiply by dt1.
        
        const tangent = {
            x: (dt1*dt1 * v1.x + dt0*dt0 * v2.x) / (dt0 * (dt0 + dt1)),
            y: (dt1*dt1 * v1.y + dt0*dt0 * v2.y) / (dt0 * (dt0 + dt1)),
            z: (dt1*dt1 * v1.z + dt0*dt0 * v2.z) / (dt0 * (dt0 + dt1)),
            w: (dt1*dt1 * v1.w + dt0*dt0 * v2.w) / (dt0 * (dt0 + dt1))
        };
        
        realTangents[i] = tangent;
    }

    // --- 3. Generate Spline Gradient Stops ---
    const gradientStops = [];
    
    // Gradient Generation Loop
    // We iterate from 0 to Steps.
    for (let i = 0; i <= steps; i++) {
        const globalT = i / steps;
        const scaled = globalT * segments;
        
        // Determine which segment we are in [0 ... segments-1]
        let segIdx = Math.floor(scaled);
        if (segIdx >= segments) segIdx = segments - 1; 
        
        const localT = scaled - segIdx; // 0..1 within segment
        
        const pIdx0 = segIdx;
        const pIdx1 = (segIdx + 1) % n;
        
        const p0 = points[pIdx0];
        const p1 = points[pIdx1];
        
        let t0, t1;
        
        if (isManual) {
            // Manual overrides from the editor handles
            t0 = effectiveTangents[pIdx0];
            t1 = effectiveTangents[pIdx1];
        } else {
            // Calculated Tangents
            // Hermite requires:
            // m0 = Tangent at Start * Duration
            // m1 = Tangent at End * Duration
            // Our realTangents[i] is already scaled for the segment *starting* at i.
            
            t0 = realTangents[pIdx0];
            
            // For the incoming tangent at P1, we have the tangent calculated for the NEXT segment.
            // We must rescale it for the CURRENT segment length.
            const t1Raw = realTangents[pIdx1];
            const distCurr = chords[pIdx0]; // Length of this segment
            const distNext = chords[pIdx1 % segments]; // Length of next segment
            
            // t1 = t1Raw * (distCurr / distNext)
            const scaleFactor = (distNext > 0.0001) ? (distCurr / distNext) : 1;
            t1 = vScale(t1Raw, scaleFactor);
        }

        const pt = hermite(p0, t0, p1, t1, localT, mathMode);
        
        let l, c, h;
        if (mathMode === 'oklch') {
            l = pt.z; c = pt.x; h = pt.y;
        } else {
            ({ l, c, h } = oklabToOklch(pt.x, pt.y, pt.z));
        }
        const alpha = pt.w ?? 1;
        
        gradientStops.push(`oklch(${l.toFixed(4)} ${c.toFixed(4)} ${h.toFixed(2)}deg / ${alpha.toFixed(3)}) ${(globalT * 100).toFixed(1)}%`);
    }
    const polyGradient = `linear-gradient(to right, ${gradientStops.join(', ')})`;

    // --- 4. Linear Reference ---
    const linStops = [];
    const listToMap = isRing ? [...colors, colors[0]] : colors;
    
    listToMap.forEach((c, i, arr) => {
        const {l, c:cn, h} = hsbToOklch(c[0], c[1], c[2]);
        const alpha = c[3] ?? 1;
        linStops.push(`oklch(${l.toFixed(4)} ${cn.toFixed(4)} ${h.toFixed(2)}deg / ${alpha.toFixed(3)}) ${(i / (arr.length - 1) * 100).toFixed(1)}%`);
    });
    const linearGradient = `linear-gradient(to right, ${linStops.join(', ')})`;

    // --- 5. Values at Position ---
    const scaledP = position * segments;
    let segP = Math.floor(scaledP); if (segP >= segments) segP = segments - 1;
    const uP = scaledP - segP;
    
    const idx0 = segP;
    const idx1 = (segP + 1) % n;

    // Recalculate T for the single point (could refactor to shared func, but inline is safe)
    let valT0, valT1;
    if (isManual) {
        valT0 = effectiveTangents[idx0];
        valT1 = effectiveTangents[idx1];
    } else {
        valT0 = realTangents[idx0];
        const distCurr = chords[idx0];
        const distNext = chords[idx1 % segments];
        const scaleFactor = (distNext > 0.0001) ? (distCurr / distNext) : 1;
        valT1 = vScale(realTangents[idx1], scaleFactor);
    }

    const polyPt = hermite(points[idx0], valT0, points[idx1], valT1, uP, mathMode);
    
    let polyL, polyC, polyH;
    if (mathMode === 'oklch') {
        polyL = polyPt.z; polyC = polyPt.x; polyH = polyPt.y;
    } else {
        ({ l: polyL, c: polyC, h: polyH } = oklabToOklch(polyPt.x, polyPt.y, polyPt.z));
    }
    const polyVal = { l: polyL, c: polyC, h: polyH, a: polyPt.w ?? 1 };
    const polyRgb = oklchToRgb(polyL, polyC, polyH);

    // Linear Value
    const pA = points[idx0]; 
    const pB = points[idx1];
    let linPt;
    
    if (mathMode === 'oklch') {
         let dh = pB.y - pA.y;
         while(dh > 180) dh -= 360; while(dh < -180) dh += 360;
         linPt = { 
             x: pA.x + (pB.x - pA.x) * uP, 
             y: pA.y + dh * uP, 
             z: pA.z + (pB.z - pA.z) * uP,
             w: (pA.w ?? 1) + ((pB.w ?? 1) - (pA.w ?? 1)) * uP
         };
         linPt.y = (linPt.y % 360 + 360) % 360;
    } else {
         linPt = { 
             x: pA.x + (pB.x - pA.x) * uP, 
             y: pA.y + (pB.y - pA.y) * uP, 
             z: pA.z + (pB.z - pA.z) * uP,
             w: (pA.w ?? 1) + ((pB.w ?? 1) - (pA.w ?? 1)) * uP
         };
    }
    
    let linL, linC, linH;
    if (mathMode === 'oklch') {
        linL = linPt.z; linC = linPt.x; linH = linPt.y;
    } else {
        ({ l: linL, c: linC, h: linH } = oklabToOklch(linPt.x, linPt.y, linPt.z));
    }
    const linearVal = { l: linL, c: linC, h: linH, a: linPt.w };
    const linearRgb = oklchToRgb(linL, linC, linH);

    return { polyGradient, linearGradient, polyVal, linearVal, polyRgb, linearRgb };
  }, [points, effectiveTangents, colors, position, steps, isRing, mathMode, isManual]);

  const copyToClipboard = () => {
    navigator.clipboard.writeText(polyGradient);
    setCopied(true);
    setTimeout(() => setCopied(false), 2000);
  };

  return (
    <div className="min-h-screen bg-gray-950 text-gray-200 py-8 px-4 font-sans flex justify-center">
      <div className="w-full max-w-6xl flex flex-col gap-6">
        
        {/* Header */}
        <header className="flex flex-col md:flex-row justify-between items-start md:items-end border-b border-gray-800 pb-4 gap-4">
          <div>
            <h1 className="text-2xl font-bold text-white mb-1 tracking-tight flex items-center gap-2">
               <GitCommit className="text-blue-500" />
               OKLCH Spline Gradient
            </h1>
            <p className="text-sm text-gray-400">Generates smooth Catmull-Rom gradients.</p>
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

            {isManual && (
                <Tooltip text="Reset to default smooth curves" side="bottom">
                    <button onClick={() => { setTangentsState([]); setIsManual(false); }} className="flex items-center gap-2 text-xs font-medium bg-gray-800 hover:bg-gray-700 px-3 py-1.5 rounded text-blue-300 transition-colors whitespace-nowrap">
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

        {/* Top Section */}
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-8">
            <div className="space-y-6">
                <div className="space-y-1">
                    <div className="flex justify-between items-center text-[10px] font-bold uppercase tracking-widest text-gray-500 px-1">
                        <span className="flex items-center gap-2">
                             Spline ({mathMode.toUpperCase()})
                             <Tooltip text="Gradient Resolution (Steps)" side="right">
                                <input type="range" min={minSteps} max="100" value={steps} onChange={(e) => setSteps(Number(e.target.value))} className="w-16 accent-blue-500 cursor-pointer"/>
                             </Tooltip>
                             <span className="text-gray-600 font-mono normal-case">{steps} stops</span>
                        </span>
                    </div>
                    
                    {/* Visual Gradient Area */}
                    <div className="relative group select-none rounded-lg shadow-lg ring-1 ring-white/10 overflow-hidden">
                        {/* Checkerboard for Alpha */}
                        <div className="absolute inset-0 bg-[url('data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxNiIgaGVpZ2h0PSIxNiI+PHBhdGggZD0iTTAgMGg4djhIMHptOCA4aDh2OEg4eiIgZmlsbD0iIzIyMiIvPjwvc3ZnPg==')] z-0" />
                        
                        <div className="h-16 w-full z-10 relative" style={{ background: polyGradient }} />
                        <div className="h-8 w-full z-10 relative border-t border-black/20" style={{ background: linearGradient }} />
                        
                        {/* Indicator Line */}
                        <div className="absolute top-0 bottom-0 z-20 pointer-events-none flex flex-col items-center justify-center w-full" style={{ left: `${position * 100}%`, transform: 'translateX(-50%)' }}>
                            <div className="w-0.5 h-full bg-white shadow-[0_0_10px_rgba(0,0,0,0.5)]" />
                            <div className="absolute top-1/2 -translate-y-1/2 bg-white text-gray-950 text-[10px] font-bold px-1.5 py-0.5 rounded-full shadow-md whitespace-nowrap">{(position * 100).toFixed(1)}%</div>
                        </div>

                        {/* Slider Input */}
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

                {polyVal && (
                <div className="grid grid-cols-2 gap-4">
                     {/* Spline Stats */}
                     <div className="bg-gray-900 p-3 rounded border border-gray-800">
                        <div className="text-[10px] text-blue-400 mb-1 uppercase tracking-wider font-bold">Spline Value</div>
                        <div className="flex items-start gap-3">
                             <div className="relative w-10 h-10 rounded shadow-sm ring-1 ring-white/10 overflow-hidden">
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
                     {/* Linear Stats */}
                     <div className="bg-gray-900 p-3 rounded border border-gray-800">
                        <div className="text-[10px] text-gray-500 mb-1 uppercase tracking-wider font-bold">Linear Value</div>
                        <div className="flex items-start gap-3">
                             <div className="relative w-10 h-10 rounded shadow-sm ring-1 ring-white/10 overflow-hidden">
                                 <div className="absolute inset-0 bg-[url('data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSI4IiBoZWlnaHQ9IjgiPjxwYXRoIGQ9Ik0wIDBoNHY0SDB6bTQgNGg0djRINHoiIGZpbGw9IiMzMzMiLz48L3N2Zz4=')] opacity-50"/>
                                 <div className="absolute inset-0" style={{ background: `oklch(${linearVal.l} ${linearVal.c} ${linearVal.h}deg / ${linearVal.a})` }} />
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
                )}
            </div>

            {/* 3D Editor */}
            <div className="space-y-2">
                 <label className="text-[10px] font-bold uppercase tracking-widest text-gray-500 flex justify-between">
                      <span className="flex items-center gap-2"><Layers size={12}/> Space Visualizer</span>
                      <span className="text-gray-600 flex items-center gap-1"><MousePointer2 size={10}/> Click stops to re-align view</span>
                 </label>
                 <div className="w-full h-64 bg-gray-900/50 rounded-xl border border-gray-800 flex items-center justify-center">
                      <SplineEditor 
                          points={points} 
                          tangents={effectiveTangents} 
                          setTangents={handleTangentUpdate}
                          activeIndex={activePointIndex}
                          setActiveIndex={setActivePointIndex}
                          viewIndex={viewPointIndex} 
                          setViewIndex={setViewPointIndex}
                          mode={editorMode}
                          setMode={setEditorMode}
                          isRing={isRing}
                          mathMode={mathMode}
                          width={320} 
                          height={250} 
                      />
                 </div>
                 <p className="text-[10px] text-gray-600 text-center">
                    Visualizes the curve in <strong>{mathMode === 'oklch' ? 'TRUE CARTESIAN (Projected from Polar)' : 'OKLAB'}</strong> space.
                 </p>
            </div>
        </div>

        {/* Color Stops */}
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

        {/* Output Area */}
        <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <CodeBlock label="Spline Gradient (CSS)" value={polyGradient} />
            <CodeBlock label="Tangent Strengths (JSON)" value={tangentStrengthList} />
            <CodeBlock label="Linear Gradient (Reference)" value={linearGradient} />
        </div>

      </div>
    </div>
  );
}