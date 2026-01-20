import React, { useState, useMemo, useEffect, useRef } from 'react';
import { X, Plus, RotateCcw, MousePointer2, GitCommit, Copy, Check, Lock, Unlock, Monitor, FileCode, Info, Braces } from 'lucide-react';

// ==========================================
// 1. Math & Color Logic
// ==========================================

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

  const fromLinear = (val) => {
    const v = val <= 0.0031308 ? 12.92 * val : 1.055 * Math.pow(val, 1/2.4) - 0.055;
    return Math.max(0, Math.min(1, v));
  };

  return { 
    r: Math.round(fromLinear(lr) * 255), 
    g: Math.round(fromLinear(lg) * 255), 
    b: Math.round(fromLinear(lb) * 255) 
  };
}

function hsbToOklch(h, s, b) {
  const c = b * s;
  const x = c * (1 - Math.abs(((h / 60) % 2) - 1));
  const m = b - c;

  let r, g, bl;
  if (h < 60) [r, g, bl] = [c, x, 0];
  else if (h < 120) [r, g, bl] = [x, c, 0];
  else if (h < 180) [r, g, bl] = [0, c, x];
  else if (h < 240) [r, g, bl] = [0, x, c];
  else if (h < 300) [r, g, bl] = [x, 0, c];
  else [r, g, bl] = [c, 0, x];

  r += m; g += m; bl += m;

  const toLinear = (val) => val <= 0.04045 ? val / 12.92 : Math.pow((val + 0.055) / 1.055, 2.4);
  const lr = toLinear(r), lg = toLinear(g), lb = toLinear(bl);

  const l = 0.4122214708 * lr + 0.5363325363 * lg + 0.0514459929 * lb;
  const m1 = 0.2119034982 * lr + 0.6806995451 * lg + 0.1073969566 * lb;
  const s1 = 0.0883024619 * lr + 0.2817188376 * lg + 0.6299787005 * lb;

  const l_ = Math.cbrt(l), m_ = Math.cbrt(m1), s_ = Math.cbrt(s1);
  const L = 0.2104542553 * l_ + 0.7936177850 * m_ - 0.0040720468 * s_;
  const a = 1.9779984951 * l_ - 2.4285922050 * m_ + 0.4505937099 * s_;
  const b1 = 0.0259040371 * l_ + 0.7827717662 * m_ - 0.8086757660 * s_; 

  const C = Math.sqrt(a * a + b1 * b1);
  const H = (Math.atan2(b1, a) * 180 / Math.PI + 360) % 360;

  return { l: L, c: C, h: H };
}

function oklchToCartesian(L, C, H) {
  const rad = H * Math.PI / 180;
  return { x: C * Math.cos(rad), y: C * Math.sin(rad), z: L };
}

function cartesianToOklch(x, y, z) {
  const c = Math.hypot(x, y);
  const h = (Math.atan2(y, x) * 180 / Math.PI + 360) % 360;
  return { l: z, c: c, h: h };
}

// --- Vector Operations ---
const vSub = (a, b) => ({ x: a.x - b.x, y: a.y - b.y, z: a.z - b.z });
const vAdd = (a, b) => ({ x: a.x + b.x, y: a.y + b.y, z: a.z + b.z });
const vScale = (v, s) => ({ x: v.x * s, y: v.y * s, z: v.z * s });
const vDot = (a, b) => a.x * b.x + a.y * b.y + a.z * b.z;
const vLen = (v) => Math.hypot(v.x, v.y, v.z);
const vNorm = (v) => {
  const l = vLen(v);
  return l === 0 ? { x: 0, y: 0, z: 0 } : vScale(v, 1 / l);
};
const vCross = (a, b) => ({
  x: a.y * b.z - a.z * b.y,
  y: a.z * b.x - a.x * b.z,
  z: a.x * b.y - a.y * b.x
});

// --- Spline Logic ---
function computeDefaultTangents(points) {
  const n = points.length;
  if (n < 2) return points.map(() => ({ x: 0, y: 0, z: 0 }));
  const tangents = new Array(n);
  for (let i = 0; i < n; i++) {
    if (i === 0) tangents[i] = vScale(vSub(points[1], points[0]), 0.5); 
    else if (i === n - 1) tangents[i] = vScale(vSub(points[n - 1], points[n - 2]), 0.5);
    else tangents[i] = vScale(vSub(points[i + 1], points[i - 1]), 0.5);
  }
  return tangents;
}

function hermite(p0, m0, p1, m1, t) {
  const t2 = t * t;
  const t3 = t2 * t;
  const h00 = 2 * t3 - 3 * t2 + 1;
  const h10 = t3 - 2 * t2 + t;
  const h01 = -2 * t3 + 3 * t2;
  const h11 = t3 - t2;
  return {
    x: h00 * p0.x + h10 * m0.x + h01 * p1.x + h11 * m1.x,
    y: h00 * p0.y + h10 * m0.y + h01 * p1.y + h11 * m1.y,
    z: h00 * p0.z + h10 * m0.z + h01 * p1.z + h11 * m1.z
  };
}

// --- Projection Logic ---
function getLocalBasis(points, index) {
  if (points.length < 3) {
    return {
      center: points[index] || { x: 0, y: 0, z: 0 },
      u: { x: 1, y: 0, z: 0 },
      v: { x: 0, y: 1, z: 0 }
    };
  }
  let pIndex = Math.max(1, Math.min(index, points.length - 2));
  const prev = points[pIndex - 1];
  const curr = points[pIndex];
  const next = points[pIndex + 1];

  let u = vNorm(vSub(next, prev));
  const v1 = vSub(curr, prev);
  const v2 = vSub(next, curr);
  let n = vCross(v1, v2);
  if (vLen(n) < 0.0001) {
    n = vCross(u, { x: 0, y: 0, z: 1 });
    if (vLen(n) < 0.0001) n = vCross(u, { x: 1, y: 0, z: 0 });
  }
  n = vNorm(n);
  const v = vNorm(vCross(n, u));
  const center = points[index];
  return { center, u, v };
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
        <div className="group relative flex items-center">
            {children}
            <div className={`
                absolute ${positionClasses[side]} 
                hidden group-hover:block z-50
                bg-gray-800 text-gray-200 text-[10px] font-medium
                px-2 py-1 rounded shadow-xl border border-gray-700
                whitespace-nowrap pointer-events-none
            `}>
                {text}
                <div className={`
                    absolute w-2 h-2 bg-gray-800 border-gray-700
                    ${side === 'top' ? 'bottom-[-5px] border-b border-r rotate-45 left-1/2 -translate-x-1/2' : ''}
                    ${side === 'bottom' ? 'top-[-5px] border-t border-l rotate-45 left-1/2 -translate-x-1/2' : ''}
                `} />
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
  width = 320, 
  height = 320 
}) {
  const svgRef = useRef(null);
  const [dragState, setDragState] = useState(null); 

  // 1. Basis
  const basis = useMemo(() => {
    const idx = viewIndex !== null ? viewIndex : Math.floor((points.length - 1) / 2);
    const safeIdx = Math.max(0, Math.min(idx, points.length - 1));
    return getLocalBasis(points, safeIdx);
  }, [points, viewIndex]);

  // 2. Scale (L2 Norm 3D)
  const scale = useMemo(() => {
    let maxDist = 0.05; 
    const centerIdx = points.indexOf(basis.center);
    const checkIndices = [];
    if(centerIdx > 0) checkIndices.push(centerIdx - 1);
    if(centerIdx < points.length - 1) checkIndices.push(centerIdx + 1);
    
    if(checkIndices.length === 0) return 150;

    checkIndices.forEach(i => {
        const dist = vLen(vSub(points[i], basis.center));
        maxDist = Math.max(maxDist, dist);
    });

    const viewportRadius = Math.min(width, height) / 2;
    const paddingFactor = 0.8; 
    return (viewportRadius * paddingFactor) / maxDist;
  }, [points, basis, width, height]);

  const offsetX = width / 2;
  const offsetY = height / 2;

  const toScreen = (pt) => {
    const rel = vSub(pt, basis.center);
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
    const handlePos3D = vAdd(pt, vScale(tangents[index], 1/3));
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

    const point = points[dragState.index];
    const newHandlePos = fromScreen(correctedX, correctedY);
    let delta = vSub(newHandlePos, point);

    if (mode === 'strength' && dragState.initialTangent) {
        const initT = dragState.initialTangent;
        const len = vLen(initT);
        if (len > 0.0001) {
            const dir = vScale(initT, 1/len);
            const projectedLen = vDot(delta, dir);
            delta = vScale(dir, projectedLen);
        }
    }

    const newTangent = vScale(delta, 3);
    const newTangents = [...tangents];
    newTangents[dragState.index] = newTangent;
    setTangents(newTangents);
  };

  const handlePointerUp = (e) => {
    if(dragState) {
        e.target.releasePointerCapture(e.pointerId);
        setDragState(null);
    }
  };

  const pathData = useMemo(() => {
    if (points.length < 2) return "";
    let d = "";
    const samples = 100;
    for (let i = 0; i < points.length - 1; i++) {
      const p0 = points[i];
      const p1 = points[i+1];
      const m0 = tangents[i];
      const m1 = tangents[i+1];
      for (let s = 0; s <= samples; s++) {
        const t = s / samples;
        const pt = hermite(p0, m0, p1, m1, t);
        const screenPt = toScreen(pt);
        if (i === 0 && s === 0) d += `M ${screenPt.x} ${screenPt.y}`;
        else d += ` L ${screenPt.x} ${screenPt.y}`;
      }
    }
    return d;
  }, [points, tangents, basis, scale]);

  return (
    <div className="bg-gray-900 rounded-lg shadow-inner inline-block select-none relative overflow-hidden border border-gray-800 group/editor">
       <div className="absolute top-3 left-3 text-[10px] text-gray-500 pointer-events-none z-10 font-mono">
         <div className="font-bold text-gray-300">SMART PROJECTION</div>
         <div>Aligned to Stop #{ (viewIndex ?? 0) + 1 }</div>
         <div className="text-[9px] opacity-50">Zoom: {scale.toFixed(0)}%</div>
       </div>

       <div className="absolute top-3 right-3 z-20 flex bg-gray-950 rounded border border-gray-800 p-0.5">
          <Tooltip text="Free Mode: Adjust angle & length" side="bottom">
              <button onClick={() => setMode('free')} className={`p-1.5 rounded ${mode === 'free' ? 'bg-blue-600 text-white' : 'text-gray-400 hover:text-gray-200'}`}>
                 <Unlock size={14} />
              </button>
          </Tooltip>
          <Tooltip text="Strength Mode: Adjust tension only" side="bottom">
              <button onClick={() => setMode('strength')} className={`p-1.5 rounded ${mode === 'strength' ? 'bg-blue-600 text-white' : 'text-gray-400 hover:text-gray-200'}`}>
                 <Lock size={14} />
              </button>
          </Tooltip>
       </div>

       <svg 
         ref={svgRef}
         width={width} 
         height={height} 
         className="block cursor-default"
         onPointerMove={handlePointerMove}
         onPointerUp={handlePointerUp}
       >
         <defs>
            <pattern id="grid" width="20" height="20" patternUnits="userSpaceOnUse">
                <path d="M 20 0 L 0 0 0 20" fill="none" stroke="#1f2937" strokeWidth="1"/>
            </pattern>
         </defs>
         <rect width="100%" height="100%" fill="url(#grid)" />
         <line x1={offsetX - 10} y1={offsetY} x2={offsetX + 10} y2={offsetY} stroke="#374151" />
         <line x1={offsetX} y1={offsetY - 10} x2={offsetX} y2={offsetY + 10} stroke="#374151" />
         <path d={pathData} stroke="white" strokeWidth="2" fill="none" strokeLinecap="round" />

         {points.map((pt, i) => {
            if (!tangents[i]) return null;
            const isActive = i === activeIndex;
            const screenPt = toScreen(pt);
            const handlePos = vAdd(pt, vScale(tangents[i], 1/3));
            const screenHandle = toScreen(handlePos);
            const backHandlePos = vSub(pt, vScale(tangents[i], 1/3));
            const screenBack = toScreen(backHandlePos);
            const isFocus = i === viewIndex;
            const opacity = isFocus ? 1 : (isActive ? 0.8 : 0.25);

            return (
              <g key={i} opacity={opacity} className="transition-opacity duration-200">
                <line x1={screenBack.x} y1={screenBack.y} x2={screenHandle.x} y2={screenHandle.y} stroke={isActive ? "#9ca3af" : "#4b5563"} strokeWidth="1" pointerEvents="none" />
                <circle 
                    cx={screenPt.x} cy={screenPt.y} r={4} 
                    fill={isActive ? "#fff" : "#9ca3af"} 
                    onPointerDown={(e) => { 
                        e.stopPropagation(); 
                        setActiveIndex(i); 
                        if(setViewIndex) setViewIndex(i);
                    }}
                    className="cursor-pointer hover:fill-white"
                >
                    <title>Stop #{i+1} (Click to Focus)</title>
                </circle>
                {isActive && mode === 'strength' && (
                     <line x1={screenPt.x} y1={screenPt.y} x2={screenPt.x + (screenHandle.x - screenPt.x) * 10} y2={screenPt.y + (screenHandle.y - screenPt.y) * 10} stroke="rgba(59, 130, 246, 0.2)" strokeDasharray="2 2" pointerEvents="none" />
                )}
                <circle 
                  cx={screenHandle.x} cy={screenHandle.y} r={6} 
                  fill={mode === 'strength' ? "#ec4899" : "#3b82f6"} 
                  stroke="white" strokeWidth={isActive ? 2 : 1}
                  className={`${isActive ? 'cursor-move' : 'cursor-pointer'} hover:scale-125`}
                  onPointerDown={(e) => handlePointerDown(e, i)}
                >
                    <title>Drag to adjust curvature</title>
                </circle>
              </g>
            );
         })}
       </svg>
    </div>
  );
}

function ColorPicker({ color, onChange, isActive, onActivate }) {
  const [h, s, b] = color;
  const hueRef = useRef(null);
  const sbRef = useRef(null);
  const handleHue = (e) => {
    const rect = hueRef.current.getBoundingClientRect();
    const y = Math.max(0, Math.min(1, (e.clientY - rect.top) / rect.height));
    onChange([y * 360, s, b]);
  };
  const handleSB = (e) => {
    const rect = sbRef.current.getBoundingClientRect();
    const x = Math.max(0, Math.min(1, (e.clientX - rect.left) / rect.width));
    const y = Math.max(0, Math.min(1, (e.clientY - rect.top) / rect.height));
    onChange([h, x, 1 - y]);
  };
  return (
    <div className={`flex gap-2 h-20 transition-all ${isActive ? 'opacity-100 scale-[1.02]' : 'opacity-60 hover:opacity-90'}`} onClick={onActivate}>
      <div ref={hueRef} className="relative w-6 h-full rounded overflow-hidden cursor-pointer shadow-sm touch-none ring-1 ring-white/10" style={{ background: 'linear-gradient(to bottom, #f00 0%, #ff0 17%, #0f0 33%, #0ff 50%, #00f 67%, #f0f 83%, #f00 100%)' }}
        onPointerDown={(e) => { e.target.setPointerCapture(e.pointerId); onActivate(); handleHue(e); }}
        onPointerMove={(e) => e.buttons === 1 && handleHue(e)}
        onPointerUp={(e) => e.target.releasePointerCapture(e.pointerId)}>
        <div className="absolute left-0 right-0 h-0.5 border border-white bg-black/20" style={{ top: `${h / 360 * 100}%` }} />
      </div>
      <div ref={sbRef} className="relative flex-1 h-full rounded overflow-hidden cursor-crosshair shadow-sm touch-none ring-1 ring-white/10" style={{ background: `linear-gradient(to top, black, transparent), linear-gradient(to right, white, hsl(${h}, 100%, 50%))` }}
        onPointerDown={(e) => { e.target.setPointerCapture(e.pointerId); onActivate(); handleSB(e); }}
        onPointerMove={(e) => e.buttons === 1 && handleSB(e)}
        onPointerUp={(e) => e.target.releasePointerCapture(e.pointerId)}>
        <div className="absolute w-3 h-3 border-2 border-white rounded-full shadow-md" style={{ left: `${s * 100}%`, top: `${(1 - b) * 100}%`, transform: 'translate(-50%, -50%)' }} />
      </div>
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
            <div className="flex justify-between items-center px-3 py-2 bg-gray-800/50 border-b border-gray-800">
                <span className="text-[10px] font-bold uppercase tracking-widest text-gray-400 flex items-center gap-2">
                    {label.includes('Strength') ? <Braces size={12}/> : <FileCode size={12}/>} 
                    {label}
                </span>
                <Tooltip text="Copy content" side="left">
                    <button onClick={handleCopy} className="text-gray-400 hover:text-white transition-colors">
                        {copied ? <Check size={14} className="text-green-400"/> : <Copy size={14}/>}
                    </button>
                </Tooltip>
            </div>
            <textarea 
                readOnly 
                value={value} 
                className="w-full h-24 bg-transparent p-3 text-[10px] font-mono text-gray-400 resize-none focus:outline-none focus:text-gray-200"
                onClick={(e) => e.target.select()}
            />
        </div>
    );
}

// ==========================================
// 3. Main Application
// ==========================================

export default function App() {
  const [colors, setColors] = useState([ [0, 1, 1], [120, 1, 1], [240, 1, 1] ]);
  const [activePointIndex, setActivePointIndex] = useState(1);
  const [viewPointIndex, setViewPointIndex] = useState(1);
  const [tangentsState, setTangentsState] = useState([]);
  const [isManual, setIsManual] = useState(false);
  const [position, setPosition] = useState(0.5); 
  const [steps, setSteps] = useState(40); 
  const [editorMode, setEditorMode] = useState('free');
  const [copied, setCopied] = useState(false);

  const points = useMemo(() => {
    return colors.map(([h, s, b]) => {
      const { l, c, h: oh } = hsbToOklch(h, s, b);
      return oklchToCartesian(l, c, oh);
    });
  }, [colors]);

  useEffect(() => {
    setTangentsState([]);
    setIsManual(false);
  }, [points.length]);

  const effectiveTangents = useMemo(() => {
    if (tangentsState.length === points.length) return tangentsState;
    return computeDefaultTangents(points);
  }, [points, tangentsState]);

  const handleTangentUpdate = (newT) => {
    setTangentsState(newT);
    setIsManual(true);
  };

  // --- Calculate Tangent Strengths (Ratio vs Default) ---
  const tangentStrengthList = useMemo(() => {
      const defaults = computeDefaultTangents(points);
      const strengths = effectiveTangents.map((t, i) => {
          const defLen = vLen(defaults[i]);
          const curLen = vLen(t);
          // If default is 0 (collinear), fallback to 1, else ratio
          if (defLen < 0.00001) return 1.0;
          return parseFloat((curLen / defLen).toFixed(3));
      });
      return JSON.stringify(strengths, null, 0).replace(/,/g, ', ');
  }, [points, effectiveTangents]);


  const { polyGradient, linearGradient, polyVal, linearVal, polyRgb, linearRgb } = useMemo(() => {
    if (points.length < 2) return { polyGradient:'', linearGradient: '', polyVal: null, linearVal: null };
    
    // --- Polynomial ---
    const stops = [];
    const n = points.length;
    for (let i = 0; i <= steps; i++) {
        const globalT = i / steps;
        const scaled = globalT * (n - 1);
        let segIdx = Math.floor(scaled);
        if (segIdx >= n - 1) segIdx = n - 2;
        const localT = scaled - segIdx;
        const pt = hermite(points[segIdx], effectiveTangents[segIdx], points[segIdx + 1], effectiveTangents[segIdx + 1], localT);
        const { l, c, h } = cartesianToOklch(pt.x, pt.y, pt.z);
        stops.push(`oklch(${l.toFixed(4)} ${c.toFixed(4)} ${h.toFixed(2)}deg) ${(globalT * 100).toFixed(1)}%`);
    }
    const polyGradient = `linear-gradient(to right, ${stops.join(', ')})`;

    // --- Linear ---
    const linStops = colors.map((_, i) => {
       const p = points[i];
       const { l, c, h } = cartesianToOklch(p.x, p.y, p.z);
       return `oklch(${l.toFixed(4)} ${c.toFixed(4)} ${h.toFixed(2)}deg) ${(i / (points.length-1) * 100).toFixed(1)}%`;
    });
    const linearGradient = `linear-gradient(to right, ${linStops.join(', ')})`;

    // --- Values at Position ---
    const scaledP = position * (n - 1);
    let segP = Math.floor(scaledP); if (segP >= n - 1) segP = n - 2;
    const uP = scaledP - segP;
    const polyPt = hermite(points[segP], effectiveTangents[segP], points[segP + 1], effectiveTangents[segP + 1], uP);
    const polyVal = cartesianToOklch(polyPt.x, polyPt.y, polyPt.z);
    const polyRgb = oklchToRgb(polyVal.l, polyVal.c, polyVal.h);

    let segL = Math.floor(scaledP); if(segL >= n-1) segL = n-2;
    const uL = scaledP - segL;
    const pA = points[segL]; const pB = points[segL+1];
    const linPt = { x: pA.x + (pB.x - pA.x) * uL, y: pA.y + (pB.y - pA.y) * uL, z: pA.z + (pB.z - pA.z) * uL };
    const linearVal = cartesianToOklch(linPt.x, linPt.y, linPt.z);
    const linearRgb = oklchToRgb(linearVal.l, linearVal.c, linearVal.h);

    return { polyGradient, linearGradient, polyVal, linearVal, polyRgb, linearRgb };
  }, [points, effectiveTangents, colors, position, steps]);

  const copyToClipboard = () => {
    navigator.clipboard.writeText(polyGradient);
    setCopied(true);
    setTimeout(() => setCopied(false), 2000);
  };

  return (
    <div className="min-h-screen bg-gray-950 text-gray-200 py-8 px-4 font-sans flex justify-center">
      <div className="w-full max-w-6xl flex flex-col gap-6">
        
        {/* Header */}
        <header className="flex justify-between items-end border-b border-gray-800 pb-4">
          <div>
            <h1 className="text-2xl font-bold text-white mb-1 tracking-tight flex items-center gap-2">
               <GitCommit className="text-blue-500" />
               OKLCH Spline Gradient
            </h1>
            <p className="text-sm text-gray-400">Generates smooth Catmull-Rom gradients in OKLCH space.</p>
          </div>
          <div className="flex gap-2">
            {isManual && (
                <Tooltip text="Reset to default smooth curves" side="bottom">
                    <button onClick={() => { setTangentsState([]); setIsManual(false); }} className="flex items-center gap-2 text-xs font-medium bg-gray-800 hover:bg-gray-700 px-3 py-1.5 rounded text-blue-300 transition-colors">
                        <RotateCcw size={14} /> Reset Curve
                    </button>
                </Tooltip>
            )}
            <Tooltip text="Copy the Spline Gradient CSS" side="bottom">
                <button onClick={copyToClipboard} className="flex items-center gap-2 text-xs font-medium bg-blue-600 hover:bg-blue-500 px-3 py-1.5 rounded text-white transition-colors">
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
                             Spline 
                             <Tooltip text="Gradient Resolution (Steps)" side="right">
                                <input type="range" min="2" max="100" value={steps} onChange={(e) => setSteps(Number(e.target.value))} className="w-16 accent-blue-500 cursor-pointer"/>
                             </Tooltip>
                             <span className="text-gray-600 font-mono normal-case">{steps} stops</span>
                        </span>
                    </div>
                    
                    <div className="relative group select-none">
                        <div className="h-16 w-full rounded-t-lg shadow-lg ring-1 ring-white/10 z-0" style={{ background: polyGradient }} />
                        <div className="h-8 w-full rounded-b-lg ring-1 ring-white/10 z-0 border-t border-black/20" style={{ background: linearGradient }} />
                        <div className="absolute top-0 bottom-0 z-10 pointer-events-none flex flex-col items-center justify-center" style={{ left: `${position * 100}%`, transform: 'translateX(-50%)' }}>
                            <div className="w-0.5 h-full bg-white shadow-[0_0_10px_rgba(0,0,0,0.5)]" />
                            <div className="absolute top-1/2 -translate-y-1/2 bg-white text-gray-950 text-[10px] font-bold px-1.5 py-0.5 rounded-full shadow-md whitespace-nowrap">{(position * 100).toFixed(1)}%</div>
                        </div>
                        <Tooltip text="Scrub to compare values" side="top">
                            <input type="range" min="0" max="1" step="0.001" value={position} onChange={(e) => setPosition(parseFloat(e.target.value))} className="absolute inset-0 w-full h-full opacity-0 cursor-ew-resize z-20"/>
                        </Tooltip>
                    </div>
                    <div className="flex justify-between text-[10px] font-bold uppercase tracking-widest text-gray-500 px-1 pt-1"><span>Linear (Reference)</span></div>
                </div>

                {polyVal && (
                <div className="grid grid-cols-2 gap-4">
                     {/* Spline Stats */}
                     <div className="bg-gray-900 p-3 rounded border border-gray-800">
                        <div className="text-[10px] text-blue-400 mb-1 uppercase tracking-wider font-bold">Spline</div>
                        <div className="flex items-start gap-3">
                             <div className="w-10 h-10 rounded shadow-sm ring-1 ring-white/10" style={{ background: `oklch(${polyVal.l} ${polyVal.c} ${polyVal.h}deg)` }} />
                             <div className="space-y-2 flex-1">
                                 <div className="text-xs font-mono text-gray-300">
                                    <div className="flex justify-between"><span>L</span> <span>{polyVal.l.toFixed(3)}</span></div>
                                    <div className="flex justify-between"><span>C</span> <span>{polyVal.c.toFixed(3)}</span></div>
                                    <div className="flex justify-between"><span>H</span> <span>{polyVal.h.toFixed(1)}</span></div>
                                 </div>
                                 <div className="border-t border-gray-800 pt-1 text-[10px] font-mono text-gray-500 flex justify-between items-center">
                                    <Monitor size={10} /><span>{polyRgb.r}, {polyRgb.g}, {polyRgb.b}</span>
                                 </div>
                             </div>
                        </div>
                     </div>
                     {/* Linear Stats */}
                     <div className="bg-gray-900 p-3 rounded border border-gray-800">
                        <div className="text-[10px] text-gray-500 mb-1 uppercase tracking-wider font-bold">Linear</div>
                        <div className="flex items-start gap-3">
                             <div className="w-10 h-10 rounded shadow-sm ring-1 ring-white/10" style={{ background: `oklch(${linearVal.l} ${linearVal.c} ${linearVal.h}deg)` }} />
                             <div className="space-y-2 flex-1">
                                 <div className="text-xs font-mono text-gray-300">
                                    <div className="flex justify-between"><span>L</span> <span>{linearVal.l.toFixed(3)}</span></div>
                                    <div className="flex justify-between"><span>C</span> <span>{linearVal.c.toFixed(3)}</span></div>
                                    <div className="flex justify-between"><span>H</span> <span>{linearVal.h.toFixed(1)}</span></div>
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
                      <span>3D Path Editor</span>
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
                          width={320} 
                          height={250} 
                      />
                 </div>
            </div>
        </div>

        {/* Color Stops */}
        <div className="border-t border-gray-800 pt-6">
            <div className="flex justify-between items-center mb-4">
                <h2 className="text-sm font-semibold uppercase tracking-wider text-gray-400">Color Stops</h2>
                <Tooltip text="Add a new color stop" side="left">
                    <button onClick={() => { const newIdx = colors.length; setColors([...colors, [Math.random() * 360, 1, 1]]); setActivePointIndex(newIdx); setViewPointIndex(newIdx); }} className="bg-blue-600 hover:bg-blue-500 text-white p-2 rounded shadow-lg shadow-blue-900/20 transition-all active:scale-95">
                        <Plus size={18} />
                    </button>
                </Tooltip>
            </div>
            <div className="grid grid-cols-2 md:grid-cols-4 lg:grid-cols-5 gap-4">
                {colors.map((color, i) => (
                    <div key={i} className={`bg-gray-900 p-3 rounded-lg border relative transition-all duration-200 group ${activePointIndex === i ? 'border-blue-500/50 shadow-md shadow-blue-900/10' : 'border-gray-800 hover:border-gray-700'}`} onClick={() => { setActivePointIndex(i); setViewPointIndex(i); }}>
                        <div className="flex justify-between items-center mb-2">
                            <span className={`font-mono text-[10px] px-1.5 py-0.5 rounded ${activePointIndex === i ? 'bg-blue-900/30 text-blue-200' : 'bg-gray-800 text-gray-500'}`}>#{i + 1}</span>
                            {colors.length > 2 && (
                                <button onClick={(e) => { e.stopPropagation(); setColors(colors.filter((_, idx) => idx !== i)); if(activePointIndex === i) { setActivePointIndex(0); setViewPointIndex(0); } }} className="text-gray-600 hover:text-red-400 p-1 opacity-0 group-hover:opacity-100 transition-opacity"><X size={14} /></button>
                            )}
                        </div>
                        <ColorPicker color={color} onChange={(c) => { const newColors = [...colors]; newColors[i] = c; setColors(newColors); }} isActive={activePointIndex === i} onActivate={() => { setActivePointIndex(i); setViewPointIndex(i); }} />
                    </div>
                ))}
            </div>
        </div>

        {/* Output Area */}
        <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <CodeBlock label="Spline Gradient (CSS)" value={polyGradient} />
            <CodeBlock label="Tangent Strengths (JSON)" value={tangentStrengthList} />
            <CodeBlock label="Linear Gradient (CSS)" value={linearGradient} />
        </div>

      </div>
    </div>
  );
}