/*jshint esversion: 6 */
// @ts-check
/* global MathJax */


// --- up-front typedefs ----
/**
 * An exact interval [u0,u1].
 * @typedef {Rational[]} Interval
 */

/**
 * One piece of a basis function:
 *  - poly: an array of exact Rational coefficients
 *  - interval: the closed interval on which this polynomial is valid
 * @typedef {{ poly: Rational[], interval: Interval }} Segment
 */

import * as T from "../../libs/CS559-Three/build/three.module.js";
import { GrWorld } from "../../libs/CS559-Framework/GrWorld.js";
import { OrbitControls } from "../../libs/CS559-Three/examples/jsm/controls/OrbitControls.js";
import { draggablePoints } from "../../libs/CS559/dragPoints.js";
import { RunCanvas } from "../../libs/CS559/runCanvas.js";

// A tiny exact rational type
class Rational {
    constructor(n, d = 1) {
        if (d < 0) { n = -n;  d = -d; }
        const g = Rational.gcd(Math.abs(n), d);
        this.n = n / g;
        this.d = d / g;
    }
    static gcd(a, b) {
        // make sure we’re working with plain non-negative ints
        a = Math.abs(Math.floor(a));
        b = Math.abs(Math.floor(b));
        while (b !== 0) {
          const t = b;
          b = a % b;
          a = t;
        }
        return a;
    }
    add(r) {
        return new Rational(
        this.n * r.d + r.n * this.d,
        this.d * r.d
        );
    }
    // subtract r: (n1/d1 - n2/d2)
    sub(r) {
        return new Rational(
        this.n * r.d - r.n * this.d,
        this.d * r.d
        );
    }
    mul(r) {
        return new Rational(this.n * r.n, this.d * r.d);
    }
    inv() {
        return new Rational(this.d, this.n);
    }
    div(r) {
        return this.mul(r.inv());
    }
    neg() {
        return new Rational(-this.n, this.d);
    }
    // allow JS to coerce Rational → Number
    valueOf() {
        return this.n / this.d;
    }
      // compare this ? r  ==> 1 if >, 0 if =, -1 if <
    compareTo(r) {
        const lhs = this.n * r.d, rhs = r.n * this.d;
        return lhs < rhs ? -1 : lhs > rhs ? 1 : 0;
    }
    static max(a, b) { return a.compareTo(b) >= 0 ? a : b; }
    static min(a, b) { return a.compareTo(b) <= 0 ? a : b; }
    // toLatex: "0", "3", "\\tfrac23", etc.
    toLatex() {
        if (this.n === 0) return "0";
        if (this.d === 1) return this.n.toString();
        // if ±1/1 handled, here proper fraction
        // If numerator larger than denom, let LaTeX convert to mixed number?
        return `\\tfrac{${this.n}}{${this.d}}`;
    }
    // strip sign for coefficient when it's exactly ±1
    isOne() { return this.n === this.d; }
    isMinusOne() { return this.n === -this.d; }
}

/**
 * Parse a JS number (possibly with decimals) into an exact Rational.
 *   0.125  → new Rational(1, 8)
 *   2.75   → new Rational(11, 4)
 */
function toRational(x) {
    const s = x.toString();
    if (!s.includes('.')) {
        // integer
        return new Rational(Number(s), 1);
    }
    const [intPart, fracPart] = s.split('.');
    const d = fracPart.length;
    const numer = Number(intPart) * 10**d + Number(fracPart);
    const denom = 10**d;
    return new Rational(numer, denom);
}
  
  
function polyAdd(a, b) {
    const m = Math.max(a.length, b.length), r = [];
    for (let i = 0; i < m; i++) {
    const ai = a[i] || new Rational(0),
            bi = b[i] || new Rational(0);
    r[i] = ai.add(bi);
    }
    return r;
}

function polyMul(a, b) {
    const N = a.length + b.length - 1;
    // make N copies of new Rational(0):
    const r = Array.from({ length: N }, () => new Rational(0));
    for (let i = 0; i < a.length; i++) {
        for (let j = 0; j < b.length; j++) {
        r[i + j] = r[i + j].add(a[i].mul(b[j]));
        }
    }
    return r;
}
  
  

// --- merge helper ---
/**
 * @param {Segment[]} segs
 * @returns {Segment[]}
 */
function mergeSegments(segs) {
    // Map from “u0|u1” → the merged Segment
    /** @type {Map<string, Segment>} */
    const map = new Map();
  
    for (let seg of segs) {
      const [u0, u1] = seg.interval;
      // key by exact Rational string (you could also do `${u0.n}/${u0.d}|${u1.n}/${u1.d}`)
      const key = u0.toLatex() + "|" + u1.toLatex();
  
      if (map.has(key)) {
        // merge into the existing segment
        const existing = map.get(key);
        existing.poly = polyAdd(existing.poly, seg.poly);
      } else {
        // clone the poly *and* reuse the exact same interval tuple
        map.set(key, {
          poly: seg.poly.slice(),
          // cast this exact 2-element Rational array to our Interval typedef
          interval: /** @type {Interval} */ ([u0, u1])
        });
      }
    }
  
    // grab all the merged Segments out of the map
    const out = Array.from(map.values());
  
    // sort by the numeric value of u0
    out.sort((A, B) => {
      const a0 = A.interval[0], b0 = B.interval[0];
      // convert Rational → Number for sort
      return (a0.n / a0.d) - (b0.n / b0.d);
    });
  
    return out;
}

// --- Cox-de Boor recursion with merging ---
/**
 * compute the blending (B-spline) basis B_{k,d}(u) on a knot vector U
 *
 * @param {Rational[]} U    an array of Rational knots (must be non-decreasing)
 * @param {number}     k    the basis index
 * @param {number}     d    the spline order (degree+1)
 * @returns {Segment[]}     a list of piecewise segments for B_{k,d}(u)
 */
function blendingFunction(U, k, d) {
    // base case d = 0
    if (d === 1) {
        const u0 = U[k];
        const u1 = U[k + 1];
        if (u0 < u1) return [{ poly: [new Rational(1)], interval: [u0, u1] }];
        else return [];
    }
    const segs = [];
    // left term
    const denL = U[k + d - 1].sub(U[k]);
    if (denL.valueOf() > 0) {
        const Lfac = [
            U[k].neg().div(denL),
            new Rational(1).div(denL)
        ]
        // const Lfac = [-U[k]/denL, 1 / denL];
        for (let seg of blendingFunction(U, k, d - 1)) {
            const lo = Rational.max(seg.interval[0], U[k]);
            const hi = Rational.min(seg.interval[1], U[k + d - 1]);
            if (lo.compareTo(hi) < 0) {
                segs.push({
                    poly: polyMul(Lfac, seg.poly),
                    interval: [lo, hi]
                });
            }
        }
    }
    // right term
    const denR = U[k + d].sub(U[k + 1]);
    // const denR = U[k + d] - U[k + 1];
    if (denR.valueOf() > 0) {

        const Rfrac = [
            U[k + d].div(denR),
            new Rational(-1).div(denR)
        ]
        // const Rfrac = [U[k + d] / denR, -1 / denR];
        for (let seg of blendingFunction(U, k + 1, d - 1)) {
            const lo = Rational.max(seg.interval[0], U[k + 1]);
            const hi = Rational.min(seg.interval[1], U[k + d]);
            if (lo.compareTo(hi) < 0) {
                segs.push({
                    poly: polyMul(Rfrac, seg.poly),
                    interval: [lo, hi]
                });
            }
        }
    }
    // merge any identical interval pieces
    return mergeSegments(segs);
}

/**
 * @param {Rational[]} U knot vector (non-decreasing)
 * @param {number} d order
 * @returns {Array<Array<{ poly: number[], interval: number[]}>>} segs
 *    an array B of length U.length - d, where
 *    B[i] is the piecewise representation of B_{i,d}(u)
*/
function calc_BlendingFunction(U, d) {
    const n = U.length - d;
    const B = new Array(n);
    for (let i = 0; i < n; i++) {
        B[i] = blendingFunction(U, i, d);
    }
    return B;
}

// /**
//  * Turn [a0, a1, a2, …] into a string like
//  *   "1 - 3u + 5u^2"
//  */
// function polyToLatex(coefs) {
//     const terms = [];
//     for (let k = 0; k < coefs.length; k++) {
//       const a = coefs[k];
//       if (Math.abs(a) < 1e-8) continue;           // skip zeros
//       const sign = a < 0 ? " - " : (terms.length ? " + " : "");
//       const absA = Math.abs(a);
//       let piece;
//       if (k === 0) {
//         piece = absA.toFixed(5).replace(/\.?0+$/, "");
//       } else {
//         const coeff =
//           Math.abs(absA - 1) < 1e-8 ? "" : absA.toFixed(5).replace(/\.?0+$/, "");
//         const power = k === 1 ? "u" : `u^{${k}}`;
//         piece = coeff + power;
//       }
//       terms.push(sign + piece);
//     }
//     return terms.join("") || "0";
//   }

// // Now render *one* big piecewise LaTeX string:
// function piecewiseToLatex(Bs, d) {
//     // Bs is an array of segment‐arrays for each i.
//     // We’ll produce one <div> per basis function.
//     const container = document.getElementById("blend-formulas");
//     container.innerHTML = "";
  
//     Bs.forEach((segs, i) => {
//       // start the cases for B_{i,d}(u)
//       let tex = `\\[ B_{${i},${d}}(u)=\\begin{cases}`;
  
//       // first: “0” outside the support
//       const umin = segs[0].interval[0],
//             umax = segs[segs.length - 1].interval[1];
//       tex += `0,&u<${umin.valueOf()}\\text{ or }u>${umax.valueOf()}\\\\`;
  
//       // then each active piece
//       segs.forEach(({ poly, interval: [u0, u1] }) => {
//         const terms = [];
  
//         for (let k = poly.length - 1; k >= 0; k--) {
//           const r = poly[k];
//           if (r.n === 0) continue;
  
//           // determine sign
//           const sign = (terms.length === 0)
//             ? (r.n < 0 ? "-" : "")
//             : (r.n < 0 ? " - " : " + ");
  
//           // **your original rational coefficient**  
//           const ar = new Rational(Math.abs(r.n), r.d);
  
//           // pull exponent into its own string to avoid nested backticks
//           const exp = (k === 1 ? "" : `^{${k}}`);
  
//           let piece;
//           if (k === 0) {
//             // constant term uses ar.toLatex()
//             piece = ar.toLatex();
//           } else if (ar.isOne()) {
//             // coefficient  ±1 just gives u^k
//             piece = `u${exp}`;
//           } else {
//             // general coefficient * u^k
//             piece = ar.toLatex() + `u${exp}`;
//           }
  
//           terms.push(sign + piece);
//         }
  
//         const polyLatex = terms.join("") || "0";
//         tex += `${polyLatex},& ${u0.valueOf()}\\le u\\le ${u1.valueOf()}\\\\`;
//       });
  
//       tex += `\\end{cases}\\]`;
  
//       const div = document.createElement("div");
//       div.innerHTML = tex;
//       container.appendChild(div);
//     });
  
//     // re-typeset all the new math
//     MathJax.typesetPromise();
// }

function piecewiseWithPlots(Bs, d) {
    const container = document.getElementById("blend-formulas");
    container.innerHTML = "";
  
    const rows = []; // to remember canvases and which i they belong to
    let globalUmin = Bs[0][0].interval[0];
    let globalUmax = Bs[Bs.length - 1][Bs[Bs.length - 1].length - 1].interval[1];
  
    Bs.forEach((segs, i) => {
      // 1. Figure out support [umin,umax]
      const umin = segs[0].interval[0].valueOf();
      const umax = segs[segs.length - 1].interval[1].valueOf();
  
      // 2. Build the LaTeX for B_{i,d}(u)
      let tex = `\\[ B_{${i},${d}}(u)=\\begin{cases}`;
      tex += `0,&u<${umin}\\text{ or }u>${umax}\\\\`;
  
      segs.forEach(({ poly, interval:[u0,u1] }) => {
        // build the “a_k u^k” pieces exactly as before…
        const terms = [];
        for (let k=poly.length-1; k>=0; k--) {
          const r = poly[k];
          if (r.n===0) continue;
          const sign = terms.length===0
            ? (r.n<0?"-":"")
            : (r.n<0?" - ":" + ");
          const ar = new Rational(Math.abs(r.n), r.d);
          const exp = k===1?"":`^{${k}}`;
          let piece;
          if (k===0)             piece = ar.toLatex();
          else if (ar.isOne())   piece = `u${exp}`;
          else                   piece = `${ar.toLatex()}u${exp}`;
          terms.push(sign+piece);
        }
        const P = terms.join("")||"0";
        tex += `${P},&${u0.valueOf()}\\le u\\le ${u1.valueOf()}\\\\`;
      });
  
      tex += `\\end{cases}\\]`;
  
      // 3. Create a flex-row container
      const row = document.createElement("div");
      row.style.display = "flex";
      row.style.alignItems = "center";
      row.style.marginBottom = "1em";
  
      // formula wrapper
      const mathDiv = document.createElement("div");
      if (degree <= 6) {
        mathDiv.style.minWidth = "600px";
        mathDiv.style.maxWidth = "600px";
      } else {
        mathDiv.style.minWidth = "800px";
        mathDiv.style.maxWidth = "800px";
      }

      mathDiv.innerHTML = tex;
      row.appendChild(mathDiv);
  
      // canvas for the plot
      const cv = document.createElement("canvas");
      cv.width = 400;
      cv.height = 100;
      cv.style.border = "1px solid #ccc";
      cv.style.marginLeft = "0.5em";
      row.appendChild(cv);
  
      container.appendChild(row);
  
      // remember for plotting
      rows.push({ canvas: cv, i, umin, umax });
    });
  
    // 4. Typeset everything, then draw plots
    // MathJax.typesetPromise().then(() => {
    //   rows.forEach(({ canvas, i, umin, umax }) => {
    //     const ctx = canvas.getContext("2d");
    //     ctx.clearRect(0,0,canvas.width,canvas.height);
    //     ctx.beginPath();
  
    //     const N = 200;
    //     for (let j=0; j<=N; j++) {
    //     //   const u = umin + (umax-umin)*(j/N);
    //         const u = globalUmin.valueOf() + (globalUmax.valueOf()-globalUmin.valueOf())*(j/N);
    //       const b = evalBasis(i, u);
    //       // map (u,b) → (x,y) in pixel coords
    //     //   const x = (u-umin)/(umax-umin)*canvas.width;
    //         const x = (u-globalUmin.valueOf())/(globalUmax.valueOf()-globalUmin.valueOf())*canvas.width;
    //       // flip y (we want b=0 at bottom)
    //       const y = canvas.height - b*canvas.height;
    //       if (j===0) ctx.moveTo(x,y);
    //       else       ctx.lineTo(x,y);
    //     }
  
    //     ctx.strokeStyle = "steelblue";
    //     ctx.lineWidth = 2;
    //     ctx.stroke();
    //   });
    // });
    MathJax.typesetPromise().then(() => {
        rows.forEach(({ canvas, i, umin, umax }) => {
          const ctx = canvas.getContext("2d");
          const W  = canvas.width, H = canvas.height;
          ctx.clearRect(0,0,W,H);
      
          // 1. Draw the axes
            ctx.beginPath();
            ctx.strokeStyle = "black";
            ctx.lineWidth   = 1;
            ctx.moveTo(0, H);     // x-axis
            ctx.lineTo(W, H);
            ctx.moveTo(0, 0);     // y-axis
            ctx.lineTo(0, H);
            ctx.stroke();

            // --- 2. draw ticks & labels ---
            let numXTicks = globalUmax + 1, numYTicks = 2;
            if (numXTicks > 10) {
                numXTicks = 10;
            }
            if (numXTicks < 2) {
                numXTicks = 2;
            }

            // make sure text is visible
            ctx.fillStyle    = "black";
            ctx.strokeStyle  = "black";
            ctx.font         = "10px sans-serif";

            // — x-ticks & labels —
            ctx.textAlign    = "center";
            ctx.textBaseline = "top";
            for (let t = 0; t <= numXTicks; t++) {
                const α    = t/numXTicks;
                const x_px = α * W;
                let u_val = 0;
                if(globalUmax < 2) {
                    u_val = globalUmin + α*(globalUmax - globalUmin);
                } else {
                    u_val = globalUmin + α*(numXTicks - globalUmin);
                }

                // tick
                ctx.beginPath();
                ctx.moveTo(x_px, H);
                ctx.lineTo(x_px, H-6);
                ctx.stroke();

                // label (slightly inside the bottom border)
                if (globalUmax < 2) {
                    ctx.fillText(u_val.toFixed(1), x_px + 2, H - 6 - 4);
                }
                else {
                    ctx.fillText(u_val.toFixed(0), x_px + 3, H - 6 - 4);
                }
            }

            // — y-ticks & labels —
            ctx.textAlign    = "right";
            ctx.textBaseline = "middle";
            for (let t = 0; t <= numYTicks; t++) {
                const β    = t/numYTicks,
                        y_px = H - β*H,
                        bVal = β;            // 0..1

                // tick
                ctx.beginPath();
                ctx.moveTo(0, y_px);
                ctx.lineTo(6, y_px);
                ctx.stroke();

                // label (pulled in 2px from the left edge)
                ctx.fillText(bVal.toFixed(1), 6 + 8, y_px + 8);
            }

      
          // --- 3. draw the blending‐curve ---
          let effectiveFlag = 0;
          if (evalBasis(i, highlightU) > 0) {
            effectiveFlag = 1;
          }
        //   ctx.beginPath();
          const N = 400;
          let flag = 0;
          for (let j = 0; j <= N - 1; j++) {
            const alpha = j/N;
            const alpha1 = (j+1)/N;
            let u = globalUmin.valueOf() + (globalUmax.valueOf()-globalUmin.valueOf())*alpha;
            let u1 = globalUmin.valueOf() + (globalUmax.valueOf()-globalUmin.valueOf())*(j+1)/N;
            let x = 0;
            let x1 = 0;
            if (globalUmax < 2) {
                if (highlightUEnabled && u + 0.01 * (globalUmax - globalUmin) > highlightU && u - 0.01 * (globalUmax - globalUmin) < highlightU) {
                    ctx.strokeStyle = "red";
                    ctx.lineWidth = 5;
                } else {
                    ctx.strokeStyle = effectiveFlag && highlightCPandBFEnabled ? "#FF8800" : (highlightCPandBFEnabled? "lightgray" : "steelblue");
                    ctx.lineWidth = 2;
                }
                x = alpha * W;
                x1 = (j+1)/N * W;
            } else {
                u = globalUmin.valueOf() + (globalUmax.valueOf() + 1 -globalUmin.valueOf())*alpha;
                u1 = globalUmin.valueOf() + (globalUmax.valueOf() + 1 - globalUmin.valueOf())*(j+1)/N;
                if (highlightUEnabled && u + 0.01 * (globalUmax - globalUmin) > highlightU && u - 0.01 * (globalUmax - globalUmin) < highlightU) {
                    ctx.strokeStyle = "red";
                    ctx.lineWidth = 5;
                } else {
                    ctx.strokeStyle = effectiveFlag && highlightCPandBFEnabled ? "#FF8800" : (highlightCPandBFEnabled? "lightgray" : "steelblue");
                    ctx.lineWidth = 2;
                }
                let old_u = globalUmin.valueOf() + (globalUmax.valueOf() + 1 - globalUmin.valueOf())* (j-1)/N;
                if(u1 > umax && u <= umax) {
                    flag = 1;
                }
                x = alpha * W;
                x1 = (j+1)/N * W;
            }
            const b = evalBasis(i, u);
            const y = H - b * H;
            const y1 = H - evalBasis(i, u1) * H;
            // if (j === 0 || flag === 1) ctx.moveTo(x,y);
            // else         ctx.lineTo(x,y);
            if(flag === 1) {
                flag = 0;
            } else {
                ctx.beginPath();
                ctx.moveTo(x, y);
                ctx.lineTo(x1, y1);
                ctx.stroke();
            }
          }
        //   ctx.strokeStyle = "steelblue";
        //   ctx.lineWidth   = 2;
        //   ctx.stroke();
        });
      });
      
}
  
  

// 2. Evaluate one basis B_i at u
function evalBasis(i, u) {
    for (let { poly, interval: [u0, u1] } of Bs[i]) {
        if( u >= u0.valueOf() && u <= u1.valueOf()) {
            return poly.reduce((sum, a_k, k) => sum + a_k.valueOf() * u**k, 0);
        }
    }
    return 0;
}

// 3. Evaluate the curve C(u) = sum_i B_i(u) P_i
function evalCurve(u) {
    let  x = 0, y = 0;
    for (let i = 0; i < thePoints.length; i++) {
        const b = evalBasis(i, u);
        x += b * thePoints[i][0];
        y += b * thePoints[i][1];
    }
    return [x, y];
}

// 4. Draw the curve inside wrapDraw (after drawing control points)
// function drawBSpline() {
//     uMin = U[order - 1].valueOf();
//     uMax = U[U.length - order].valueOf();
//     const step = (uMax - uMin) / 400;
//     let prev = null;    // { u, x, y }
    
//     for (let u = uMin; u <= uMax + 1e-9; u += step) {
//       const [x, y] = evalCurve(u);
//       if (prev) {
//         // decide whether this little piece crosses the highlightU
//         const inSegment =
//           highlightUEnabled && 
//           highlightU + step * 4 >= prev.u &&
//           highlightU - step * 4 <= u;
        
//         // draw that tiny piece
//         context.beginPath();
//         context.moveTo(prev.x, prev.y);
//         context.lineTo(x, y);
//         context.strokeStyle = inSegment ? "red" : "#0000AA";
//         context.lineWidth   = inSegment ? 5 : 2;
//         context.stroke();
//       }
//       prev = { u, x, y };
//     }
//   }

function drawBSpline() {
    uMin = U[order - 1].valueOf();
    uMax = U[U.length - order].valueOf();
    const step = (uMax - uMin) / 400;
    let first = true;
    
    for (let u = uMin; u <= uMax - step; u += step) {
        const [x, y] = evalCurve(u);
        const [x1, y1] = evalCurve(u + step);
        if (first) {
            context?.beginPath();
            context.moveTo(x, y);
            context.lineTo(x1, y1);
            context.strokeStyle = "blue";
            first = false;
        } else {

            if (highlightUEnabled && u > highlightU - step * 4 && u < highlightU + step * 4) {
                context.strokeStyle = "red";
                context.lineWidth = 5;
            } else {
                context.strokeStyle = "blue";
                context.lineWidth = 2;
            }
            context?.beginPath();
            context.lineTo(x1, y1);
            context.lineTo(x, y);
        }
        context.stroke();
    }
}
  
// function drawBSpline() {
//     uMin = U[order - 1].valueOf();
//     uMax = U[U.length - order].valueOf();
//     const step = (uMax - uMin) / 200;
//     let first = true;
//     context?.beginPath();
//     for (let u = uMin; u <= uMax; u += step) {
//         const [x, y] = evalCurve(u);
//         if (first) {
//             context.moveTo(x, y);
//             first = false;
//         } else {
//             context.lineTo(x, y);
//         }
//     }
//     context.strokeStyle = "blue";
//     context.lineWidth = 2;
//     context.stroke();
// }

/**
 * Build a clamped, uniform knot vector for a B-spline
 * with 'nPts' control points and spline order = 'd'.
 * Returns an array of length nPts + d.
 */
function makeOpenUniformKnots(nPts, d) {
    const m = nPts + d;
    const U = new Array(m);
    // clamped at start
    for (let i = 0; i < d; i++) {
        U[i] = 0;
    }
    // interior: 1, 2, ..., nPts - d
    for (let i = d; i < nPts; i++) {
        U[i] = i - d + 1;
    }
    // clamped at end
    for (let i = nPts; i < m; i++) {
        U[i] = nPts - d + 1;
    }
    return U;
}

function updateKnotsAndBasis() {
    // 1) rebuild the plain-number knot vector
    if (num_points !== thePoints.length) {
        knot_vectors = makeOpenUniformKnots(thePoints.length, order);
        num_points = thePoints.length;
        U = knot_vectors.map(toRational);
        Bs = calc_BlendingFunction(U, order);
        uMin = U[order - 1].valueOf();
        uMax = U[U.length - order].valueOf();
        slider2.min   = uMin;
        slider2.max   = uMax;
        slider2.value = uMin;
        slider2.step  = (uMax - uMin) / 200;
        highlightU = parseFloat(slider2.value);
        uDisplay.value = highlightU.toFixed(2);
    }

    // 2) convert to Rationals
    U = knot_vectors.map(toRational);

    // 3) re-compute the blending functions
    Bs = calc_BlendingFunction(U, order);

    // 4) re-render the LaTeX
    // piecewiseToLatex(Bs, order);
    piecewiseWithPlots(Bs, order);
    renderKnotUI();
}

/**
 * Based on the order d of the B-spline, create a default set of control points of minimum number d.
 * @param {number} n 
 * @returns {Array<number[]>} an array of control points
 */
function makeDefaultControlPoints(n) {
    let points = [];
    for (let i = 1; i <= n; i++) {
        let l1 = 50;
        let l2 = canvas.width - 100;
        let h = canvas.height / 2;
        points.push([l1 + l2 / (n - 1) * (i - 1), h + (-1)**i * 100]);
    }
    return points;
}

/**
 * Show the current knot vector as text,
 * and build one <input> per knot.
 */
function renderKnotUI() {
    const display = document.getElementById("knot-vector-display");
    const inputs  = document.getElementById("knot-inputs");
    // show it as a JS array literal
    display.textContent = "Knot vector U = [" + knot_vectors.join(", ") + "]";
    
    // clear out old inputs
    inputs.innerHTML = "Specify my own knot vector (non-decreasing):";
    // for each knot, make a number‐input
    knot_vectors.forEach((u,i) => {
      const inp = document.createElement("input");
      inp.type    = "number";
      inp.step    = "any";
      inp.value   = u;
      inp.size    = 4;
      inp.title   = `U[${i}]`;
      inp.dataset.idx = i;
      inp.style = "text-align: center;"
      // when user types a new value:
      inp.addEventListener("change", onUserKnotChange);
      inputs.appendChild(inp);
    });
}

/**
 * Called when the user edits one of the knot inputs.
 * We overwrite `knot_vectors[idx]`, rebuild U, recompute Bs, re‐render formulas+plots+curve.
 */
function onUserKnotChange(evt) {
    const idx = +evt.target.dataset.idx;
    const val = parseFloat(evt.target.value);
    // update, recompute, redraw
    knot_vectors[idx] = val;
    U  = knot_vectors.map(toRational);
    Bs = calc_BlendingFunction(U, order);
    uMin = U[order - 1].valueOf();
    uMax = U[U.length - order].valueOf();
    slider2.min   = uMin;
    slider2.max   = uMax;
    slider2.value = uMin;
    slider2.step  = (uMax - uMin) / 200;
    highlightU = parseFloat(slider2.value);
    uDisplay.value = highlightU.toFixed(2);
    piecewiseWithPlots(Bs, order);
    wrapDraw();
    // update the textual display, in case you want it live
    document.getElementById("knot-vector-display")
            .textContent = "Knot vector U = [" + knot_vectors.join(", ") + "]";
}


// helper function - set the slider to have max = # of control points
function wrapDraw() {
    // 1. Clear the canvas
    context.clearRect(0, 0, canvas.width, canvas.height);
    makeEffectiveFlag();
    // 2. Draw each control point
    let i = 0;
    for (let [x, y] of thePoints) {
        context.beginPath();
        context.arc(x, y, 5, 0, Math.PI * 2);  // radius = 5
        if (effectiveFlag[i] === 1 && highlightCPandBFEnabled) {
            context.fillStyle = "#FF8800";
        }
        else {
            context.fillStyle = "black";
        }
        context.fill();
        i++;
    }

    // 3. Draw the B-spline curve
    drawBSpline();
}

function setNumPoints() {
    runcanvas.setupSlider(1, 10, 1);
}

function makeEffectiveFlag() {
    for (let i = 0; i < thePoints.length; i++) {
        let value = evalBasis(i, highlightU);
        if (value > 0) {
            effectiveFlag[i] = 1;
        }
        else {
            effectiveFlag[i] = 0;
        }
    }
}




let canvas = /** @type {HTMLCanvasElement} */ (document.getElementById("canvas1"));
let context = canvas.getContext("2d");

/** @type Array<number[]> */
let thePoints = [
    [150, 150],
    [150, 450],
    [450, 450],
    [450, 150]
];
let effectiveFlag = []; // 1: effective, 0: ineffective
let highlightU = 0;           // current slider value
let highlightUEnabled = false; // checkbox state
let highlightCPandBFEnabled = false; // checkbox state, whether to highlight control points and blending functions
// grab the DOM nodes
const slider = document.getElementById("degree-slider");
const degreeDisplay = document.getElementById("degree-value");
const slider2 = document.getElementById("u-slider");
const uDisplay = document.getElementById("u-value");
const toggle1 = document.getElementById("highlight-toggle");
const toggle2 = document.getElementById("highlight-effective-toggle");

// whenever the slider moves…
slider2.addEventListener("input", () => {
  highlightU = parseFloat(slider2.value);
  uDisplay.value = highlightU.toFixed(2);
  wrapDraw();
  updateKnotsAndBasis();
});

// whenever the checkbox flips…
toggle1.addEventListener("change", () => {
  highlightUEnabled = toggle1.checked;
  wrapDraw();
  updateKnotsAndBasis();
});
toggle2.addEventListener("change", () => {
  highlightCPandBFEnabled = toggle2.checked;
  wrapDraw();
  updateKnotsAndBasis();
});
  
let num_points = thePoints.length;
// let knot_vectors = [0, 1, 2, 3, 4, 5];
let knot_vectors = [0, 0, 1, 2, 3, 3];
// let knot_vectors = [0, 1, 2, 3, 4, 5, 6, 7];
let U = knot_vectors.map(toRational);
let degree = 1;
let order = degree + 1;
let uMin = U[order - 1].valueOf();
let uMax = U[U.length - order].valueOf();
let Bs = calc_BlendingFunction(U, order);

let runcanvas = new RunCanvas(canvas, wrapDraw);

setNumPoints();
runcanvas.setValue(1);
runcanvas.digits = 0;



// let slider = runcanvas.range;
degree = Number(slider.value);
order = degree + 1;
thePoints = makeDefaultControlPoints(order)
num_points = thePoints.length;
knot_vectors = makeOpenUniformKnots(thePoints.length, order);
U = knot_vectors.map(toRational);
Bs = calc_BlendingFunction(U, order);

// 1. Determine the active u-range
uMin = U[order - 1].valueOf();
uMax = U[U.length - order].valueOf();
slider2.min   = uMin;
slider2.max   = uMax;
slider2.value = uMin;
slider2.step  = (uMax - uMin) / 200;
highlightU = parseFloat(slider2.value);
uDisplay.value = highlightU.toFixed(2);
wrapDraw();
renderKnotUI();

// add the point dragging UI
draggablePoints(canvas, thePoints, () => {
    // whenever the user adds or drags a point…
    updateKnotsAndBasis();
    wrapDraw();
    setNumPoints();
}, 10, setNumPoints);

slider.oninput = function() {
    degree = Number(slider.value);
    degreeDisplay.value = degree;
}

let prevDegree = degree;
slider.onchange = function() {
    // whenever the user moves the slider to change the degree of the B-spline…
    if (prevDegree !== Number(slider.value)) {
        prevDegree = Number(slider.value);
        degreeDisplay.value = degree;
        degree = Number(slider.value);
        order = degree + 1;
        thePoints = makeDefaultControlPoints(order);
        knot_vectors = makeOpenUniformKnots(thePoints.length, order);
        U = knot_vectors.map(toRational);
        Bs = calc_BlendingFunction(U, order);
        uMin = U[order - 1].valueOf();
        uMax = U[U.length - order].valueOf();
        slider2.min   = uMin;
        slider2.max   = uMax;
        slider2.value = uMin;
        slider2.step  = (uMax - uMin) / 200;
        highlightU = parseFloat(slider2.value);
        uDisplay.value = highlightU.toFixed(2);
        updateKnotsAndBasis();
        wrapDraw();
        renderKnotUI();
        draggablePoints(canvas, thePoints, () => {
            // whenever the user adds or drags a point…
            updateKnotsAndBasis();
            wrapDraw();
            renderKnotUI();
            setNumPoints();
        }, 10, setNumPoints);
    }
}

// let slider = runcanvas.range;
// let degree = Number(slider.value);
// let order = degree + 1;
// let knot_vectors = make(thePoints.length, order);
// console.log(knot_vectors);
// let U = knot_vectors.map(toRational);
// let Bs = calc_BlendingFunction(U, order);
// let uMin = U[order - 1].valueOf();
// let uMax = U[U.length - order].valueOf();

// let slider = runcanvas.range;
// const n = Number(slider.value);
// runcanvas.setValue(n);
// thePoints = makeDefaultControlPoints(n);
// updateKnotsAndBasis();
// wrapDraw();
// slider.onchange = function() {
//     // whenever the user moves the slider…
//     const n = Number(slider.value);
//     runcanvas.setValue(n);
//     thePoints = makeDefaultControlPoints(n);
//     updateKnotsAndBasis();
//     wrapDraw();
// };
  
// Bs[i] is now an array of segments; each segment has
//    - .poly: coefficients [a0, a1, ...] => a0 + a1 * u + a2 * u^2 + ...
//    - .interval: [u_lo, u_hi] on which that polynomial is valid

// Bs is Array of basis‐function segments; for clamped p=2 on knots=[0,1,2,3,4,5]
const container = document.getElementById("blend-formulas");
container.innerHTML = ""; 
// after you compute Bs = calc_BlendingFunction(…), pass them in:
// piecewiseToLatex(Bs, order);
piecewiseWithPlots(Bs, order);
