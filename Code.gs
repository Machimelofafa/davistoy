/* ====================================================================================================
   Code.gs  –  server-side maths & HTML bootstrapping
   ==================================================================================================== */

/**
 * Global constant used by the solver to guard against runaway values.
 */
const RTH_LIMIT = 100; // abort if cumulative Rth exceeds this

/**
 * Validate the incoming payload.
 * Throws descriptive errors when required fields are missing or invalid.
 */
function validatePayload(p) {
  var numericFields = ['srcLen', 'srcWid', 'dies', 'hConv', 'coolerRth', 'spacingX', 'spacingY', 'diePower'];
  numericFields.forEach(function(field) {
    if (typeof p[field] !== 'number' || isNaN(p[field])) {
      throw new Error('Invalid payload: "' + field + '" must be a number');
    }
  });
  if (p.srcLen <= 0 || p.srcWid <= 0 || p.dies <= 0) {
    throw new Error('Source dimensions and dies must be positive');
  }
  if (p.coolerMode === 'conv' && p.hConv <= 0) {
    throw new Error('hConv must be positive for convection mode');
  }
  if (p.coolerMode === 'direct' && p.coolerRth <= 0) {
    throw new Error('coolerRth must be positive for direct mode');
  }
  if (typeof p.diePower === 'number' && p.diePower <= 0) {
    throw new Error('diePower must be positive');
  }
  if (typeof p.spacingX === 'number' && p.spacingX < 0) {
    throw new Error('spacingX must be non-negative');
  }
  if (typeof p.spacingY === 'number' && p.spacingY < 0) {
    throw new Error('spacingY must be non-negative');
  }
  if (!Array.isArray(p.layers)) {
    throw new Error('Invalid payload: "layers" must be an array');
  }
  if (typeof p.layout !== 'string') p.layout = 'line';
  if (typeof p.coords !== 'string') p.coords = '';
  p.layers.forEach(function(L, idx) {
    if (typeof L !== 'object' || L === null) {
      throw new Error('Invalid layer at index ' + idx + ': must be an object');
    }
    ['t', 'kxy', 'kz'].forEach(function(prop) {
      if (typeof L[prop] !== 'number' || isNaN(L[prop])) {
        throw new Error('Invalid layer at index ' + idx + ': "' + prop + '" must be a number');
      }
      if (L[prop] <= 0) {
        throw new Error('Layer ' + (idx + 1) + ' property "' + prop + '" must be positive');
      }
    });
  });
}

/**
 * Calculate vertical resistance and final footprint for a single die.
 */
function solveSingleDieStack(p1) {
  let len = p1.srcLen / 1000;
  let wid = p1.srcWid / 1000;

  const rEach = [];
  const rCum  = [];
  const lengths = [len];
  const widths  = [wid];
  const widthsX = [wid];
  const widthsY = [len];

  p1.layers.forEach(L => {
    let current_layer_total_R = 0;
    const layer_thickness_in_microns = L.t;

    const ratioIso = L.kxy / L.kz;
    const tanIso   = Math.tan(Math.atan(Math.sqrt(ratioIso)));
    const tanXY    = Math.tan(Math.atan(Math.sqrt(L.kxy / L.kz)));

    let slice_iter_len_iso = len;
    let slice_iter_wid_iso = wid;
    let slice_iter_widX_aniso = widthsX[widthsX.length - 1];
    let slice_iter_widY_aniso = widthsY[widthsY.length - 1];

    if (layer_thickness_in_microns > 0) {
      const step_um = Math.max(1, layer_thickness_in_microns / 100);
      for (let processed = 0; processed < layer_thickness_in_microns; processed += step_um) {
        const slice_um = Math.min(step_um, layer_thickness_in_microns - processed);
        const t_micro_slice_m = slice_um * 1e-6;

        const A_slice = slice_iter_len_iso * slice_iter_wid_iso;
        let Ri_slice = Infinity;
        if (L.kz > 0 && A_slice > 0) {
          Ri_slice = t_micro_slice_m / (L.kz * A_slice);
        }
        current_layer_total_R += Ri_slice;
        const runningTotal = (rCum.length > 0 ? rCum[rCum.length - 1] : 0) + current_layer_total_R;
        if (runningTotal > RTH_LIMIT) {
          throw new Error('Cumulative Rth exceeds ' + RTH_LIMIT + ' \xB0C/W');
        }

        const delta_iso_slice = 2 * t_micro_slice_m * tanIso;
        slice_iter_len_iso += delta_iso_slice;
        slice_iter_wid_iso += delta_iso_slice;

        const delta_aniso = 2 * t_micro_slice_m * tanXY;
        slice_iter_widX_aniso += delta_aniso;
        slice_iter_widY_aniso += delta_aniso;
      }
    } else {
      current_layer_total_R = 0;
    }

    rEach.push(current_layer_total_R);
    const lastRCum = rCum.length > 0 ? rCum[rCum.length - 1] : 0;
    rCum.push(lastRCum + current_layer_total_R);

    len = slice_iter_len_iso;
    wid = slice_iter_wid_iso;

    lengths.push(len);
    widths.push(wid);
    widthsX.push(slice_iter_widX_aniso);
    widthsY.push(slice_iter_widY_aniso);
  });

  let rCool = 0;
  const final_area_for_cooler = len * wid;
  if (p1.coolerMode === 'conv') {
    if (p1.hConv > 0 && final_area_for_cooler > 0) {
      rCool = 1 / (p1.hConv * final_area_for_cooler);
    } else {
      rCool = Infinity;
    }
  } else if (p1.coolerMode === 'direct') {
    rCool = p1.coolerRth;
  }

  const rStack = rCum.length > 0 ? rCum[rCum.length - 1] : 0;
  const rVert = rStack + rCool;

  return { rEach, rCum, lengths, widths, widthsX, widthsY, rCool, rStack, rVert };
}

/**
 * Calculate die coordinates for a given layout.
 * @param {string} layout - The name of the layout ('line', '2-lines', '2-linesQ', 'custom').
 * @param {number} n - The number of dies.
 * @param {number} spacingX - The horizontal center-to-center spacing in mm.
 * @param {number} spacingY - The vertical center-to-center spacing in mm.
 * @param {string} custom - The custom coordinate string.
 * @returns {Array<Array<number>>} An array of [x,y] coordinates for each die.
 */
function getCoords(layout, n, spacingX, spacingY, custom) {
  const arr = [];
  const sX = typeof spacingX === 'number' ? spacingX : 0;
  const sY = typeof spacingY === 'number' ? spacingY : 0;
  const halfX = sX / 2;
  const halfY = sY / 2;

  switch (layout) {
    case '2-lines':
      for (let i = 0; i < n; i++) {
        const row = Math.floor(i / 2); // 0 for first column, 1 for second
        const col = i % 2;
        const x = row * sX;
        const y = col * sY - halfY;
        arr.push([x, y]);
      }
      break;

    case '2-linesQ': // Quincunx/staggered
      for (let i = 0; i < n; i++) {
        const row = Math.floor(i / 2);
        const col = i % 2;
        const x = row * sX + (col * halfX); // Stagger x
        const y = col * sY - halfY;
        arr.push([x, y]);
      }
      break;

    case 'custom':
      if (typeof custom === 'string' && custom.length > 0) {
        custom.split(';').forEach(p => {
          const parts = p.split(',');
          if (parts.length === 2) {
            arr.push([parseFloat(parts[0]) || 0, parseFloat(parts[1]) || 0]);
          }
        });
      }
      // Pad with default coords if not enough were provided
      while (arr.length < n) arr.push([arr.length * sX, 0]);
      break;

    case 'line':
    default:
      for (let i = 0; i < n; i++) {
        arr.push([i * sX, 0]);
      }
      break;
  }
  return arr.slice(0, n);
}

/** Utility to compute distance between two coords in metres. */
function distanceInMeters(a,b){
  const dx = (a[0]-b[0]) / 1000;
  const dy = (a[1]-b[1]) / 1000;
  return Math.sqrt(dx*dx + dy*dy);
}

/** Area overlap of two axis-aligned rectangles. */
function rectangleOverlap(r1, r2) {
  const x_overlap = Math.max(0, Math.min(r1.x1, r2.x1) - Math.max(r1.x0, r2.x0));
  const y_overlap = Math.max(0, Math.min(r1.y1, r2.y1) - Math.max(r1.y0, r2.y0));
  return x_overlap * y_overlap;
}

/** Compute the union area of a set of axis-aligned rectangles. */
function computeUnionArea(rects) {
  let events = [];
  for (let r of rects) {
    events.push({ x: r.x0, y0: r.y0, y1: r.y1, add: 1 });
    events.push({ x: r.x1, y0: r.y0, y1: r.y1, add: -1 });
  }
  events.sort((a, b) => a.x - b.x);
  let active = new Map();
  function modify(y0, y1, delta) {
    const key = y0 + ',' + y1;
    const prev = active.get(key) || 0;
    const next = prev + delta;
    if (next <= 0) active.delete(key);
    else active.set(key, next);
  }
  function totalYcovered() {
    if (active.size === 0) return 0;
    let segs = [];
    for (let key of active.keys()) {
      const [a, b] = key.split(',').map(Number);
      segs.push([a, b]);
    }
    segs.sort((s1, s2) => s1[0] - s2[0]);
    let length = 0;
    let [cur0, cur1] = segs[0];
    for (let i = 1; i < segs.length; i++) {
      const [n0, n1] = segs[i];
      if (n0 <= cur1) {
        cur1 = Math.max(cur1, n1);
      } else {
        length += (cur1 - cur0);
        [cur0, cur1] = [n0, n1];
      }
    }
    length += (cur1 - cur0);
    return length;
  }
  let prevX = events[0].x;
  let area = 0;
  for (let ev of events) {
    const curX = ev.x;
    const dx = curX - prevX;
    if (dx > 0) {
      const yCover = totalYcovered();
      area += dx * yCover;
      prevX = curX;
    }
    modify(ev.y0, ev.y1, ev.add);
  }
  return area;
}

/**
 * Build the final conductance matrix including lateral coupling and cooler.
 */
function buildConductanceMatrix(p, singleDieResults) {
  const N = singleDieResults.length;
  const coords = getCoords(p.layout, N, p.spacingX, p.spacingY, p.coords);

  const numLayers = p.layers.length;
  let G = Array.from({ length: N }, () => Array(N).fill(0));

  for (let i = 0; i < N; i++) {
    for (let j = i + 1; j < N; j++) {
      const d_ij = distanceInMeters(coords[i], coords[j]);
      let R_ij = 0;
      let anyOverlap = false;

      for (let l = 0; l < numLayers; l++) {
        const rect_i = {
          x0: coords[i][0]/1000 - singleDieResults[i].widthsX[l + 1]/2,
          x1: coords[i][0]/1000 + singleDieResults[i].widthsX[l + 1]/2,
          y0: coords[i][1]/1000 - singleDieResults[i].widthsY[l + 1]/2,
          y1: coords[i][1]/1000 + singleDieResults[i].widthsY[l + 1]/2
        };
        const rect_j = {
          x0: coords[j][0]/1000 - singleDieResults[j].widthsX[l + 1]/2,
          x1: coords[j][0]/1000 + singleDieResults[j].widthsX[l + 1]/2,
          y0: coords[j][1]/1000 - singleDieResults[j].widthsY[l + 1]/2,
          y1: coords[j][1]/1000 + singleDieResults[j].widthsY[l + 1]/2
        };
        const area = rectangleOverlap(rect_i, rect_j);
        if (area > 0) {
          anyOverlap = true;
          const R_l = (p.layers[l].t * 1e-6) / (p.layers[l].kxy * area);
          R_ij += R_l;
        }
      }

      if (anyOverlap) {
        const gij = 1 / R_ij;
        G[i][i] += gij;
        G[j][j] += gij;
        G[i][j] -= gij;
        G[j][i] -= gij;
      }
    }
  }

  const Gold = G;
  const Nplus1 = N + 1;
  let Gnew = new Array(Nplus1);
  for (let i = 0; i < Nplus1; i++) {
    Gnew[i] = new Array(Nplus1).fill(0);
  }
  for (let i = 0; i < N; i++) {
    for (let j = 0; j < N; j++) {
      Gnew[i][j] = Gold[i][j];
    }
  }

  for (let i = 0; i < N; i++) {
    const R_stack_i = singleDieResults[i].rStack;
    const Gstack_i  = R_stack_i > 0 ? 1 / R_stack_i : 0;
    Gnew[i][N]    -= Gstack_i;
    Gnew[N][i]    -= Gstack_i;
    Gnew[i][i]    += Gstack_i;
    Gnew[N][N]    += Gstack_i;
  }

  const bottomIndex = numLayers;
  let footprints = [];
  for (let i = 0; i < N; i++) {
    const cx = coords[i][0] / 1000.0;
    const cy = coords[i][1] / 1000.0;
    const wX = singleDieResults[i].widthsX[bottomIndex];
    const wY = singleDieResults[i].widthsY[bottomIndex];
    footprints.push({ x0: cx - wX/2, x1: cx + wX/2, y0: cy - wY/2, y1: cy + wY/2 });
  }

  const A_union = computeUnionArea(footprints);
  let G_cool_tot = 0;
  if (p.coolerMode === 'conv') {
    if (p.hConv > 0 && A_union > 0) {
      G_cool_tot = p.hConv * A_union;
    }
  } else if (p.coolerMode === 'direct') {
    if (p.coolerRth > 0) {
      G_cool_tot = 1 / p.coolerRth;
    }
  }
  Gnew[N][N] += G_cool_tot;

  return { G: Gnew, coords };
}

/**
 * Solve the matrix system and interpret the temperature results.
 */
function solveSystem(G, Pvec) {
  function solveMatrix(A,b){
    const n=A.length; const B=b.slice();
    const M=A.map(r=>r.slice());
    for(let i=0;i<n;i++){
      let max=i; for(let j=i+1;j<n;j++) if(Math.abs(M[j][i])>Math.abs(M[max][i])) max=j;
      const tmp=M[i]; M[i]=M[max]; M[max]=tmp; const t=B[i]; B[i]=B[max]; B[max]=t;
      const piv=M[i][i]; if(Math.abs(piv)<1e-12) return Array(n).fill(NaN);
      for(let j=i+1;j<n;j++){ const f=M[j][i]/piv; for(let k=i;k<n;k++) M[j][k]-=f*M[i][k]; B[j]-=f*B[i]; }
    }
    const x=Array(n).fill(0);
    for(let i=n-1;i>=0;i--){ let sum=B[i]; for(let j=i+1;j<n;j++) sum-=M[i][j]*x[j]; x[i]=sum/M[i][i]; }
    return x;
  }


  const temps = solveMatrix(G,Pvec);
  const N = Pvec.length - 1;
  const dieTemps = temps.slice(0, N);
  const maxTemp = Math.max.apply(null, dieTemps.filter(v=>typeof v==='number'));
  const avgTemp = dieTemps.reduce((a,b)=>a+b,0)/dieTemps.length;

  return { dieTemps, maxTemp, avgTemp };
}

/**
 * Solve the layer stack:
 * • conduction per layer (Ri) and cumulative, with 1µm discretization
 * • isotropic cone footprint  (lengths[], widths[])
 * • orthotropic footprints    (widthsX[], widthsY[])
 * • cooler term (direct or convection)
 * • per-die Rth and total-stack Rth
 *
 * @param {Object} p  payload from the browser (see ui.html for structure)
 * @return {Object}   results used by ui.html
 */
function coreSolve(p) {
  validatePayload(p);

  const N = Math.max(1, Math.floor(p.dies));
  const results = [];
  let tplResult = null;
  for (let i = 0; i < N; i++) {
    const p_i = Object.assign({}, p, { dies: 1 });
    const res_i = solveSingleDieStack(p_i);
    if (i === 0) tplResult = res_i;
    results[i] = res_i;
  }

  const matrixInfo = buildConductanceMatrix(p, results);
  const NN = matrixInfo.G.length;
  let Pvec = new Array(NN).fill(0);
  for (let i = 0; i < N; i++) Pvec[i] = p.diePower || 1.0;
  Pvec[N] = 0;

  const sys = solveSystem(matrixInfo.G, Pvec);

  // ----- Compute Crosstalk Matrix -----
  const rthMatrix = Array.from({ length: N }, () => Array(N).fill(0));
  for (let j = 0; j < N; j++) {
    const Ptest = new Array(NN).fill(0);
    Ptest[j] = 1; // unit power for source die j
    const temps = solveSystem(matrixInfo.G, Ptest).dieTemps;
    for (let i = 0; i < N; i++) {
      rthMatrix[i][j] = temps[i];
    }
  }
  const tempCrosstalkMatrix = rthMatrix.map(row => row.map(v => v * (p.diePower || 1)));

  const rDie = sys.maxTemp;
  const P_total = N * (p.diePower || 1.0);
  const rTotal = P_total > 0 ? sys.maxTemp / P_total : 0;

  const { rEach, rCum, lengths, widths, widthsX, widthsY, rCool, rStack } = tplResult;
  const materialNames = p.layers.map(function(L){ return typeof L.mat === 'string' ? L.mat : ''; });

  return {
    rEach,
    rCum,
    widths,
    lengths,
    widthsX,
    widthsY,
    materialNames,
    coords: matrixInfo.coords,
    rDie,
    rTotal,
    numDies: N,
    rCoolPerDie: rCool,
    rStack,
    rDieList: sys.dieTemps,
    rDieAvg: sys.avgTemp,
    rthMatrix,
    tempCrosstalkMatrix
  };
}
function computeSensitivity(p, baseR) {
  var frac = 0.01;
  var layerSens = p.layers.map(function(L, idx) {
    var res = {};
    ['t','kxy','kz'].forEach(function(par){
      var delta = (L[par] || 0) * frac || frac;
      var plus = JSON.parse(JSON.stringify(p));
      plus.sensitivity = false;
      plus.layers[idx][par] = L[par] + delta;
      var minus = JSON.parse(JSON.stringify(p));
      minus.sensitivity = false;
      minus.layers[idx][par] = Math.max(L[par] - delta, 1e-9);
      var rPlus = coreSolve(plus).rDie;
      var rMinus = coreSolve(minus).rDie;
      res[par] = ((rPlus - rMinus) / (2 * delta)) * (L[par] / baseR);
    });
    return res;
  });

  var cooler = null;
  if (p.coolerMode === 'conv') {
    var deltaC = (p.hConv || 0) * frac || frac;
    var plusC = Object.assign({}, p, { hConv: p.hConv + deltaC, sensitivity:false });
    var minusC = Object.assign({}, p, { hConv: Math.max(p.hConv - deltaC, 1e-9), sensitivity:false });
    var rPlusC = coreSolve(plusC).rDie;
    var rMinusC = coreSolve(minusC).rDie;
    cooler = { hConv: ((rPlusC - rMinusC) / (2 * deltaC)) * (p.hConv / baseR) };
  } else if (p.coolerMode === 'direct') {
    var deltaD = (p.coolerRth || 0) * frac || frac;
    var plusD = Object.assign({}, p, { coolerRth: p.coolerRth + deltaD, sensitivity:false });
    var minusD = Object.assign({}, p, { coolerRth: Math.max(p.coolerRth - deltaD, 1e-9), sensitivity:false });
    var rPlusD = coreSolve(plusD).rDie;
    var rMinusD = coreSolve(minusD).rDie;
    cooler = { coolerRth: ((rPlusD - rMinusD) / (2 * deltaD)) * (p.coolerRth / baseR) };
  }

  return { layers: layerSens, cooler: cooler };
}

function solve(p) {
  var base = coreSolve(p);
  if (p && p.sensitivity) {
    base.sensitivity = computeSensitivity(p, base.rDie);
  }
  return base;
}

/* ==================== Monte Carlo solver ==================== */
function randomNormal() {
  let u = 0, v = 0;
  while (u === 0) u = Math.random();
  while (v === 0) v = Math.random();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}

function solveMonteCarlo(p) {
  const iter = Math.max(1, Math.floor(p.iterations || 1));
  const uncT = typeof p.uncT === 'number' ? p.uncT : 0;
  const uncK = typeof p.uncK === 'number' ? p.uncK : 0;

  const results = [];
  const critical = Array(p.layers.length).fill(0);

  for (let i = 0; i < iter; i++) {
    const perturbed = p.layers.map(L => ({
      t:   L.t   * (1 + randomNormal() * uncT),
      kxy: L.kxy * (1 + randomNormal() * uncK),
      kz:  L.kz  * (1 + randomNormal() * uncK)
    }));

    try {
      const r = solve({
        srcLen: p.srcLen,
        srcWid: p.srcWid,
        dies:   p.dies,
        spacingX: p.spacingX,
        spacingY: p.spacingY,
        layout:  p.layout,
        coords:  p.coords,
        coolerMode: p.coolerMode,
        coolerRth:  p.coolerRth,
        hConv:      p.hConv,
        diePower:  p.diePower,
        layers: perturbed
      });
      results.push(r.rDie);
      let maxV = -Infinity, maxI = 0;
      r.rEach.forEach((val, idx) => {
        if (typeof val === 'number' && val > maxV) { maxV = val; maxI = idx; }
      });
      critical[maxI]++;
    } catch (err) {
      // skip failed iteration
    }
  }

  const n = results.length;
  if (n === 0) {
    return { iterations: 0, results: [], critical };
  }

  const mean = results.reduce((a, b) => a + b, 0) / n;
  const variance = results.reduce((a, b) => a + Math.pow(b - mean, 2), 0) / n;
  const sorted = results.slice().sort((a, b) => a - b);
  const median = n % 2 ? sorted[(n - 1) / 2] : (sorted[n / 2 - 1] + sorted[n / 2]) / 2;

  return {
    iterations: n,
    results,
    mean,
    stdev: Math.sqrt(variance),
    min: sorted[0],
    max: sorted[n - 1],
    median,
    critical
  };
}

function solveSweep(p) {
  validatePayload(p);
  const layerIdx = Math.max(0, Math.floor(p.targetLayer || 0));
  const param = String(p.targetParam || 't');
  const start = Number(p.startVal);
  const end = Number(p.endVal);
  const steps = Math.max(2, Math.floor(p.steps || 2));
  const results = [];
  for (let i = 0; i < steps; i++) {
    const val = start + (end - start) * i / (steps - 1);
    const layersCopy = p.layers.map(l => ({ mat:l.mat, t:l.t, kxy:l.kxy, kz:l.kz }));
    if (layersCopy[layerIdx]) layersCopy[layerIdx][param] = val;
    const r = coreSolve({
      srcLen: p.srcLen,
      srcWid: p.srcWid,
      dies:   p.dies,
      spacingX: p.spacingX,
      spacingY: p.spacingY,
      layout:  p.layout,
      coords:  p.coords,
      coolerMode: p.coolerMode,
      coolerRth:  p.coolerRth,
      hConv:      p.hConv,
      diePower:  p.diePower,
      layers: layersCopy
    });
    results.push({ sweptValue: val, resultValue: r.rTotal });
  }
  return results;
}

/* ====================================================================================================
   UI bootstrap helpers (doGet, inc) - These functions remain unchanged.
   ==================================================================================================== */
function doGet() { //
  return HtmlService
           .createTemplateFromFile('index') //
           .evaluate() //
           .setTitle('Stack-up Rth Tool'); //
}

function inc(name) { //
  return HtmlService
           .createHtmlOutputFromFile(name) //
           .getContent(); //
}
