/* ====================================================================================================
   Code.gs  –  server-side maths & HTML bootstrapping
   ==================================================================================================== */

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
function solve(p) {

  const RTH_LIMIT = 100; // abort if cumulative Rth exceeds this

  // Validate payload structure
  var numericFields = ['srcLen', 'srcWid', 'dies', 'hConv', 'coolerRth'];
  numericFields.forEach(function(field) {
    if (typeof p[field] !== 'number' || isNaN(p[field])) {
      throw new Error('Invalid payload: "' + field + '" must be a number');
    }
  });

  // Guard against non-positive values
  if (p.srcLen <= 0 || p.srcWid <= 0 || p.dies <= 0) {
    throw new Error('Source dimensions and dies must be positive');
  }
  if (p.coolerMode === 'conv' && p.hConv <= 0) {
    throw new Error('hConv must be positive for convection mode');
  }
  if (p.coolerMode === 'direct' && p.coolerRth <= 0) {
    throw new Error('coolerRth must be positive for direct mode');
  }

  if (!Array.isArray(p.layers)) {
    throw new Error('Invalid payload: "layers" must be an array');
  }
  p.layers.forEach(function(L, idx) {
    if (typeof L !== 'object' || L === null) {
      throw new Error('Invalid layer at index ' + idx + ': must be an object');
    }
    ['t', 'kx', 'ky', 'kz'].forEach(function(prop) {
      if (typeof L[prop] !== 'number' || isNaN(L[prop])) {
        throw new Error('Invalid layer at index ' + idx + ': "' + prop + '" must be a number');
      }
      if (L[prop] <= 0) {
        throw new Error('Layer ' + (idx + 1) + ' property "' + prop + '" must be positive');
      }
    });
  });

  /* ---------- initial footprint (m) ----------------------------- */
  let len = p.srcLen / 1000;   // depth  (Y)  — mm → m
  let wid = p.srcWid / 1000;   // width  (X)

  /* ---------- accumulators -------------------------------------- */
  const rEach   = [];          // individual Ri for each original layer
  const rCum    = [];          // cumulative Ri after each original layer
  
  const lengths = [len];       // Isotropic length (Y) after each layer (starts with initial)
  const widths  = [wid];       // Isotropic width (X) after each layer (starts with initial)

  const widthsX = [wid];       // Anisotropic width (X) after each layer (starts with initial)
  const widthsY = [len];       // Anisotropic length (Y) after each layer (starts with initial)

  /* ---------- layer loop ---------------------------------------- */
  p.layers.forEach(L => { //
    let current_layer_total_R = 0; //
    const layer_thickness_in_microns = L.t; //

    // Pre-compute α angles and their tangents for this layer
    const ratioIso = (L.kx + L.ky) / (2 * L.kz);
    const alphaIso = Math.atan(Math.sqrt(ratioIso));
    const tanIso   = Math.tan(alphaIso);

    const alphaX = Math.atan(Math.sqrt(L.kx / L.kz));
    const tanX   = Math.tan(alphaX);

    const alphaY = Math.atan(Math.sqrt(L.ky / L.kz));
    const tanY   = Math.tan(alphaY);

    let slice_iter_len_iso = len; //
    let slice_iter_wid_iso = wid; //
    
    let slice_iter_widX_aniso = widthsX[widthsX.length - 1]; //
    let slice_iter_widY_aniso = widthsY[widthsY.length - 1]; //

    if (layer_thickness_in_microns > 0) { //
      const step_um = Math.max(1, layer_thickness_in_microns / 100); // adaptive step
      for (let processed = 0; processed < layer_thickness_in_microns; processed += step_um) { //
        const slice_um = Math.min(step_um, layer_thickness_in_microns - processed); //
        const t_micro_slice_m = slice_um * 1e-6; //

        const A_slice = slice_iter_len_iso * slice_iter_wid_iso; //
        let Ri_slice = Infinity; //

        if (L.kz > 0 && A_slice > 0) { //
          Ri_slice = t_micro_slice_m / (L.kz * A_slice); //
        }
        current_layer_total_R += Ri_slice; //
        const runningTotal = (rCum.length > 0 ? rCum[rCum.length - 1] : 0) + current_layer_total_R;
        if (runningTotal > RTH_LIMIT) {
          throw new Error('Cumulative Rth exceeds ' + RTH_LIMIT + ' \xB0C/W');
        }

        const delta_iso_slice = 2 * t_micro_slice_m * tanIso;
        slice_iter_len_iso += delta_iso_slice; //
        slice_iter_wid_iso += delta_iso_slice; //

        const delta_aniso_X_slice = 2 * t_micro_slice_m * tanX;
        const delta_aniso_Y_slice = 2 * t_micro_slice_m * tanY;
        slice_iter_widX_aniso += delta_aniso_X_slice; //
        slice_iter_widY_aniso += delta_aniso_Y_slice; //
      }
    } else { 
        current_layer_total_R = 0; //
    }

    rEach.push(current_layer_total_R); //
    const lastRCum = rCum.length > 0 ? rCum[rCum.length - 1] : 0; //
    rCum.push(lastRCum + current_layer_total_R); //

    len = slice_iter_len_iso; //
    wid = slice_iter_wid_iso; //

    lengths.push(len); //
    widths.push(wid);  //
    widthsX.push(slice_iter_widX_aniso); //
    widthsY.push(slice_iter_widY_aniso); //
  });

  /* ---------- cooler term --------------------------------------- */
  let rCool = 0; //
  const final_area_for_cooler = len * wid; //
  if (p.coolerMode === 'conv') { //
    if (p.hConv > 0 && final_area_for_cooler > 0) { //
        rCool = 1 / (p.hConv * final_area_for_cooler); //
    } else {
        rCool = Infinity; //
    }
  } else if (p.coolerMode === 'direct') { //
    rCool = p.coolerRth; //
  }
  // If p.coolerMode is 'none', rCool remains 0 as initialized.

  /* ---------- per-die & total ----------------------------------- */
  const rStack = rCum.length > 0 ? rCum[rCum.length - 1] : 0; //
  const rDie = rStack + rCool; //
  
  let rTotal = Infinity; //
  if (p.dies > 0) { //
      rTotal = rDie / p.dies; //
  }

  /* ---------- back to the browser ------------------------------- */
  return {
    rEach,      // Array of R_th for each layer (per-die stack component)
    rCum,       // Array of cumulative R_th after each layer (per-die stack component)
    widths,     // Array of isotropic widths (X) after each layer (including initial)
    lengths,    // Array of isotropic lengths (Y) after each layer (including initial)
    widthsX,    // Array of anisotropic widths (X) after each layer (including initial)
    widthsY,    // Array of anisotropic lengths (Y) after each layer (including initial)
    rDie,       // Total R_th per die (stack + cooler)
    rTotal,     // Overall R_th for all dies combined
    numDies: p.dies, // Pass number of dies to client
    rCoolPerDie: rCool, // ADDED: Pass per-die cooler resistance to client
    rStack      // Send cumulative stack Rth back for percentage calculations
  };
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
      t:  L.t  * (1 + randomNormal() * uncT),
      kx: L.kx * (1 + randomNormal() * uncK),
      ky: L.ky * (1 + randomNormal() * uncK),
      kz: L.kz * (1 + randomNormal() * uncK)
    }));

    try {
      const r = solve({
        srcLen: p.srcLen,
        srcWid: p.srcWid,
        dies:   p.dies,
        coolerMode: p.coolerMode,
        coolerRth:  p.coolerRth,
        hConv:      p.hConv,
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

/* ==================== Thickness optimization ==================== */
function optimizeThicknesses(p) {
  if (!Array.isArray(p.layers) || p.layers.length === 0) {
    throw new Error('No layers provided');
  }
  let base;
  try {
    base = solve(p);
  } catch (err) {
    throw new Error('Invalid starting configuration: ' + err.message);
  }

  let current = p.layers.map(l => Object.assign({}, l));
  let bestR = base.rDie;
  let bestLayers = current.map(l => Object.assign({}, l));
  const maxIter = 20;
  const lr = 0.2;
  const eps = 1e-4;

  for (let iter = 0; iter < maxIter; iter++) {
    const grads = [];
    for (let i = 0; i < current.length; i++) {
      const dt = Math.max(1, current[i].t * 0.05);
      const orig = current[i].t;
      current[i].t = orig + dt;
      let plus;
      try { plus = solve(Object.assign({}, p, { layers: current })); }
      catch (e) { plus = { rDie: Infinity }; }
      current[i].t = Math.max(0.01, orig - dt);
      let minus;
      try { minus = solve(Object.assign({}, p, { layers: current })); }
      catch (e) { minus = { rDie: Infinity }; }
      current[i].t = orig;
      grads[i] = (plus.rDie - minus.rDie) / (2 * dt);
    }

    const proposed = current.map((l, i) => ({
      kx: l.kx, ky: l.ky, kz: l.kz,
      mat: l.mat,
      t: Math.max(0.1, l.t - lr * grads[i])
    }));

    let res;
    try { res = solve(Object.assign({}, p, { layers: proposed })); }
    catch (err) { break; }

    if (res.rDie < bestR - eps) {
      bestR = res.rDie;
      bestLayers = proposed.map(o => Object.assign({}, o));
      current = proposed;
    } else {
      break;
    }
  }

  const recs = bestLayers.map((l, i) => ({
    layer: i + 1,
    tOriginal: p.layers[i].t,
    tOptimized: l.t
  }));

  return {
    rDieStart: base.rDie,
    rDieOptimized: bestR,
    recommendations: recs
  };
}

/* ==================== Material substitution ==================== */
function analyzeMaterialSubstitution(p) {
  const db = {
    'Copper':   {kx:393, ky:393, kz:393, costFactor:3.5,  maxTemp:200},
    'Aluminum': {kx:205, ky:205, kz:205, costFactor:1.0,  maxTemp:150},
    'Silver':   {kx:419, ky:419, kz:419, costFactor:15.0, maxTemp:250},
    'Diamond':  {kx:2000,ky:2000,kz:2000,costFactor:100.0,maxTemp:400},
    'Silicon':  {kx:148, ky:148, kz:148, costFactor:2.0,  maxTemp:300},
    'Graphite': {kx:150, ky:150, kz:1000,costFactor:5.0,  maxTemp:350},
    'TIM Standard': {kx:5, ky:5, kz:5, costFactor:8.0, maxTemp:200},
    'TIM Premium':  {kx:15,ky:15,kz:15,costFactor:25.0,maxTemp:250}
  };

  if (!Array.isArray(p.layers) || p.layers.length === 0) {
    throw new Error('No layers provided');
  }

  let base;
  try { base = solve(p); } catch (err) { throw new Error('Invalid configuration'); }

  const suggestions = [];

  p.layers.forEach((layer, idx) => {
    for (const [name, mat] of Object.entries(db)) {
      const modified = p.layers.map((l, i) => i === idx ? Object.assign({}, l, { kx: mat.kx, ky: mat.ky, kz: mat.kz }) : Object.assign({}, l));
      let res;
      try { res = solve(Object.assign({}, p, { layers: modified })); }
      catch (err) { continue; }
      if (!isFinite(res.rDie)) continue;
      const improvement = base.rDie - res.rDie;
      if (improvement <= 0) continue;
      suggestions.push({
        layer: idx + 1,
        material: name,
        newRth: res.rDie,
        improvement,
        costFactor: mat.costFactor,
        ratio: improvement / mat.costFactor
      });
    }
  });

  suggestions.sort((a,b)=> b.ratio - a.ratio);

  return { baseRth: base.rDie, suggestions };
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
