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
function coreSolve(p) {

  const RTH_LIMIT = 100; // abort if cumulative Rth exceeds this

  // Validate payload structure
  var numericFields = ['srcLen', 'srcWid', 'dies', 'hConv', 'coolerRth', 'spacing'];
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
  if (typeof p.spacing === 'number' && p.spacing < 0) {
    throw new Error('spacing must be non-negative');
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
    const ratioIso = L.kxy / L.kz;
    const alphaIso = Math.atan(Math.sqrt(ratioIso));
    const tanIso   = Math.tan(alphaIso);

    const alphaXY = Math.atan(Math.sqrt(L.kxy / L.kz));
    const tanXY   = Math.tan(alphaXY);

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

        const delta_aniso_X_slice = 2 * t_micro_slice_m * tanXY;
        const delta_aniso_Y_slice = 2 * t_micro_slice_m * tanXY;
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
  const rVert = rStack + rCool; // vertical path for a single die

  let rTotal = Infinity; // legacy total without coupling
  if (p.dies > 0) {
      rTotal = rVert / p.dies;
  }

  /* ---------- thermal coupling network ------------------------- */
  const N = Math.max(1, Math.floor(p.dies));
  function getCoords(layout, n, spacing, custom) {
    const arr = [];
    const s = typeof spacing === 'number' ? spacing : 0;
    if (layout === 'square' && n >= 4) {
      const d = s / 2;
      arr.push([-d,-d],[d,-d],[-d,d],[d,d]);
      for (let i=4;i<n;i++) arr.push([0,0]);
    } else if (layout === 'diamond' && n >= 4) {
      const d = s / 2;
      arr.push([0,-d],[-d,0],[0,d],[d,0]);
      for (let i=4;i<n;i++) arr.push([0,0]);
    } else if (layout === 'custom' && typeof custom === 'string') {
      custom.split(';').forEach(p=>{
        const parts = p.split(',');
        if (parts.length===2) arr.push([parseFloat(parts[0])||0, parseFloat(parts[1])||0]);
      });
      while (arr.length < n) arr.push([0,0]);
    } else { // line
      for (let i=0;i<n;i++) arr.push([i*s,0]);
    }
    return arr.slice(0,n);
  }

  const coords = getCoords(p.layout, N, p.spacing, p.coords);
  const Ga = rVert > 0 ? 1/rVert : 0;
  const totalThick = p.layers.reduce((a,L)=>a + L.t*1e-6,0);
  const kCouple = p.layers[0] ? p.layers[0].kxy : 1;
  const crossArea = totalThick * wid; // rough approximation
  const G = Array.from({length:N},()=>Array(N).fill(0));
  for (let i=0;i<N;i++) {
    for (let j=i+1;j<N;j++) {
      const dx = (coords[i][0]-coords[j][0]) / 1000;
      const dy = (coords[i][1]-coords[j][1]) / 1000;
      const dist = Math.sqrt(dx*dx + dy*dy);
      if (dist>0 && crossArea>0 && kCouple>0) {
        const Rl = dist/(kCouple*crossArea);
        const gij = 1/Rl;
        G[i][i] += gij; G[j][j] += gij; G[i][j] -= gij; G[j][i] -= gij;
      }
    }
    G[i][i] += Ga;
  }

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

  const Pvec = Array(N).fill(1); // 1W each
  const temps = solveMatrix(G,Pvec);
  const maxTemp = Math.max.apply(null, temps.filter(v=>typeof v==='number'));
  const avgTemp = temps.reduce((a,b)=>a+b,0)/temps.length;
  const rDie = maxTemp;
  rTotal = maxTemp/ N;

  // build full resistance matrix via unit heat inputs
  function buildResistanceMatrix(M){
    const n=M.length;
    const res=[];
    for(let i=0;i<n;i++){
      const b=Array(n).fill(0); b[i]=1;
      res.push(solveMatrix(M,b));
    }
    return res;
  }

  const rMatrix = buildResistanceMatrix(G);
  const rPerDie = rMatrix.map((row,i)=> row[i]);

  /* ---------- back to the browser ------------------------------- */
  return {
    rEach,      // Array of R_th for each layer (per-die stack component)
    rCum,       // Array of cumulative R_th after each layer (per-die stack component)
    widths,     // Array of isotropic widths (X) after each layer (including initial)
    lengths,    // Array of isotropic lengths (Y) after each layer (including initial)
    widthsX,    // Array of anisotropic widths (X) after each layer (including initial)
    widthsY,    // Array of anisotropic lengths (Y) after each layer (including initial)
    coords,     // Die coordinates used for visualisation
    rDie,       // Worst-case R_th per die including coupling
    rTotal,     // Overall R_th for all dies combined
    numDies: p.dies, // Pass number of dies to client
    rCoolPerDie: rCool, // ADDED: Pass per-die cooler resistance to client
    rStack,     // Send cumulative stack Rth back for percentage calculations
    rDieList: temps,     // Temperatures for simultaneous 1W per die
    rDieAvg: avgTemp,
    rMatrix,   // Full resistance matrix
    rPerDie    // Individual die-to-ambient resistances
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
        spacing: p.spacing,
        layout:  p.layout,
        coords:  p.coords,
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
