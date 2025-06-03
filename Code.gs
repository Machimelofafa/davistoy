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
    rCoolPerDie: rCool // ADDED: Pass per-die cooler resistance to client
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
