# Stack-up Rth & DeltaT Calculator

This repository hosts a Google Apps Script web application designed to estimate the thermal resistance (Rth) and temperature rise (ΔT) of multi-layer material stack-ups, commonly found in electronic assemblies. It provides a comprehensive tool for thermal analysis, accounting for complex heat spreading phenomena and multi-die configurations.

## Key Features

*   **Anisotropic Heat Spreading:** Accurately models heat flow considering different thermal conductivities in the XY-plane and Z-direction.
*   **Advanced Multi-Die Support:** Analyze thermal interactions for various die layouts, including 'Line', '2-Lines', 'Quincunx', and 'Custom' arrangements, to determine maximum and average die temperatures.
*   **Temperature Rise Calculation:** Directly computes temperature increase (ΔT in °C) based on per-die power dissipation.
*   **Detailed Visualizations:** Offers insightful graphical representations, including X and Y-axis side views of heat dissipation and a top-down view of die layouts with final heat spread areas.
*   **Cooling Model:** Incorporates a cooling model using a convection coefficient based on the final spread area.
*   **Sensitivity Analysis:** Identifies critical layers and parameters that significantly impact thermal performance.
*   **Monte Carlo Simulation:** Quantifies uncertainty by assessing the impact of manufacturing tolerances and material variations on Rth.
*   **Parameter Sweep:** Allows sweeping a chosen layer parameter over a range to observe its effect on total Rth.
*   **Save/Load Configurations:** Enables saving and loading of stack-up configurations directly within the browser for easy design iteration and comparison.
*   **Progress and Error Reporting:** Provides clear feedback on calculation progress and highlights invalid inputs.

## Getting Started

This web application is built using Google Apps Script. To use it, you need to import the code into your Google Drive and deploy it as a web app.

### Import into Google Apps Script

1.  Open [script.new](https://script.new) to create a new, blank Google Apps Script project.
2.  Replace the default `Code.gs` content with the content from this repository's `Code.gs` file.
3.  Create new HTML files in your Apps Script project named `index.html`, `controls.html`, `ui.html`, and `styles.html`. Copy the respective contents from this repository into these new files.
4.  Open the project manifest file (`appsscript.json`) and replace its contents with the version provided in this repository.

### Deploy as a Web App

1.  In the Google Apps Script editor, navigate to **Deploy > New deployment**.
2.  Select **Web app** as the deployment type.
3.  Provide a descriptive name for your deployment and configure the execution and access settings as needed (e.g., 


Execute as: `Me` (your Google account) and Who has access: `Anyone` (or `Anyone with Google account` if you want to restrict access).
4.  Click **Deploy** and authorize the script if prompted.

After successful deployment, you will receive a web app URL. Visiting this URL in your browser will load the `index.html` and launch the Rth Calculator.

## Basic Usage

1.  **Define Heat Source:** Input the heat source dimensions (length and width), number of dies, power per die, and select a die layout (Line, 2-Lines, Quincunx, or Custom). For custom layouts, provide semicolon-separated coordinates.
2.  **Build Stack-Up:** Add material layers sequentially from the heat source to the cooler. For each layer, specify its material name, thickness (µm), in-plane thermal conductivity (kxy in W/(m·K)), and through-plane thermal conductivity (kz in W/(m·K)).
3.  **Configure Cooling:** Choose your cooling method: `None`, `Direct Rth` (specify a fixed thermal resistance), or `Convection` (specify a convection coefficient).
4.  **Calculate:** Click the "Calculate" button to perform the thermal analysis.
5.  **Analyze Results:** Review the primary outputs (ΔT per die, total stack Rth), detailed visualizations (X/Y axis thermal dissipation, top view layout), and the summary table which breaks down Rth contributions, cumulative Rth, and sensitivity indices for each layer.
6.  **Explore Advanced Features:** Use the "Monte Carlo" button to run uncertainty quantification or "Run Sweep Analysis" to plot the effect of varying a single parameter.
7.  **Save/Load:** Utilize the "💾 Save Stack" and "📂 Load Stack" buttons to save and retrieve your current stack-up configurations.

## Calculation Approach

The solver employs a detailed layer-by-layer approach:

*   **Slice-based Calculation:** Each material layer is divided into 1 µm slices. The vertical resistance is calculated for each slice, and the heat-spreading footprint is updated based on angles derived from the conductivity ratios.
*   **Anisotropic Spreading:** The model accounts for anisotropic heat spreading by tracking widths in the X and Y directions separately.
*   **Cumulative Resistance:** Resistance is accumulated slice by slice. The calculation aborts if the cumulative Rth exceeds `RTH_LIMIT` (100 °C/W) to prevent runaway values.
*   **Cooler Model Integration:** After processing all layers, the selected cooler model (direct resistance or convection based on the final spread area) is added to determine the per-die and total stack resistance.
*   **Lateral Coupling:** For multi-die configurations, the tool computes overlap areas between die footprints across layers to build a conductance matrix, accounting for lateral thermal coupling.
*   **System Solution:** A matrix system is solved to obtain die temperatures, from which the maximum and average temperature rises, and overall stack resistance, are derived.

## Assumptions & Limitations

*   **Units:**
    *   Heat source length and width are entered in **millimeters (mm)**.
    *   Layer thicknesses are in **micrometers (µm)**.
    *   Thermal conductivities (kxy, kz) are in **W/(m·K)**.
    *   Results (Rth) are reported in **°C/W**.
*   **Homogeneity:** Each layer is considered homogeneous.
*   **Idealized Spreading:** Heat spreading is idealized as described in the calculation approach.
*   **Rth Limit:** Calculations stop if cumulative Rth exceeds 100 °C/W.

## Contributing

Contributions are welcome! If you have improvements or bug fixes, please fork this repository and submit a pull request. You can manually test your updates using the "Run" button in the Google Apps Script editor. If you add unit tests, please ensure they pass before submitting your PR.
