# Stack-up Rth Calculator
 
This repository contains a Google Apps Script web app that estimates the thermal resistance (Rth) of a material stack-up. The UI is delivered from `index.html` and the logic resides in `Code.gs`.

## Key Features

- **Anisotropic Heat Spreading** between XY and Z directions
- **Advanced Multi-Die Support** for common and custom layouts
- **Temperature Rise Calculation** directly from per‑die power
- **Detailed Visualizations** including side views and top layouts
- **Cooling Model** using a convection coefficient and spread area
- **Sensitivity Analysis** highlighting critical layers
- **Monte Carlo Simulation** for uncertainty quantification
- **Parameter Sweep** plots the effect of varying layer properties
- **Improved Cooler Model** using the final spread area
- **Save/Load Configurations** and progress/error reporting
- **Light/Dark Theme** for comfortable viewing
 
 ## Import into Google Apps Script
 
 1. Open [script.new](https://script.new) to create a blank project.
 2. Replace the default `Code.gs` with the contents of this repo's `Code.gs`.
 3. Add new files named `index.html`, `controls.html`, `ui.html` and `styles.html`; paste in the matching contents from this repo.
 4. Open the project manifest (`appsscript.json`) and replace its contents with the version provided here.
 
 ## Deploy as a Web App
 
 1. In the editor choose **Deploy > New deployment**.
 2. Select **Web app** as the deployment type.
 3. Provide a description and confirm the execution and access settings.
 4. Click **Deploy** and authorize the script.
 
 After deployment you will receive a web app URL. Visiting that URL loads `index.html` and starts the tool.
 
 ## Basic Usage
 
1. Enter the heat source length and width.
2. Add material layers in the table.
   Each layer requires thickness **t**, in‑plane conductivity **kxy** and through‑plane conductivity **kz**.
3. Click **Run** to calculate thermal resistance.
4. Use **Monte Carlo** or **Parameter Sweep** to explore uncertainties and trends.
5. Hover over form labels to view tooltips explaining the parameters.
6. Tooltips on the results table and graphs explain the values shown.
 
## Calculation Approach

The solver divides each layer into 1&nbsp;µm slices and updates the heat-spreading
footprint after every slice using angles derived from the conductivity ratios.
Resistance is accumulated slice by slice; if the cumulative value ever exceeds
`RTH_LIMIT` (100&nbsp;°C/W) the calculation aborts. Once all layers are processed, the
selected cooler model—direct resistance or a convection term based on the final
area—is added to obtain the per-die and total stack resistance.

## Assumptions & Limitations

* Source length and width are entered in **millimetres**.
* Layer thicknesses use **micrometres** and conductivities are in **W/(m·K)**.
* Results (and limits) are reported in **°C/W**.
* Each layer is considered homogeneous; spreading is idealised as described
  above.
* Calculations stop if cumulative Rth exceeds 100&nbsp;°C/W.

## Contributing
 
Improvements are welcome! Fork this repo and open a pull request with your changes.
Use the **Run** button in the Apps Script editor to manually test your updates.
If you add unit tests later, run them before submitting your PR.
