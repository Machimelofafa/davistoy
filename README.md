# Stack-up Rth Calculator

This repository contains a Google Apps Script web app that estimates the thermal resistance (Rth) of a material stack-up. The UI is delivered from `index.html` and the logic resides in `Code.gs`.

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
3. Click **Run** to calculate thermal resistance.
4. Use **Monte Carlo** to simulate uncertainty; adjust iterations and deviation in the Monte Carlo section.
5. Hover over form labels to view tooltips explaining the parameters.
6. Tooltips on the results table and graphs explain the values shown.

## Contributing

Improvements are welcome! Fork this repo and open a pull request with your changes.
Use the **Run** button in the Apps Script editor to manually test your updates.
If you add unit tests later, run them before submitting your PR.

