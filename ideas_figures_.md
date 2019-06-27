## Molecules
- Iopamidol
- Polymer
- Monomer

## Main Points
- This approach could be applied to molecules with a single peak, for example: Ioxilan (Joey did it) but we don't do it here.
- The current approach does not necessarily apply to experimental conditions beyond those used durign the calibration
- The ML approach is more robust to varibel experimental conditions.
- Classification could be enough for many clinical applications.
- The image processing time is much faster because it does not require complicated Bloch fitting.
- ML can extended the range of pH because it dependson the total spectral features and not only on being able to observer distinct peaks.

## Theme 01: Previous vs new 
1. Estimate Lorentzian properties of @4.2 and @4.6 at **fixed** exp. conditions
2. Build calibration curve between peak ratio and pH by linear regression
3. Apply calibration to all conditions
4. Mesure performance

## Theme 02: New way with Lorentzian
1. Estimate Lorentzian properties of @4.2 and @4.6 at *ALL* exp. conditions
2. Split data into training and testing (70%, 30%)
2. Build calibration curve between peak ratio and pH by linear regression
3. Apply calibration to all conditions
4. Mesure performance

## Theme 03: New way without Lorentzian
1. Split data into training and testing (70%, 30%)
2. Learn dimensionality reduction (PCA)
3. Build calibration curve between reduced spectra and pH by linear regression
4. Apply calibration to all conditions
5. Mesure performance

# Figures

## Figure 01
- Principles of ML. Simlar to Joey's figure

## Figure 02
- Effect of pH on the "shape" of a Z spectra for Iopamidol
  - Include other molecules in supplemental info.
## Figure 03
x = measured pH
y = Predicted pH
color coded by concentration
- pH prediction performance for Iopamidol Regression (RMSE = metrics)
  - A) 

