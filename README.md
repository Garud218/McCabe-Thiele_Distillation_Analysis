# McCabe-Thiele Distillation Analysis

## Overview
This project is a course assignment for CHE213 (Chemical Engineering) under the guidance of Prof. Soumik Das. It involves designing and analyzing a distillation column for separating a methanol-water mixture using the McCabe-Thiele graphical method, incorporating a side stream withdrawal.

**Author:** Tushar Verma<br>
**Date:** Mar'24 - Apr'24  

## Aim
To determine the number of ideal trays required for the separation of a methanol-water mixture in a distillation column with a side stream, as well as the locations of the feed tray and the side stream withdrawal tray, based on the given specifications.

## Objectives
- Designed an advanced distillation column framework using the McCabe-Thiele graphical method for multi-section columns with side stream extraction.
- Developed precise vapor-liquid equilibrium curves by fitting experimental VLE data with advanced curve-fitting techniques.
- Optimized operating line construction to determine ideal reflux ratios, enhancing process efficiency and separation performance.

## Materials Required
- MATLAB software for executing the scripts.
- Provided VLE data for methanol-water binary mixture.

## Procedure
1. Fit the experimental VLE data to the equation \( y = \frac{a x}{1 + b x + c x^2} \) using nonlinear least squares.
2. Solve overall mass and component balances to find distillate and bottoms flow rates.
3. Calculate the minimum reflux ratio using the pinch point at the side stream.
4. Set the actual reflux ratio to 2.5 times the minimum.
5. Divide the column into three sections (above side stream, between side stream and feed, below feed) and compute flow rates and operating line equations for each.
6. Plot the McCabe-Thiele diagram including equilibrium curve, operating lines, feed line, and side stream line.
7. Implement staircase construction to count theoretical stages and identify feed and side stream locations.

## Experiments and Observations
The project uses computational methods in MATLAB to model the distillation process. Key calculated values include:

- Fitted VLE parameters: \( a = 7.9673 \), \( b = 9.4109 \), \( c = -2.4383 \)
- Distillate (D): 190.53 kmol/hr at 97% methanol
- Bottoms (W): 259.47 kmol/hr at 2% methanol
- Side stream (S): 50 kmol/hr at 70% methanol
- Minimum reflux ratio: 0.57
- Actual reflux ratio: 1.42
- Feed line intersection: (0.3829, 0.7185)
- Side stream pinch: 0.8724
- Section II intersection: (0.3959, 0.6666)

The staircase construction yields linear trends in operating lines with high accuracy in intersections.

Raw data and equations are in the problem statement; scripts handle all computations.

## Analysis and Results
- The VLE curve fits well to the data, enabling accurate equilibrium predictions.
- Operating lines for the three sections are calculated with slopes based on flow ratios: Section I (m1 ≈ 0.586), Section II (m2 ≈ 0.560), Section III (m3 ≈ 1.239).
- The reaction follows the specified separations with the given feed (500 kmol/hr, 45% methanol, 80% liquid).
- Calculated number of theoretical stages: 10.82 (approximately 11 ideal trays).
- Side stream withdrawal: Tray 4 (from the top).
- Feed tray: Tray 7 (from the top).
- The method confirms efficient separation with optimized reflux.

## Possible Sources of Errors
- Numerical inaccuracies in solving nonlinear equations or curve fitting.
- Assumptions of constant molar overflow and ideal VLE behavior.
- Discretization in plotting and staircase construction.
- Potential wind or environmental factors if experimental, but here computational.

## Repository Structure
- `Problem Statement.pdf`: Detailed problem description, VLE data, and requirements.
- `script_lab4.m`: Primary script for VLE fitting, mass balances, operating lines, plotting, and stage calculation.
- `Feedeqll.m`: Helper function for solving feed line and equilibrium intersection.
- `LAB4.m`: Supplementary script for equilibrium curve plotting, operating lines, and stage counting.
- `feed_eq.m`: Alternative helper function for feed equation.

## How to Replicate
1. Open MATLAB and load the scripts.
2. Run `script_lab4.m` or `LAB4.m` to perform calculations and generate plots.
3. Inspect the output for flow rates, reflux ratios, and number of stages.
4. For Python alternative, use NumPy and SciPy to replicate the fitting (curve_fit), solving (fsolve, root_scalar), and computations as demonstrated in equivalent code.

For further details, refer to the attached PDF and MATLAB scripts.
