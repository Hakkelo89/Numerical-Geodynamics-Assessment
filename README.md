# Numerical-Geodynamics-Assessmen
# 2D Heat Transport Model

This repository contains MATLAB scripts for simulating and analyzing 2D heat transport in a geological setting. The model utilizes advection-diffusion equations solved with a fourth-order Runge-Kutta (RK4) time integration scheme. Input data, including geological units, is derived from an image file, and material properties are parameterized for various rock and sediment types.

---

## Files Included

### 1. `helms_parameters2.m`
- **Purpose**: Sets up the parameters, material properties, and simulation settings for the heat transport model. This script prepares the input grid, computes material properties, and calls the main simulation script `helmsdale_take_2.m`.
- **Key Features**:
  - Defines material properties (density, thermal conductivity, heat capacity) for each geological unit.
  - Supports different grid resolutions for spatial convergence testing.
  - Configurable geothermal gradients and radiogenic heat production.
  - Handles both constant and variable coefficient simulations (`TEST = 'YES'` or `TEST = 'NO'`).
  - Plots and displays simulation results.

### 2. `helmsdale_take_2.m`
- **Purpose**: Implements the 2D advection-diffusion heat transport model using the RK4 method for time integration. It computes the heat transport and outputs numerical error norms for convergence validation.
- **Key Features**:
  - Initializes the computational grid and temperature field.
  - Solves heat transport equations iteratively until the simulation end time.
  - Supports custom visualization of temperature fields, including contour plotting.
  - Includes numerical error analysis for model accuracy.

---

## Prerequisites

- MATLAB (R2020b or later recommended).
- A valid geological image file (`section.tiff`) that represents the cross-sectional domain.
- Knowledge of MATLAB scripting for parameter customization.

---

## Usage

### Step 1: Prepare Input Image
- Ensure the geological input image file (`section.tiff`) is located in the same directory as the scripts.
- This image should define the spatial distribution of geological units.

### Step 2: Customize Parameters
- Open `helms_parameters2.m` and configure the following:
  - **TEST mode**: Set `TEST = 'YES'` for simplified testing with constant coefficients, or `TEST = 'NO'` for variable coefficients.
  - **Grid Resolutions**: Modify `NN` to define grid sizes for convergence testing.
  - **Geothermal Gradient**: Adjust `t_gradient` to set the gradient of the temperature with depth.

### Step 3: Run the Simulation
- Execute `helms_parameters2.m` in MATLAB.
  ```bash
  >> helms_parameters2
  ```
- The script will process the input image, compute material properties, and call `helmsdale_take_2.m` to solve the heat transport model.

### Step 4: View Results
- Temperature distributions and contours will be plotted at regular intervals.
- Numerical error norms (`Errx` and `Errz`) will be displayed at the end of the simulation for convergence analysis.

---

## Output

### Plots
- **Temperature Distribution**: Visualized as a 2D color map with temperature contours (e.g., 100°C and 150°C).
- **Time Annotation**: Displays the elapsed simulation time in years on the plot.

### Numerical Error
- **Error Norms**:
  - `Errx`: Normalized error in the x-direction.
  - `Errz`: Normalized error in the z-direction.

---

## Function Overview

### `helms_parameters2.m`
- **Main Script**:
  - Reads and interpolates the geological model from the image file.
  - Sets up simulation parameters and material properties.
  - Calls `helmsdale_take_2.m` to run the simulation.

### `helmsdale_take_2.m`
- **Main Simulation Code**:
  - Solves the 2D advection-diffusion equation for heat transport.
  - Implements RK4 for time integration.
  - Visualizes the temperature field during the simulation.
  - Outputs numerical error norms for validation.

### Utility Functions
1. **`makefig`**:
   - Plots the temperature field as a color map with contours.
   - Annotates the plot with the simulation time in years.
2. **`diffusion`**:
   - Computes the heat flux divergence and radiogenic heat source term.
   - Returns the rate of temperature change at each grid cell.

---

## Customization

- **Material Properties**:
  - Update the `matprop` table in `helms_parameters2.m` to modify density, thermal conductivity, heat capacity, or heat production for each rock unit.
- **Boundary Conditions**:
  - Adjust the geothermal gradient and surface temperature in `helms_parameters2.m`.
- **Grid Size**:
  - Modify the `NN` array to define the number of grid cells for spatial convergence studies.

---

## Example Input and Output

### Input
- Geological cross-section (`section.tiff`) with 9 rock units.
- Simulation parameters:
  - Surface temperature: `Ttop = 8°C`.
  - Geothermal gradient: `t_gradient = 35°C per km`.

### Output
- Temperature distribution plot with key contours.
- Error norms for numerical validation.

---

## Contact

For questions or support, contact with me or open an issue in this repository.
