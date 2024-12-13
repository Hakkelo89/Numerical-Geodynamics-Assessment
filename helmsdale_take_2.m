%***** 2D ADVECTION DIFFUSION MODEL OF HEAT TRANSPORT *******************
% This script simulates the 2D heat transport using advection-diffusion
% equations and solves it numerically using the Runge-Kutta 4th-order (RK4) method.

%***** Initialise Model Setup ******************************************

% Define grid coordinates for cell centers and faces based on grid spacing.
xc = h/2:h:W-h/2; % x-coordinate vector for cell center positions [m].
zc = h/2:h:D-h/2; % z-coordinate vector for cell center positions [m].
xf = 0:h:W; % x-coordinate vector for cell face positions [m].
zf = 0:h:D; % z-coordinate vector for cell face positions [m].

% Create 2D coordinate arrays for the grid.
[Xc, Zc] = meshgrid(xc, zc); % Meshgrid defines the full grid in x and z directions.

% Insulating boundaries:
ix3 = [1, 1:Nx, Nx]; % Indices for insulating left and right sides of the domain.
iz3 = [1, 1:Nz, Nz]; % Indices for insulating top and bottom boundaries.

% Set the initial temperature distribution:
T = Ttop + geotherm(2) .* Zc; % Initialize with a linear geothermal gradient.
T(air) = Ttop; % Set temperature of air cells to the surface temperature.

% Initialize material properties:
rho = rho0 .* (1 - aT .* (T - Ttop)); % Density field accounting for thermal expansion.
kT = kT0 .* ones(Nz, Nx); % Thermal conductivity field (initialized to uniform values).

% Initialize the figure to visualize the initial condition:
figure(1); clf; % Create and clear the figure window.
makefig(xc, zc, T); % Plot the initial temperature distribution.

%***** Solve Model Equations *******************************************

% Compute the initial time step using the CFL condition:
dt = CFL * (h/2)^2 / max(kT(:)); % CFL ensures numerical stability.

% Initialize time-stepping variables:
t = 0; % Current simulation time [s].
k = 0; % Time step counter.
dTdt = 0; % Rate of temperature change.

% Time-stepping loop: Solve the equations until the simulation end time.
while t <= tend
    % Increment the time and step counter:
    t = t + dt; % Update time by adding the time step.
    k = k + 1; % Increment the step count.

    % Store the previous temperature field and rate of change for RK4:
    dTdto = dTdt; % Save the previous rate of change.
    To = T; % Save the current temperature field.

    % Calculate the temperature rate of change using RK4 stages:
    dTdt1 = diffusion(T, kT, h, ix3, iz3, geotherm, Hr, rho, Cp); % Stage 1.
    dTdt2 = diffusion(T + dTdt1 / 2 * dt, kT, h, ix3, iz3, geotherm, Hr, rho, Cp); % Stage 2.
    dTdt3 = diffusion(T + dTdt2 / 2 * dt, kT, h, ix3, iz3, geotherm, Hr, rho, Cp); % Stage 3.
    dTdt4 = diffusion(T + dTdt3 * dt, kT, h, ix3, iz3, geotherm, Hr, rho, Cp); % Stage 4.

    % Add the heat source term (Hr scaled to Watts).
    Hs = (Hr * 10^-6) ./ (rho .* Cp);

    % Update the temperature field using the RK4 method:
    T = T + (dTdt1 + 2 * dTdt2 + 2 * dTdt3 + dTdt4) / 6 * dt + Hs;

    % Enforce boundary condition: Set air cells to surface temperature.
    T(air) = Ttop;

    % Analytical solution for validation (optional).
    Tana = T + geotherm(2) .* Zc;

    % Plot the model's progress every 'nop' time steps:
    if ~mod(k, nop)
        makefig(xc, zc, T); % Visualize the current temperature distribution.
    end
end

%***** Calculate Numerical Error Norms *********************************
% Evaluate the error between the numerical solution and analytical solution.
Errx = norm(T - Tana, 2) ./ norm(Tana, 2); % Error in x-direction.
Errz = norm(T - Tana, 1) ./ norm(Tana, 1); % Error in z-direction.

% Display the errors and numerical scheme used:
disp(' ');
disp('Time integration scheme: RK4');
disp(['Numerical error in x = ', num2str(Errx)]);
disp(['Numerical error in z = ', num2str(Errz)]);
disp(' ');

%***** Utility Functions ************************************************

% Function to visualize the temperature field:
function makefig(x, z, T)
    % Plot the temperature field as a color map:
    imagesc(x, z, T); axis equal; % Display the temperature field and set aspect ratio.
    c = colorbar; % Add a color bar to indicate temperature values.
    hold on;

    % Add contour lines for key temperature levels:
    contour(x, z, T, [100, 150, 150], 'k'); % Contours for 100°C and 150°C.

    % Highlight specific contours with red lines and labels:
    [C, h] = contour(x, z, T, [150, 150], 'r', 'LineWidth', 2); % 150°C in red.
    clabel(C, h, 'FontSize', 12, 'Color', 'r'); % Optional: Add labels to the 150°C contour.

    [C, h] = contour(x, z, T, [100, 100], 'r', 'LineWidth', 2); % 100°C in red.
    clabel(C, h, 'FontSize', 12, 'Color', 'r'); % Optional: Add labels to the 100°C contour.

    % Label the axes and title the figure:
    xlabel('Horizontal Distance [m]', 'FontSize', 15); % Label for x-axis.
    ylabel('Depth [m]', 'FontSize', 15); % Label for y-axis.
    ylabel(c, 'Temperature [°C]', 'FontSize', 15); % Label for the color bar.
    title('Temperature Distribution', 'FontSize', 18); % Add a title to the plot.

    drawnow; % Refresh the plot to update changes.
end

% Function to calculate diffusion rate:
function [dTdt] = diffusion(f, k, h, ix, iz, geotherm, Hr, rho, Cp)
    % Compute thermal conductivity at cell faces:
    kx = (k(:, ix(1:end-1)) + k(:, ix(2:end))) / 2; % x-direction conductivity.
    kz = (k(iz(1:end-1), :) + k(iz(2:end), :)) / 2; % z-direction conductivity.

    % Compute heat flux at cell faces:
    qx = -kx .* diff(f(:, ix), 1, 2) / h; % Heat flux in x-direction.
    qz = -kz .* diff(f(iz, :), 1, 1) / h; % Heat flux in z-direction.

    % Apply basal boundary condition for geothermal flux:
    qz(end, :) = -kz(end, :) .* geotherm(2); % Enforce geothermal gradient.

    % Calculate temperature change due to flux divergence:
    dTdt_diffusion = -(diff(qx, 1, 2) / h + diff(qz, 1, 1) / h);

    % Include radiogenic heat source:
    heat_source = Hr ./ (rho .* Cp); % Convert radiogenic heat into temperature rate.

    % Total temperature rate of change:
    dTdt = dTdt_diffusion + heat_source; % Combine diffusion and source terms.
end
% Validate against observed data
fprintf('Validating against observed data...\n');
temp_sim = interp1(zc, T, depth_obs);
MAE = mean(abs(temp_obs - temp_sim));
RMSE = sqrt(mean((temp_obs - temp_sim).^2));
R2 = 1 - sum((temp_obs - temp_sim).^2) / sum((temp_obs - mean(temp_obs)).^2);
fprintf('MAE: %.2f, RMSE: %.2f, R-squared: %.2f\n', MAE, RMSE, R2);

% Plot results
figure;
plot(temp_obs, depth_obs, 'ro', 'DisplayName', 'Observed');
hold on;
plot(temp_sim, depth_obs, 'b-', 'DisplayName', 'Simulated');
xlabel('Temperature (°C)');
ylabel('Depth (m)');
legend;
title('Observed vs. Simulated Temperatures');
set(gca, 'YDir', 'reverse');
grid on;
