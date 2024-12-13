%***** RUN 2D MODEL FROM IMAGE ***********************************
% This section defines key parameters and runs the heat transport model
% for a geological cross-section defined by an input image.

% Indicate whether this is a test run.
% TEST = 'YES' will use simplified, constant coefficients to validate the model.
TEST = 'NO';

% Define grid sizes for the x-direction to test spatial convergence.
NN = [100, 200, 400]; % Number of grid cells for different resolutions.
convergence = 'Space'; % Indicates the convergence type (e.g., Space).

% Define the CFL condition parameter to ensure numerical stability.
CFL = 0.5; % The CFL number governs the time step size for the simulation.

% Loop through the grid resolutions to simulate the heat transport model
for nn = 1:length(NN)

    % Define model domain and grid parameters
    W = 16e3;                      % Domain width (horizontal extent) in meters.
    Nx = NN(nn);                   % Number of grid cells in the x-direction.
    h = W / Nx;                    % Grid spacing (cell size) in meters.
    n_units = 9;                   % Number of geological rock units defined in the image.

    % Load the geological model from the image and interpolate to the grid.
    try
        [units, D, Nz] = ModelFromImage('section.tiff', n_units, W, Nx);
    catch ME
        % Handle errors if the model setup cannot be loaded from the image.
        error('Error loading model from image: %s', ME.message);
    end

    %***** Define Sediment Parameters ***********************************

    % Define densities and porosities for different sediment types.
    rho_particle = 2650;  % Particle density of solid sediments (kg/m^3).
    rho_air = 1.225;      % Air density (kg/m^3).
    sa_phi = 0.253;       % Porosity of sand (fraction).
    gr_phi = 0.325;       % Porosity of gravel (fraction).
    si_phi = 0.17;        % Porosity of silt (fraction).

    % Thermal conductivities for different materials.
    sa_particle_kT = 8;   % Sand thermal conductivity (W/m·K).
    gr_particle_kT = 5;   % Gravel thermal conductivity (W/m·K).
    si_particle_kT = 3.5; % Silt thermal conductivity (W/m·K).
    air_kT = 0.025;       % Air thermal conductivity (W/m·K).

    % Calculate the bulk density of sediments, accounting for air in pores.
    rho_sand = (1 - sa_phi) * rho_particle + (sa_phi * rho_air); % Bulk density of sand.
    rho_grav = (1 - gr_phi) * rho_particle + (gr_phi * rho_air); % Bulk density of gravel.
    rho_silt = (1 - si_phi) * rho_particle + (si_phi * rho_air); % Bulk density of silt.

    % Define material properties for each rock unit
    matprop = [
                 % unit conductivity(kT)     density(rho0)   sp. heat capacity(Cp)       heat production(Hr)
        1               3.678                   2697.6              845                         4.172e-6           % HE1
        2               2.467                   2750                775                         2.9e-6             % Gneiss
        3               3.218                   2703.5              845                         5.575e-6           % HE2
        4               0.272                   rho_sand            830                         1e-6               % Sand
        5               1.075                   rho_grav            1000                        1e-6               % Gravel
        6               1.3                     2000                878                         1e-6               % Clay (sea)
        7               2.49                    rho_silt            1000                        1e-6               % Silt
        8               0.61                    2000                1510                        1e-6                % Mud (Sea)
        9               1e-6                    rho_air             1012                        0                   % Air/water
    ];

   %***** Generate Coefficient Fields **********************************
    % Define thermal properties for the simulation based on rock units.
    switch TEST
        case 'YES' % Constant coefficient variant for testing.
            rho0 = 2400 * ones(Nz, Nx); % Constant density field.
            Cp = 1000 * ones(Nz, Nx);  % Constant specific heat capacity.
            sigma = 1000 * ones(Nz, Nx); % Constant thermal conductivity.
            Hr = ones(Nz, Nx);         % Uniform radiogenic heating.
            air = units == 9;          % Identify air cells.
            t_gradient = 35 / 1000;    % Geothermal gradient (°C/m).
            geotherm = [0, t_gradient]; % Define the geothermal profile.
        case 'NO' % Spatially varying coefficients based on material properties.
            rho0 = reshape(matprop(units, 3), Nz, Nx); % Density field (kg/m^3).
            Cp = reshape(matprop(units, 4), Nz, Nx);  % Heat capacity (J/kg·K).
            sigma = reshape(matprop(units, 2), Nz, Nx); % Thermal conductivity (W/m·K).
            Hr = reshape(matprop(units, 5), Nz, Nx);  % Radiogenic heating (W/m^3).
            air = units == 9;          % Identify air cells.
            t_gradient = 35 / 1000;    % Geothermal gradient (°C/m).
            geotherm = [0, t_gradient]; % Define the geothermal profile.
    end

    % Calculate volumetric heat capacity and normalized conductivity.
    a = rho0 .* Cp;                    % Volumetric heat capacity (J/m^3·K).
    kT0 = sigma * 1000 ./ a;           % Normalized thermal conductivity (s·K/m^2).

    %***** Set Simulation Parameters ************************************
    Ttop = 8;                          % Surface temperature (°C).
    Tbot = Ttop + D * geotherm(2);     % Bottom temperature (°C).
    cT = 1e-9;                         % Conductivity temperature dependence.
    mT = 2;                            % Conductivity power-law exponent.
    g0 = 9.8;                          % Gravitational acceleration (m/s^2).
    aT = 1e-4;                         % Thermal expansivity (1/°C).
    yr = 60 * 60 * 24 * 365.25;        % Seconds per year (s).
    tend = 1e6 * yr;                   % Total simulation time (1 million years in seconds).
    nop = 100;                         % Output every 'nop' steps for visualization.

    %***** Run the Main Model *******************************************
    % The simulation solves heat transport using the initialized parameters.
    try
        run('./helmsdale_take_2.m'); % Call the heat transport model script.
    catch ME
        % Error handling for issues in the main model script.
        error('Error running the main model script: %s', ME.message);
    end
end

%***** MAKEFIG FUNCTION ***********************************
% Function to plot the temperature distribution
function makefig(xc, zc, T)
    % Convert simulation time to years
    timeInYears = t / yr;

    % Plot the temperature field as a color map.
    imagesc(xc, zc, T);
    axis equal; % Ensure uniform axis scaling.
    colorbar; % Add color bar to indicate temperature scale.
    hold on;

    % Overlay contour lines for 100°C and 150°C.
    contour(xc, zc, T, [100, 150], 'k', 'LineWidth', 1.5);

    % Label the axes and add a title.
    xlabel('Horizontal Distance [m]', 'FontSize', 15);
    ylabel('Depth [m]', 'FontSize', 15);
    title('Temperature Distribution', 'FontSize', 18);

    % Annotate the plot with the current simulation time in years.
    annotationText = sprintf('Time: %.2f years', timeInYears);
    text(0.02, 0.98, annotationText, 'Units', 'normalized', ...
        'FontSize', 14, 'Color', 'w', 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'k', 'EdgeColor', 'w');
end
