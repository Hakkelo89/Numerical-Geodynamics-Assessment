%***** RUN 2D MODEL FROM IMAGE ***********************************

% clear workspace
clear all; close all; %clc;

TEST = 'real'; % choose from 'calibration' and 'real'

% load model setup from image, interpolate to target grid size
W       = 16e3;     % domain width (must correspond to width of image) [m]
Nx      = 200;      % target grid size z-direction
h       = W/Nx;     % grid spacing based on image width and target grid size
n_units = 9;        % number of rock units contained in image
[units,D,Nz] = ModelFromImage('section.tiff',n_units,W,Nx);

% material properties for each rock unit (update based on your calibration)

% calculate sediment parameters
% set variables (values found online)
rho_particle = 2650; % density of sand, gravel & silt
rho_air  = 1.225;    % density of air
sa_phi   = 0.253;    % sand porosity fraction
gr_phi   = 0.325;    % gravel porosity fraction
si_phi   = 0.17;     % silt porosity fraction

sa_particle_kT = 8;      % sand particle thermal conductivity 
gr_particle_kT = 5;      % gravel particle thermal conductivity
si_particle_kT = 3.5;    % silt particle thermal conductivity
air_kT         = 0.025;  % air thermal conductivity

% calculate bulk density 
rho_sand   = sa_phi * rho_air + (1-sa_phi) * rho_particle;
rho_gravel = gr_phi * rho_air + (1-gr_phi) * rho_particle;
rho_silt   = si_phi * rho_air + (1-si_phi) * rho_particle;


matprop = [
% unit conductivity(kT) density(rho)  heat capacity(Cp)  heat production(Hr)
   1	      3.678       2697.6	       845	              4.172         %HE1
   2	      1	          2000	           1000	              1             %Gneiss
   3	      1	          rho_sand         830	              1             %Sand 1
   4	      3.218       2703.5	       845	              5.575         %HE2
   5	      1	          rho_gravel       1000	              0.6           %Gravel
   6	      1	          2000	           1381	              2.75          %Clay
   7	      1	          rho_silt         1000	              1.4           %Silt
   8	      1	          2000	           2512	              3.5           %Mud
   9	      1e-6        rho_air	       1000	              0];           % air/water

switch TEST

    case 'real'

    % get coefficient fields based on spatial distribution of rock units from image
    % pay attention if any unit conversion is required!
    rho0    = reshape(matprop(units,3),Nz,Nx);
    Cp      = reshape(matprop(units,4),Nz,Nx);
    sigma   = reshape(matprop(units,2),Nz,Nx);
    Hr      = reshape(matprop(units,5),Nz,Nx);

    case 'calibration'
    % calibration. set uniform parameters 
    rho0 = 2400*ones(Nz,Nx);
    Cp   = 1000*ones(Nz,Nx);
    sigma  = ones(Nz,Nx);
    Hr   = ones(Nz,Nx);

end

kappa = sigma ./ rho0 .* Cp;

air = units == 9;
geotherm = 35/1000;

%Test in constant coefficient unit value case (repeat for other variables)
%rho = 2400*ones(Nz,Nx);

% continue setting remaining model parameters, then call model routine

%*****  RUN MODEL
%run('./ModelFromImage.m');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TGrad = 35;                    % Reference Temperature Gradient from National Survey
Ttop  = 5;                     % surface temperature
Tbot  = Ttop + TGrad*(D/1000); % Temperature at the bottom of the domain using reference gradient


%rho0  = 1000;         % reference density [kg/m3]
%kT0   = 1e-7;         % heat diffusivity [m2/s]
cT    = 1e-9;         % kT T-dependence prefactor
mT    = 2;            % kT T-dependence powerlaw
g0    = 9.8;          % gravity [m/s2]
aT    = 1e-4;         % thermal expansivity
ADVN  = 'WENO5';      % advection scheme ('UPW1', 'CFD2', 'UPW3', 'WENO5')

yr    = 3600*24*365;  % seconds per year [s]
tend  = 3e6*yr;       % stopping time [s]
CFL   = 0.95;         % Time step limiter. Decreasing this increases the timestep. Less accuracy, quicker
nop   = 100;          % output figure produced every 'nop' steps
dTdto = 0;

run('./helmsdale_take_2.m');