% remember to CLEAR ALL if starting from scratch
usesparse = true;

%% seawater eos and thermodynamics
str_EOS = 'gsw302';
swEOS = swEOS_chooser(str_EOS); % 'mgso4' or 'gsw302' or 'gsw305

%% model properties
dz = 2500; % grid spacing in m
Nz = 40; %total number of vertical volumes

% deltat = 5*24*3600; % time step in s
deltat = 5*3.16e7;
tend = 14e6*3.16e7; % total simulationtend = 14e6*3.16e7; % total simulation time in s
 
%% planetary properties
g = 1.3;
Tsurf = 102;  % McFadden et al. 2007
%Hseafloor = 0.075;  % Chosen to give equilibrium ice thickness of 10 km, using conduction
                    % equation below.H
%  Hseafloor = 0.05;  % thicker: W m-2 minimum (larger than radiogenic only value, but equivalent to a 12.5 km thick ice shell---capturestidal heating above mantle in that sense)
%    Hseafloor = 0.10;  % nominal: W m-2 minimum (larger than radiogenic only value, but equivalent to a 12.5 km thick ice shell---capturestidal heating above mantle in that sense)
%   Hseafloor = 0.05;  % anomaly: W m-2 minimum (larger than radiogenic only value, but equivalent to a 12.5 km thick ice shell---capturestidal heating above mantle in that sense)
%Hseafloor = 0.190;  % W m-2 peak (tidal + radiogenic; O'Brien et al. 2002)
Hseafloor = 0.375;  % Chosen to correspond to 100 MW plume heating an area 
                    % 20 km in diameter
% turned on .375 flux at t=1363399 after running at 0.019

% H ramping over time
DHDT = 0; % make heat increase by a factor of 10 over 125Myr due to changing eccentricity(Hussmann et al. 2004)
t_hosc=125e6*3.16e7;


%S0 = 0.01;  % molal.  
% S0 = 0.3;  % molal.  Roughly 35 g/kgish
S0 = 2*gsw_SSO; % standard seawater salinity of 36.1650 in g/kg

h0=6.2e3;  % Initial ice thickness

% Sfluxseafloor = 0;

%  Sfluxseafloor = 3e-8;
Sfluxseafloor = 1.9e-7   % mol/(m^2 s) % Chosen to give 0.06 moles/sec 
                          % delivered over an area 20 km in diameter
%% averaging properties
averaging_time = 1000*3.16e7; % Average properties over n years
%averaging_time = deltat; % Don't average`

%% continuation
continuecurrentrun  = 0; % set this to restart after using control-c to cancel the previous run
loadfromfile = 1; % load a start.mat file stored in the same directory

%% run the model
column_convect_implicitmodel

%% plot the results 
%  plotlist = {'T','S','rho','nsquared','q','h','v','dz0dz'};
 plotlist = {'T','S','rho','q','h','v','c','kt'};
% plotlist = {'T','c','nsquared','h'}
% surf_type = 'imgsc';
surf_type = 'sf';
figure(23); clf
plot_column_convection

 %% save vectors of T, S, and h as a starting file
 save('start.mat','T','S','h'); % create a start.mat file
