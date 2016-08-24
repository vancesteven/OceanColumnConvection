function T_K = mgso4_t_freezing(S,P)

% MGSO4_T_FREEZING    Freezing temperature of aqueous MgSO4 
%=========================================================================
%
% USAGE:  T_K = mgso4_t_freezing(S,P)
%
% DESCRIPTION:
%    Freezing temperature of aqueous magnesium sulfate using interpolation on Vance
%    et al. (2014) data
%
% INPUT: 
%   S = salinity    [molal]
%   T = temperature [Kelvin]
%   P = pressure    [bars]
%
% OUTPUT:
%   dens = density  [kg/m^3] 
% 
% AUTHOR:  Steve Vance (svance@jpl.nasa.gov)
%
% REFERENCES:  
% Marion et al (2005), Geochimica et Cosmochimica Acta, Vol. 69, No. 2, pp. 259?274
% Vance and Brown (2005), doi:10.1016/j.icarus.2005.06.005
% Vance and Brown (2013), Geochimica et Cosmochimica Acta, 110:176?189
% Vance et al. (2014), Planetary And Space Science, 96:62?70
%=========================================================================

% CALLER: general purpose
% CALLEE:  mgso4_loader polyvaln

global tfreezeinterpolant

if isempty(tfreezeinterpolant)
    disp('Loading MgSO4 data')
    mgso4_loader
end

T_K = tfreezeinterpolant(P,S);