function L = mgso4_latent_heat_melting(S,P)

% MGSO4_LATENT_HEAT_MELTING    Latent heat of melting ice in aqueous MgSO4 
%=========================================================================
%
% USAGE:  L = mgso4_latent_heat_melting(S,P)
%
% DESCRIPTION:
%    Latent heat of aqueous magnesium sulfate using interpolation on Vance
%    et al. data, computed from the Clausius-Clapeyron relationship
%
% INPUT: 
%   S = salinity    [molal]
%   P = pressure    [bars]
%
% OUTPUT:
%   L = latent heat of melting  [J/kg/K] 
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

global Linterpolant

if isempty(Linterpolant)
    disp('Loading MgSO4 data')
    mgso4_loader
end

L = Linterpolant(P,S);