function dens = mgso4_dens(S,T,P)

% MGSO4_DENS    Density of MgSO4 brine
%=========================================================================
%
% USAGE:  dens = mgso4_dens(S,T,P)
%
% DESCRIPTION:
%    Density of magnesium sulfate brine using interpolation on Vance & Brown data
%
% INPUT: 
%   S = salinity    [molal]
%   T = temperature [Kelvin]
%   P = pressure    [bars]
%
% OUTPUT:
%   dens = density  [kg/m^3] 
% 
% AUTHOR:  Jason Goodman (goodman_jason@wheatonma.edu)
%
% REFERENCES:  
% Marion et al (2005), Geochimica et Cosmochimica Acta, Vol. 69, No. 2, pp. 259?274
% Vance and Brown (2005), doi:10.1016/j.icarus.2005.06.005
% Vance and Brown (2013), Geochimica et Cosmochimica Acta, 110:176?189
%=========================================================================

% CALLER: general purpose
% CALLEE:  mgso4_loader polyvaln

global rhointerpolant

if isempty(rhointerpolant)
    disp('Loading MgSO4 data')
    mgso4_loader
end

dens = rhointerpolant(P,S,T);