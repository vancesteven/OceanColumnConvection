function adtg = mgso4_adtg(S,T,P)

% MGSO4_ADTG    Adiabatic temperature gradient of MgSO4 brines
%=========================================================================
%
% USAGE:  adtg = mgso4_adtg(S,T,P)
%
% DESCRIPTION:
%    Adiabatic temperature gradient of magnesium sulfate brine using 
%    Vance & Brown data
%
% INPUT: 
%   S = salinity    [molal]
%   T = temperature [Kelvin]
%   P = pressure    [bars]
%
% OUTPUT:
%   adtg = Adiabatic temp gradient  [K/bar] 
%
% Definition:
%   adtg = alpha T / (rho Cp)
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

global adtginterpolant

if isempty(adtginterpolant)
    disp('Loading MgSO4 data')
    mgso4_loader
end

adtg = adtginterpolant(P,S,T);