function beta = mgso4_beta(S,T,P)

% MGSO4_beta    Coefficient of haline contraction of MgSO4 brines
%=========================================================================
%
% USAGE:  beta = mgso4_beta(S,T,P)
%
% DESCRIPTION:
%    Haline contraction coefficient of magnesium sulfate brine using 
%    Vance & Brown data
%
% INPUT: 
%   S = salinity    [molal]
%   T = temperature [Kelvin]
%   P = pressure    [bars]
%
% OUTPUT:
%   beta = Coefficient of haline contraction  [1/mol] 
%
% Definition:
%   beta = - (1/rho) d(rho)/ds at constant P,T.
% 
% AUTHOR:  Jason Goodman (goodman_jason@wheatonma.edu)
%
% REFERENCES:  
% Marion et al (2005), Geochimica et Cosmochimica Acta, Vol. 69, No. 2, pp. 259?274
% Vance and Brown (2005), doi:10.1016/j.icarus.2005.06.005
% Vance and Brown (2013), Geochimica et Cosmochimica Acta, 110:176?189
%=========================================================================

% CALLER: general purpose
% CALLEE:  mgso4_loader polydern polyvaln

global betainterpolant

if isempty(betainterpolant)
    disp('Loading MgSO4 data')
    mgso4_loader
end

beta = betainterpolant(P,S,T);