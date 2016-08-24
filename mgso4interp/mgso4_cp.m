function cp = mgso4_cp(S,T,P)

% MGSO4_CP    Heat capacity of MgSO4 brines
%=========================================================================
%
% USAGE:  Cp = mgso4_cp(S,T,P)
%
% DESCRIPTION:
%    Heat capacity at constant pressure of magnesium sulfate brine using 
%    Vance & Brown data
%
% INPUT: 
%   S = salinity    [molal]
%   T = temperature [Kelvin]
%   P = pressure    [bars]
%
% OUTPUT:
%   cp = Heat capacity  [J/kg-K] 
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

global cpinterpolant

if isempty(cpinterpolant)
    disp('Loading MgSO4 data')
    mgso4_loader
end

cp = cpinterpolant(P,S,T);