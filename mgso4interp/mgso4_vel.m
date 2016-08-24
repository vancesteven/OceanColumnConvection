function vel = mgso4_vel(S,T,P)

% MGSO4_VEL    Sound Speed of MgSO4 brines
%=========================================================================
%
% USAGE:  vel = mgso4_vel(S,T,P)
%
% DESCRIPTION:
%    Sound speed of magnesium sulfate brine using 
%    Vance and Brown data
%
% INPUT: 
%   S = salinity    [molal]
%   T = temperature [Kelvin]
%   P = pressure    [bars]
%
% OUTPUT:
%   vel = Coefficient of thermal expansion  [km/s]
%
% Definition:
%   vel = sound speed [km/s]
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

global velinterpolant

if isempty(velinterpolant)
    disp('Loading MgSO4 data')
    mgso4_loader
end

vel = velinterpolant(P,S,T);