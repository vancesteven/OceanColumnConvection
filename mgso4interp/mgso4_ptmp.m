
function PT = mgso4_ptmp(S,T,P,PR)

% MGSO4_PTMP    Potential temperature
%===========================================================================
%
% USAGE:  ptmp = mgso4_ptmp(S,T,P,PR) 
%
% DESCRIPTION:
%    Calculates potential temperature for MgSO4 brines, using UNESCO 1983
%    method from terrestrial literature.
%   
% INPUT:  (all must have same dimensions)
%   S  = salinity    [molar]
%   T  = temperature [Kelvin]
%   P  = pressure    [bar]
%   PR = Reference pressure  [bar]
%
% OUTPUT:
%   ptmp = Potential temperature relative to PR [Kelvin]
%
% AUTHOR:  Jason Goodman (goodman_jason@wheatonma.edu) using code by
%          Phil Morgan 92-04-06  (morgan@ml.csiro.au)
%
% REFERENCES:
%    Fofonoff, P. and Millard, R.C. Jr
%    Unesco 1983. Algorithms for computation of fundamental properties of 
%    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
%    Eqn.(31) p.39
%
%=========================================================================

% CALLER:  general purpose
% CALLEE:  mgso4_adtg.m

%------
% BEGIN
%------

% theta1
del_P  = PR - P;
del_th = del_P.*mgso4_adtg(S,T,P);
th     = T + 0.5*del_th;
q      = del_th;

% theta2
del_th = del_P.*mgso4_adtg(S,th,P+0.5*del_P);
th     = th + (1 - 1/sqrt(2))*(del_th - q);
q      = (2-sqrt(2))*del_th + (-2+3/sqrt(2))*q;

% theta3
del_th = del_P.*mgso4_adtg(S,th,P+0.5*del_P);
th     = th + (1 + 1/sqrt(2))*(del_th - q);
q      = (2 + sqrt(2))*del_th + (-2-3/sqrt(2))*q;

% theta4
del_th = del_P.*mgso4_adtg(S,th,P+del_P);
PT     = th + (del_th - 2*q)/6;
  
%=========================================================================
