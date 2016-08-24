function pden = mgso4_pden(S,T,P,PR)

% MGSO4_PDEN    Potential density
%===========================================================================
%
% USAGE:  pden = mgso4_pden(S,T,P,PR) 
%
% DESCRIPTION:
%    Calculates potential density of water mass relative to the specified
%    reference pressure by pden = mgso4_dens(S,ptmp,PR).
%   
% INPUT:  (all must have same dimensions)
%   S  = salinity    [molal ]
%   T  = temperature [K]
%   P  = pressure    [bar]
%   PR = Reference pressure  [bar]
%
% OUTPUT:
%   pden = Potential density relative to the ref. pressure [kg/m^3] 
%
% AUTHOR:  Jason Goodman (goodman_jason@wheatonma.edu) using code by
%          Phil Morgan 1992/04/06  (morgan@ml.csiro.au)
%
% REFERENCES:
%   A.E. Gill 1982. p.54
%   "Atmosphere-Ocean Dynamics"
%   Academic Press: New York.  ISBN: 0-12-283522-0
%=========================================================================

% CALLER:  general purpose
% CALLEE:  mgso4_ptmp.m mgso4_dens.m

ptmp = mgso4_ptmp(S,T,P,PR);
pden = mgso4_dens(S,ptmp,PR);

return      
%=========================================================================

