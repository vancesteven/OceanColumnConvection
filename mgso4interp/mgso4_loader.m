% MGSO4_LOADER    Load MgSO4 data
%=========================================================================
%
% USAGE:  script
%
% DESCRIPTION:
%    Load MgSO4 equation of state data from .mat file.  Shared by all mgso4
%    functions to ensure data consistency.
%
% INPUT: none
%
% OUTPUT:
%   Sets values for:
%     mgso4_adtgpoly, mgso4_alphapoly, mgso4_cppoly, mgso4_rhopoly
% 
% AUTHOR:  Jason Goodman (goodman_jason@wheatonma.edu)
%
% REFERENCES:  
% Marion et al (2005), Geochimica et Cosmochimica Acta, Vol. 69, No. 2, pp. 259?274
% Vance and Brown (2005), doi:10.1016/j.icarus.2005.06.005
% Vance and Brown. (2013), Geochimica et Cosmochimica Acta, 110:176?189.
%=========================================================================

% CALLER: all other mgso4_* functions
% CALLEE: loads mgso4poly.mat

global rhointerpolant cpinterpolant...
    alphainterpolant adtginterpolant betainterpolant velinterpolant Linterpolant tfreezeinterpolant

method = 'linear';

load MgSO4_EOS_parms_2012_26_17_LT.mat
load TfreezeMgSO4_Vance2014.mat

[Pgrid,Sgrid,Tgrid] = meshgrid(Pg_MPa*10,mg,Tg_C+273.15);
[Pgl,Sgl] = ndgrid(P_bar,S_mol);
pa_to_bar = 1e5;
% calculate interpolants

Pgrid = permute(Pgrid,[2 1 3]);
Sgrid = permute(Sgrid,[2 1 3]);
Tgrid = permute(Tgrid,[2 1 3]);
rhog = permute(rhog,[2 1 3]);
Cpg = permute(Cpg,[2 1 3]);
alphag = permute(alphag,[2 1 3]);
velg = permute(velg,[2 1 3]);

rhointerpolant = griddedInterpolant(Pgrid,Sgrid,Tgrid,rhog,method);
cpinterpolant = griddedInterpolant(Pgrid,Sgrid,Tgrid,Cpg,method);
alphainterpolant = griddedInterpolant(Pgrid,Sgrid,Tgrid,alphag,method);
velinterpolant = griddedInterpolant(Pgrid,Sgrid,Tgrid,velg,method);
Linterpolant = griddedInterpolant(Pgl,Sgl,L_Jkg,method);
tfreezeinterpolant = griddedInterpolant(Pgl,Sgl,T_freeze_K,method);


%   adtg = Adiabatic temp gradient  [K/bar] 
adtgg = pa_to_bar * alphag.*Tgrid./(Cpg.*rhog);
adtginterpolant = griddedInterpolant(Pgrid,Sgrid,Tgrid,adtgg,method);

% Beta calc, centered difference
Sforbeta = [Sgrid(:,1,:) (Sgrid(:,1:end-1,:)+Sgrid(:,2:end,:))/2];
betag = (diff(rhog,1,2)./(diff(Sgrid,1,2)))./((rhog(:,1:end-1,:)+rhog(:,2:end,:))/2);
% Linearly extrapolate down minimum S
S0 = Sforbeta(1,1,1); S1 = Sforbeta(1,2,1);S2 = Sforbeta(1,3,1);
beta_at_S0 = betag(:,1,:) + (betag(:,1,:)-betag(:,2,:)).*(S1-S0)./(S2-S1);
betag = [beta_at_S0 betag];

betainterpolant = griddedInterpolant(Pgrid,Sforbeta,Tgrid,betag,method);