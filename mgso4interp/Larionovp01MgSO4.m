function C_Sm = Larionovp01MgSO4

LarionovMgSO4 = getLarionov1984MgSO4Conductivities;
Textrap = LarionovMgSO4.T_K(1:3);
Pextrap = 0.1:100:2000;
for iT = 1:nTbs 
    k_S_m_extrap(iT,:) = interp1(LarionovMgSO4.P_MPa,LarionovMgSO4.Larionov_p01m.k_S_m(:,iT),Pextrap,'pchip','extrap');
end

[Pgg,Tgg]=meshgrid(Pextrap,Textrap);
k_S_mMgSO4p01Planet = NaN*ones(size(P_MPa)); 