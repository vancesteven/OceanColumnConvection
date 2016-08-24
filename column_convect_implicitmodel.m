%% column_convect_implicitmodel
% solve for 1D temperature-salinity structure in an ice covered fluid with
% specified heat and salt fluxes
%
% INPUT FILE
% column_convect_run...
%
% VISUALIZATION
% plot_column_convect

%% unit conversions
patobar = 1e-5;

if ~exist('usesparse')
    error('''usesparse'' and possibly other inputs are not defined. call ''column_convect_implicit'' from input file');
end
if ~exist('DHDT')
    DHDT = 0;
else
    HS0 = Hseafloor;
end

% Ice properties
Tmelt = 273; % K
%Andersson et al. 2005
ice_cond_b1 = 632;   % W/m
ice_cond_b0 = 0.38; % W/m-K
ice_cond_c = -0.00197; % W/m-K2

%Hobbs 1974
% ice_cond_b1 = 488;   % W/m
% ice_cond_b0 = 0.468; % W/m-K
rhoi = 920; % kg/m^3;
L = 3.3e5;  % J/kg


bar2pa = 1.0137e5;                    
pscale = g*rhoi/bar2pa;
                 
% Water properties -- these are only used as a starting point
% rhow = mgso4_dens(S0,Tmelt,0); % kg/m^3
% cpw = mgso4_cp(S0,Tmelt,0);  % J/kg-K
rhow = swEOS.dens(S0,Tmelt,0); % kg/m^3
cpw = swEOS.cp(S0,Tmelt,0);  % J/kg-K


ia = 0; % iteration # within average
ja = 0; % current average being computed

if (continuecurrentrun)
    ja = length(havg);
    t = t(end);
else
    if (loadfromfile)
        disp('Loading from file');
        load start.mat T S h
        t=0;
    else
         T = ones(Nz,1)*Tmelt;
         S = ones(Nz,1)*S0;
         h = h0;
         t = 3.16e7*0;  % Start time
%         kt = zeros(size(dpint));  % Needed for first timestep
        %h = 10.12e3;
        %S = ones(Nz,1)*S0;
        %T = mgso4_ptmp(S,ones(Nz,1)*Tmelt,ones(Nz,1)*p(end),p);
        %h = 1e4;
    end
end


% heatflowparameter = (ice_cond_b1*log(Tmelt/Tsurf) + ice_cond_b0*(Tmelt - Tsurf) + ice_cond_c/2*(Tmelt^2 - Tsurf^2));
% h0 = heatflowparameter / Hseafloor;
p = ((rhoi*h*g) + rhow*g*((Nz-1):-1:0)*dz)'*patobar;  % pressure in bar
pint = ((rhoi*h*g) + rhow*g*((Nz-1.5):-1:0)*dz)'*patobar;  % pressure in bar at interfaces
dpint = diff(p);  % pressure change across interfaces;
pref = p(round(Nz/2));

% Convective parameterization properties
k_conv = .01 * 50; % Eddy velocity 1 cm/s (Goodman and Lenferenk 2012), eddy mixing length 50 m (smaller than the grid size but large enough to ensure mixing happens at the grid resolution)

% Double diffusive parameterization properties
ddflag = 1;  % Do double diffusion?
%     salt fingering
Kstar = (.01)^2; % 1 cm^2/s
Rc = 1.6; % Buoyancy ratio at which salt fingering becomes important, unitless
n = 6; % Salt fingering exponent
%     diffusive layering
lilkt = 1.3916e-7; % m^2/s  (molecular diffusivity of heat in water)

cfl = k_conv * deltat/dz^2;

% Set up averages
Saccum = zeros(size(S));
Taccum = zeros(size(T));
qaccum = zeros(size(T));  % heat flux
caccum = zeros(Nz-1,1);   % convection frequency
ktaccum = zeros(Nz-1,1);  % salt fingering frequency
ksaccum = zeros(Nz-1,1);  % diffusive layering frequency
haccum = 0;               % ice thickness



while (t < tend)
    
    t = t + deltat;
    
    %% Equation of State
    % Compute various equation-of-state properties.  Minimize number of calls
    % to swEOS library to gain speed.    

    % T, S at interfaces
    Sint = (S(2:end)+S(1:end-1))/2;
    Tint = (T(2:end)+T(1:end-1))/2;
    % density of each gridbox
%     rho  = mgso4_dens(S,T,p);
    rho  = swEOS.dens(S,T,p);

    % Thermal expansion coefficient at interfaces
    alphaint = swEOS.alpha(Sint,Tint,pint);
    % Haline expansion coefficient
    betaint = swEOS.beta(Sint,Tint,pint);
    % Adiabatic temperature gradient, K/bar
    adtgint = swEOS.adtg(Sint,Tint,pint);
    % Adiabatic temperature jump between levels
    dt_adiab = adtgint.*dpint;
    
    % use a P-m dependent melting temperature and latent heat
    Lprior = L;
    Tmeltprior=Tmelt;
    L = swEOS.latent_heat_melting(Sint(end),pint(end)); %scalar for latent heat at the ice-ocean interface
    Tmelt = swEOS.tfreezing(Sint(end),pint(end));
    if isnan(Tmelt)
       Tmelt = Tmeltprior;
    end
    if isnan(L)
        L = Lprior;
    end

        
    %% Top boundary conditions - Implicit
    
    % Heat flux through ice depends on ice thickness.  Thermal conductivity is
    % k = b1/T - b0; integrating through the ice shell gives a steady-state
    % heat flux of
    %    Htop = heatflowparameter./h;    
    % where
    heatflowparameter = (ice_cond_b1*log(Tmelt/Tsurf) + ice_cond_b0*(Tmelt - Tsurf) + ice_cond_c/2*(Tmelt^2 - Tsurf^2));
%         heatflowparameter = (ice_cond_b1*log(Tmelt/Tsurf) + ice_cond_b0*(Tmelt - Tsurf));

    % The steady-state approximation is only technically valid for timescales longer than
    % the diffusive relaxation timescale of the ice: h^2/nu = h^2 (rhoi Cp)/k =
    % millions of years.  The goal is to have an interactive ice shell that
    % allows a reasonable steady-state, not one that accurately represents
    % transients.
    
    % Since we're doing an implicit scheme, we can't calculate Htop
    % explicitly yet.
    
    % Top BC: Calculate heat flow needed to reset temperature in topmost 
    % gridpoint to Tmelt over the course of 1 timestep.  Our energy balance
    % for the ice base is:
    % heat conducted into ice = heat lost by ocean + latent heat released by freezing
    % Htop = lostbyocean + latentheatrelease
    
%     cptop = mgso4_cp(S(end),T(end),p(end));  % Heat capacity, J/kg-K
    cptop = swEOS.cp(S(end),T(end),p(end));  % Heat capacity, J/kg-K
    joulesabovefreezing = rhow*cptop*(T(end)-Tmelt)*dz;  % J/m2
    lostbyocean = joulesabovefreezing/deltat;
    
    % latentheatrelease is related to melting/freezing rate as follows:
    % latentheatrelease = (dh/dt) * rhoi * L
    %  (dh/dt)*rhoi*L = Htop - lostbyocean
    %   dh/dt = Htop/(rhoi*L) - lostbyocean/(rhoi*L)
    % Now here's where things get awkward because of the implicit scheme:
    %   (hnew - hold)/deltat = (heatflowparameter/hnew)/(rhoi*L) - lostbyocean/(rhoi*L)
    % Solve this for hnew to get a quadratic equation:
    %   hnew^2 + (lostbyocean*deltat/(rhoi*L) - hold)*hnew - heatflowparameter*deltat/(rhoi*L) = 0
    %   hnew^2 + Bh*hnew + Ch = 0    where
    Bh = (lostbyocean*deltat/(rhoi*L) - h);          % meters (really!) 
    Ch = - heatflowparameter*deltat/(rhoi*L);        % meters^2
    % Solve via quadratic equation to get new h:
    hnew = (- Bh + sqrt(Bh^2 - 4 * Ch))/2;           % meters
    
    % Upper ocean temperature change:
    surface_dtdt = (Tmelt - T(end))/deltat;    % K/s
            
    freshwaterflux = (rhoi/rhow)*(h - hnew)/deltat;   % meters of water/s
    surface_dsdt = - S(end) * freshwaterflux / dz + Sfluxsurface/(rhow*dz);   % molal/s

    h = hnew;
    
    %% Recalculate pressure field
    p = ((rhoi*h*g) + rhow*g*((Nz-1):-1:0)*dz)'*patobar;  % pressure in bar
    pint = ((rhoi*h*g) + rhow*g*((Nz-1.5):-1:0)*dz)'*patobar;  % pressure in bar at interfaces
    dpint = diff(p);  % pressure change across interfaces;
    pref = p(round(Nz/2));
    
    %% Bottom boundary conditions
    % Constant heat flux.
    
%     cpbot = mgso4_cp(S(1),T(1),p(1));  % Heat capacity, J/kg-K
    cpbot = swEOS.cp(S(1),T(1),p(1));  % Heat capacity, J/kg-K
    basal_dtdt = Hseafloor/(rhow*cpw*dz);

    % Constant salt flux (in mol/(m^2 s))
    
    basal_dsdt = Sfluxseafloor/(rhow*dz);

    
    %% Interior solution
    
    ptmp = T(2:end) - dt_adiab;
    %potential density of the water one gridbox up, if moved to this level.
%     pden = mgso4_dens(S(2:end),ptmp,p(1:end-1));
    pden = swEOS.dens(S(2:end),ptmp,p(1:end-1));
    
    % Convection if potential density of overlying water is greater than
    % density here.
    
    convecting = (pden - rho(1:end-1)) > 0;  % Convecting if potential density is increasing as we go up
    k = convecting.*k_conv;
    
    % Double diffusion
    if (ddflag)
        % Adiabatic temperature gradient (used by double diffusion)
        dtdz_adiab = (diff(T)-dt_adiab)/dz;
        % Salinity gradient
        dsdz = diff(S)/dz;

        % Density gradient ratio alpha dtdz / beta dsdz determines double
        % diffusion

        Rrho = alphaint.*dtdz_adiab./(betaint.*dsdz); % Unitless
           
        Kt = zeros(size(Rrho));
        Ks = zeros(size(Rrho));
        Rf = zeros(size(Rrho));
        
        sf = ((dtdz_adiab > 0) & (Rrho > 1) & (Rrho < 100));  % Salt fingering yes or no
        dl = ((dtdz_adiab < 0) & (Rrho < 1) & (Rrho > .01));  % Diffusive layering yes or no

        Ks(sf) = Kstar./( 1 + (Rrho(sf)/Rc).^n);
        Kt(sf) = 0.7*Kstar./(Rrho(sf).*(1+(Rrho(sf)/Rc).^n));

        C = 0.0032.*exp(4.8.*Rrho(dl).^.72);  % Unitless
        Ra = 0.25e9.*(Rrho(dl).^(-1.1));      % Unitless

        Rf(dl) = (1./Rrho(dl) + 1.4.*(1./Rrho(dl) - 1).^(1.5)) ...
              ./ (1 + 14*(1./Rrho(dl) - 1).^1.5); % Unitless

        Kt(dl) = C.*Ra.^(1/3).*lilkt;
        Ks(dl) = Rf(dl).*Rrho(dl).*Kt(dl);

        kt_dd = (dl | sf).*Kt;
        ks_dd = (dl | sf).*Ks;        
    end

    ks = k + ks_dd;
    kt = k + kt_dd;
    
    %% Implicit Scheme

    % I T(i+1) = T(i) + DT*T(i+1) + b
    % T(i+1) = (I - DT)\(T(i) + b)
    
    % Top = linear operator for flux through top of box, 
    % Top * T = (rho*Cp*deltat*ktop/dz^2)*(T(j+1) - T(j))
    % Bot = linear operator for flux through bottom of box,
    % Bot * T = (rho*Cp*deltat*kbot/dz^2)*(T(j) - T(j-1))
    % Fluxes through top of top box and bottom of bottom box are zero for
    % now
     CS = (deltat*ks)/dz^2;  % CFL number for salinity and temp
     CT = (deltat*kt)/dz^2;
     if (usesparse)
        STop = (spdiags([0;CS],1,Nz,Nz)-spdiags([CS;0],0,Nz,Nz));
        TTop = (spdiags([0;CT],1,Nz,Nz)-spdiags([CT;0],0,Nz,Nz));
        SBot = (spdiags([0;CS],0,Nz,Nz)-spdiags([CS;0],-1,Nz,Nz));
        TBot = (spdiags([0;CT],0,Nz,Nz)-spdiags([CT;0],-1,Nz,Nz));
     else
         STop = (diag(CS,1)-diag([CS;0]));
         TTop = (diag(CT,1)-diag([CT;0]));
         SBot = (diag([0; CS])-diag(CS,-1));
         TBot = (diag([0; CT])-diag(CT ,-1));
     end

    % Sparse

    
    DS = STop - SBot; % Diffusion operators
    DT = TTop - TBot; 

    % forcing / boundary terms
    % Salinity has a boundary term for flux at top and bottom
    bS = zeros(size(S));
    bS(1) = basal_dsdt*deltat;
    bS(end) = surface_dsdt*deltat;

    % Temperature has a forcing term to ensure that an unforced adiabatic
    % temperature profile has zero tendency.
    adiabatic_profile = [0;cumsum(dt_adiab)];
%    adiabatic_profile = adiabatic_profile - adiabatic_profile(end)+Tmelt;
    bT = -DT*adiabatic_profile;
    % Top and bottom prescribed heat flow
    % Heat flux from top
    bT(end) = bT(end) + surface_dtdt*deltat;
    % Flux from bottom
    bT(1) = bT(1) + basal_dtdt*deltat;
    
    % Calculate heat flux for diagnostics before updating temp
    q = -rhow*cpw*(TTop*T-TTop*adiabatic_profile)*dz/deltat;
    q(end) = lostbyocean;
    
    % Forward scheme, just to test
%    T = T + DT*T + bT;
%    S = S + DS*S + bS;
    if (usesparse)
        T = (speye(Nz,Nz)-DT)\(T+bT);
        S = (speye(Nz,Nz)-DS)\(S+bS);
    else
        T = (eye(Nz)-DT)\(T+bT);
        S = (eye(Nz)-DS)\(S+bS);
    end
% Sparse

%     S(S>2.7) = 2.7;
    
    Hseafloor=HS0-DHDT*2*HS0*sin(pi/2*t/t_hosc);%from Hussmann et al. 2004 eccentricity oscillations Fig. 8
       
    %% Averaging
    if (ia > max(averaging_time/deltat,1))
        ja = ja + 1;
        tavg(ja) = t;
        Savg(:,ja) = Saccum/ia;
        Tavg(:,ja) = Taccum/ia;
        cavg(:,ja) = caccum/ia;
        ktavg(:,ja) = ktaccum/ia;
        ksavg(:,ja) = ksaccum/ia;
        havg(:,ja) = haccum/ia;
        qavg(:,ja) = qaccum/ia;
        yearavg(:,ja) = t/3.16e7;
        Saccum = zeros(size(S));
        Taccum = zeros(size(T));
        qaccum = zeros(size(S));
        caccum = zeros(Nz-1,1);
        ktaccum = zeros(Nz-1,1);
        ksaccum = zeros(Nz-1,1);
        haccum = 0;
        fprintf('%06d %d %f\n',round(t/(86400*365.25)),round(tend/(86400*365.25)),Hseafloor);
        ia = 0;
        if isinf(Hseafloor)
            error('problem with Hseafloor ramping')
        end
    %% Save every 1000 avg steps
      if (~mod(ja,1000))
%         save(sprintf('checkpoint%d.mat',ja))                                    
        
        if exist('checkpoint_file')
              system(['rm ' checkpoint_file]);
        end
          japrev = ja;
  
        S0str = num2str(S0);S0str(S0str=='.')='p';
        Sfstr = num2str(Sfluxseafloor); Sfstr(Sfstr=='-')='m';
        Hstr = num2str(Hseafloor); Hstr(Hstr=='.')='p';
        Tstr = num2str(t/3.16e13); Tstr(Tstr=='.') = 'p';
        checkpoint_file = [str_EOS '_S0_' S0str '_T' Tstr 'Myr_H' Hstr 'Wm2_Sf' Sfstr '_' datestr(now,30) '.mat'];       
    	disp('Saving checkpoint');
        save(checkpoint_file)   

      end   
    else
        ia = ia + 1;
        Saccum = Saccum + S;
        Taccum = Taccum + T;
        caccum = caccum + convecting;
        ktaccum = ktaccum + Ks;
        ksaccum = ksaccum + Kt;
        haccum = haccum + h;
        qaccum = qaccum + q;
    end
end
    
