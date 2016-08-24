% %  plotlist = {'T','S','rho','nsquared','q','h','v','dz0dz'};
%  plotlist = {'T','S','rho','q','h','v','dz0dz'};
% % plotlist = {'T','c','nsquared','h'}
n = length(plotlist);
switch n
    case 1,
        figr = 1; figc = 1;
    case 2,
        figr = 2; figc = 1;
    case 3
        figr = 3; figc = 1;
    case 4
        figr = 4; figc = 1;
    case {5,6}
        figr = 3; figc = 2;
    case {7,8}
        figr = 4; figc = 2;
end

if strcmp(str_EOS,'mgso4')
    str_Sunits = 'molal';
elseif strcmp(str_EOS,'gsw302') | strcmp(str_EOS,'gsw305')
     str_Sunits = 'SA';
end

% tstart = 
[ro,co] = size(Tavg);

% tp = (1:co).*tavg/(1e6*3.16e7); % time in Myr
tp = tavg/(1e6*3.16e7); % time in Myr; s. vance revised this on June 28 2016

z = ((ro:-1:1)-.5).*dz/1000; % 
pp = repmat(p,[1 co]);

clf



for i=1:n
    switch plotlist{i}
        case 'T'
            subplot(figr,figc,i,'align')
            switch surf_type
                case 'imgsc'
                     imagesc(tp,z,Tavg);
                case 'sf'
                    hs = surf(tp,z,Tavg);shading interp;view(0,90); axis tight;set(gca,'ydir','reverse');
            end
            title('Temperature (K)');
        case 'S'
            subplot(figr,figc,i,'align')
            switch surf_type
                case 'imgsc'
                    imagesc(tp,z,Savg);
                case 'sf'
                    surf(tp,z,Savg);shading interp;view(0,90); axis tight;set(gca,'ydir','reverse');
            end
            title(['Salinity (' str_Sunits ')']);
        case 'c'
            subplot(figr,figc,i,'align')
            switch surf_type
                case 'imgsc'
                    imagesc(tp,z,cavg);
                case 'sf'
                    surf(tp,z(2:end),cavg);shading interp;view(0,90); axis tight;set(gca,'ydir','reverse');
            end
            title('Convective activity');
        case 'nsquared'
%             de = mgso4_dens(Savg(1:end-1,:),Tavg(1:end-1,:),pp(1:end-1,:));
%             pde = mgso4_pden(Savg(2:end,:),Tavg(2:end,:),pp(2:end,:),pp(1:end-1,:));
            de = swEOS.dens(Savg(1:end-1,:),Tavg(1:end-1,:),pp(1:end-1,:));
            [np,nt]=size(Tavg);
            switch str_EOS 
                case 'gsw302' 
                for ip = 1:np-1
                    parfor it = 1:nt
                        pde(ip,it) = swEOS.pden(Savg(ip+1,it),Tavg(ip+1,it),pp(ip+1,it),pp(ip,it));
                    end
                end
                case 'mgso4'
                pde = mgso4_pden(Savg(2:end,:),Tavg(2:end,:),pp(2:end,:),pp(1:end-1,:));
            end
  
            nsquared = g.*(de-pde)./(dz.*de);
            
            subplot(figr,figc,i,'align')
            switch surf_type
                case 'imgsc'
                    imagesc(tp,z,nsquared);
                case 'sf'
                    imagesc(tp,z,nsquared);
            end
            %             surf(t,z,nsquared);shading interp;view(0,90); axis tight
            title('Stratification (N^2), s^{-2}');
        case 'rho'
            de = swEOS.dens(Savg(1:end-1,:),Tavg(1:end-1,:),pp(1:end-1,:));
            
            subplot(figr,figc,i,'align')
            surf(tp,z(1:end-1),de);shading interp;view(0,90); axis tight;set(gca,'ydir','reverse');
            title('Density (kg m^{-3})')
        case 'v'
%             vel = mgso4_vel(Savg(1:end-1,:),Tavg(1:end-1,:),pp(1:end-1,:));
            vel = swEOS.vel(Savg(1:end-1,:),Tavg(1:end-1,:),pp(1:end-1,:));
            
            subplot(figr,figc,i,'align')
%             imagesc(t,z,vel);
            surf(tp,z(1:end-1),vel);shading interp;view(0,90); axis tight;set(gca,'ydir','reverse');
            title('Sound Speed (kms s^{-1})')
        case 'dz0dz'
%             vel = mgso4_vel(Savg(1:end-1,:),Tavg(1:end-1,:),pp(1:end-1,:));
            vel = swEOS.vel(Savg(1:end-1,:),Tavg(1:end-1,:),pp(1:end-1,:));
            if ~exist('de')
                de = swEOS.dens(Savg(1:end-1,:),Tavg(1:end-1,:),pp(1:end-1,:));
            end
            if ~exist('vel')
                vel = swEOS.vel(Savg(1:end-1,:),Tavg(1:end-1,:),pp(1:end-1,:));
            end
            delz = 1e3*gradient(z(1:end-1))'*ones(1,length(tp));
            z0 = de.*vel*1000;
            dz0dz = gradient(z0)./delz;
%             subplot(figr,figc,i,'align')
%             toffset = 1000;
            imagesc(tp,z,dz0dz);
            title('Vertical Acoustic Impedance Gradient (kg m^{-3} s^{-1})')
        case 'q'
            subplot(figr,figc,i,'align')
%             imagesc(t,z,qavg,[-.4 .1]);
            surf(tp,z,qavg);shading interp;view(0,90); axis tight;set(gca,'ydir','reverse');
            title('Heat flux (W/m^2)');
        case 'h'
            subplot(figr,figc,i,'align')
            plot(tp,havg/1e3);
            ylabel('Thickness (km)');
            title('Ice thickness (km)');
        case 'kt'
            subplot(figr,figc,i,'align')
%             imagesc(t,z,ktavg);
            surf(tp,z(2:end),ktavg);shading interp;view(0,90); axis tight;set(gca,'ydir','reverse');
            title('Double diffusive mixing K_T (m^2 s^{-1})');
        case 'ks'
            subplot(figr,figc,i,'align')
%             imagesc(t,z,ksavg);
            surf(tp,z(2:end),ksavg);shading interp;view(0,90); axis tight;set(gca,'ydir','reverse');
            title('Double diffusive mixing K_S (m^2 s^{-1})');
    end
%    set(gca,'XLim',[0 600])
    if (i > n - figc)
        xlabel('Time (Myr)');
    end
    if ((mod(i,figc)==1 || figc == 1) && ~strcmp(plotlist{i},'h'))
        ylabel('Depth (km)');
    end
    if plotlist{i}~='h'
        colorbar
    end
end
fh = gcf;
fh.Name = ['S0:' num2str(S0) ' ' str_Sunits ' Hsf:' num2str(Hseafloor) ' Ssf:' num2str(Sfluxseafloor) ' dt:' num2str(deltat/3.16e7) ' yr ' ...
' dz:' num2str(dz/1000) ' km'];