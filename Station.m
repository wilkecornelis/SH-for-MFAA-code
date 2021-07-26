classdef Station < dynamicprops
    % Station class
    
    properties
        d                       % Inter-element spacing in m
        N                       % Number of elements in station
        M                       % Number of elements in a tile
        P                       % Number of receive paths (Also number of tiles)
        f                       % Current operating frequency (Hz)
        lambda                  % Current operating wavelength (m)
        Tsys                    % System temperature (Kelvin)
        tau                     % Integration time in seconds
        B                       % Integration bandwidth in Hz
        S_stat                  % Sensitivity of station (A/T)
        SEFD_stat               % SEFD of station
        SEFD_elem               % SEFD of individual element
        
        Dir_elem                % Directivity of a single element
        Ae_elem                 % Effective area of a single element (m^2)
        Ae_tile                 % Effective area of a tile (m^2)
        Ae_stat                 % Effective area of the station (m^2)
        
        HPBW_elem               % HPBW of an individual element
        HPBW_tile               % HPBW of a tile beam
        HPBW_stat               % HPBW of station beam
        FNBW_stat               % FNBW of station beam
        
        FoV_elem                % FoV of a single element
        FoV_tile                % FoV of a tile
        FoV_stat                % FoV of the station
        
        lat                     % station latitude in decimal degrees
        long                    % station longitude in decimal degrees
        normal                  % station normal vector
    end
    
    %% methods and functions
    methods
        % Constructor
        function obj = Station(d,N,M,Dir,Tsys,tau,B,lat,long)
            obj.d = d;
            obj.N = N;
            obj.M = M;
            obj.P = N/M;
            obj.Dir_elem = Dir;
            obj.Tsys = Tsys;
            obj.tau = tau;
            obj.B = B;
            obj.lat = lat;
            obj.long = long;
            obj.normal = [cosd(lat)*cosd(long) cosd(lat)*sind(long) sind(lat)];
        end
        
        function setF(obj,f_obs)
            c = 2.99792e8;
            kB = 1.38e-23;
            lambda = c/f_obs;
            
            obj.f = f_obs;
            obj.lambda = lambda;
            obj.Ae_elem = min(obj.lambda^2 * obj.Dir_elem / (4 * pi), obj.d^2);
            obj.Ae_tile = obj.M * obj.Ae_elem;
            obj.Ae_stat = obj.N * obj.Ae_elem;
            obj.HPBW_stat = 2*(pi/2 - acos(1.391*obj.lambda/(pi*sqrt(obj.N)*obj.d)));
            obj.FNBW_stat = 2*(pi/2 - acos(obj.lambda/(sqrt(obj.N)*obj.d)));
            
            if obj.M < 2
                obj.SEFD_stat = 1e26 * 2 * kB * obj.Tsys / obj.Ae_stat;
                obj.SEFD_elem = 1e26 * 2 * kB * obj.Tsys / obj.Ae_elem;
                obj.S_stat = obj.Ae_stat/obj.Tsys;
                
                obj.HPBW_elem =  acos(0.5);
                obj.FoV_elem = 0.25 * pi * obj.HPBW_elem^2;
            else
                obj.SEFD_stat = 1e26 * 2 * kB * obj.Tsys / sqrt(obj.Ae_tile * obj.Ae_stat);
                obj.SEFD_elem = 1e26 * 2 * kB * obj.Tsys / obj.Ae_elem;
                obj.S_stat = sqrt(obj.Ae_elem * obj.Ae_stat)/obj.Tsys;
                
                obj.HPBW_elem =  obj.lambda/(sqrt(obj.M)*obj.d);
                obj.FoV_elem = 0.25 * pi * obj.HPBW_elem^2;
            end
        end
        
        function changeM(obj,M)
            obj.M = M;
            obj.P = obj.N/obj.M;
        end
        
        function changeStatCoord(obj,lat,long)
            obj.lat = lat;
            obj.long = long;
        end
              
        function [lmn,bl,flux] = srcCat(obs,cat,srcsel,t_delta)
            % Evaluates a single source in provided catalogue 
            
            % inputs: 
            % obs      - observer object 
            % cat      - PSC catalogue 
            % srcsel   - index of source in PSC to evaluate
            % t_delta  - local time at which to evaluate
            
            % outputs:
            % lmn      - lmn coordinates of source
            % bl       - galactic coordinates of source
            % flux     - flux of source
            
            % Get RaDec of source
            ra = [cat(srcsel).alpha];
            dec = [cat(srcsel).delta];
                     
            % Calculate source position at time t_delta
            JD = datenum(t_delta) - datenum(1998, 2, 13, 12,0,0) + 2450858;
            srcpos_app(:,:,1) = radectolm(ra,dec,JD,obs.long,obs.lat).';
            
            flux(:,1) = [cat(srcsel).flux];
            bl(:,1) = [cat(srcsel).gal_lat];
            bl(:,2) = [cat(srcsel).gal_long];
            
            
            lmn = srcpos_app.';
        end
        
        function [lmn,sigmas] = genSkyBack(obs,fluxlim_low)
            % Place background sources randomly from source statistics
            % First, calculate station SEFD after integration to get lower
            % flux limit of source statistics
            
            % unused 

            fluxlim = logspace(log10(fluxlim_low),0,1000);                  % Flux limits of sources (max limited to the min of foreground flux)
            
            rhosrc = srcstat(fluxlim,obs.f);                                % Source density
            Nsrc = 2*pi*rhosrc;                                             % Source statistics for hemispherical FoV
            Nsrc_tot = floor(max(Nsrc))-ceil(min(Nsrc));                    % Total number of sources in FoV_elem
            % Use Archimedes Hat-Box theorem to project random cylindircal
            % coords to hemisphere. This ensures uniform distribution of
            % sources on hemisphere.
            phis = (2*pi*rand(Nsrc_tot,1));    % Random phi cylindrical coords
            % Calculate lmn coordinates for background sources (cylinder projected on hemisphere sphere)
            z = rand(Nsrc_tot,1);
            x = sqrt(1-z.^2).*cos(phis);
            y = sqrt(1-z.^2).*sin(phis);
            lmn = [x y z];
            
            % Calculate EEP towards background interferers
            thetas_back = acosd(z);
            EEP_back(:,1) = cosd(thetas_back);
            sigmas = (EEP_back.^2).*interp1(Nsrc,fluxlim,linspace(ceil(min(Nsrc)),Nsrc_tot,Nsrc_tot)).';  % Source powers
            sigmas = sigmas.';
         end
        
        function [lmn,bl,sigmas] = skyCat(obs,cat,srcsel,t_delta)
            % Place foreground sources from catalogue
            [lmn,bl,flux] = srcCat(obs,cat,srcsel,t_delta);
            thetas = acosd(lmn(3,:,:)./sqrt(lmn(1,:,:).^2 + lmn(2,:,:).^2 + lmn(3,:,:).^2));
            alt_app = asind(lmn(3,:));
            up = alt_app>0;
            lmn = lmn(:,up).';

            thetas = thetas(up);
            % Calculate apparent source powers (attenuated by EEP)  
            EEP = cosd((thetas));
%             sigmas = (EEP.^2).*flux(up).';
            sigmas = flux(up).';
            bl = bl(up,:);
        end
        
        function [SIR,SNR,int_count,sigma_cal,lm_intsel,flux_intsel] = calcSIRSNR(obs,lmn,sigmas,cal_sel)
            % Define antenna coordinates with reference to station
            r_stat = zeros(obs.N, 2);
            elem_ind = 1;
            for m = 0:sqrt(obs.N)-1
                for n = 0:sqrt(obs.N)-1
                    r_stat(elem_ind, :) = obs.d .* [m n];
                    elem_ind = elem_ind + 1;
                end
            end
            
            % Define antenna coordinates on a single tile
            r_tile = zeros(obs.M, 2);
            elem_ind = 1;
            for m = 0:sqrt(obs.M)-1
                for n = 0:sqrt(obs.M)-1
                    r_tile(elem_ind,:) = obs.d.*[m n];
                    elem_ind = elem_ind +1;
                end
            end
            
            % HPBW radius
            r_HPBW = sin(obs.HPBW_stat/2);
            r_FNBW = sin(obs.FNBW_stat/2);
            
            % Initial calibration source flux
            cal_flux(1) = sigmas(cal_sel);
              
            disp(['Cal source alt: ', num2str(acosd(sqrt(lmn(cal_sel,1).^2 + lmn(cal_sel,2).^2)))]);
            disp(['Initial cal source flux (apparent): ', num2str(cal_flux(1))]);
            disp(['Total source count is: ', num2str(length(lmn))]);
            
            % Geometric delay vector of cal source for each tile
            a_cal_tile = exp(-1j*(2*pi/obs.lambda)*r_tile*lmn(cal_sel,1:2).');   % delay vector
            a_cal_stat = exp(-1j*(2*pi/obs.lambda)*r_stat*lmn(cal_sel,1:2).');   % delay vector
            
            % Loop over interferers and determine denominator of SIR
            % expression
            Nsrc_tot = length(lmn);
            den = 0;
            sigma_app = zeros(length(Nsrc_tot));
            int_count = 0;
            cal_count = 2;
            for int_idx = 1:Nsrc_tot
                % Calculate tile beam towards current source (pointed at brightest source)
                if obs.M>1
                    AF_tile = abs(sum(exp(-1i*(2*pi/obs.lambda)*r_tile*lmn(int_idx,1:2).').*(conj(a_cal_tile))));
                else
                    AF_tile = obs.M;
                end
                % phase delay vector of source after beamforming to cal
                % source
                a = exp(-1i*(2*pi/obs.lambda)*(r_stat*lmn(int_idx,1:2).')).*conj(a_cal_stat);

                % Apparent source power after tile beamforming
                sigma_app(int_idx) = (AF_tile/obs.M)^2.*sigmas(int_idx);
                
                % Distance between current source and cal source
                dist = sqrt((lmn(cal_sel,1)-lmn(int_idx,1))^2 + (lmn(cal_sel,2)-lmn(int_idx,2))^2);
                
                % Denominator SIR expression
                if sigma_app(int_idx)>(obs.SEFD_stat/sqrt(obs.B*obs.tau))
                    if int_idx ~= cal_sel && dist >r_FNBW
                        int_count = int_count+1;
                        arr_at = (a*(a'*ones(obs.N,1)));
                        den = den + (sigma_app(int_idx)*(arr_at-ones(obs.N,1)));
          
                        lm_intsel(int_count,:) = lmn(int_idx,1:2);
                        flux_intsel(int_count) = sigma_app(int_idx)*(ones(1,obs.N)*arr_at);
                        
                    % If source is inside HPBW then it is added to the
                    % calibration source
                    elseif int_idx ~= cal_sel
                        % Apparent power of calibration source taking into
                        % account attenuation by station beam
                        cal_flux(cal_count) = sigmas(int_idx)*abs(sum(a*(a'*ones(obs.N,1))))/(obs.N)^2;
                        cal_count = cal_count+1;
                    end
                end    
            end

            % Apparent power of calibration source
            sigma_cal = sum(cal_flux);
            disp(['Final cal source flux (apparent): ', num2str(sigma_cal)]);

            % Numerator of SIR expression
            num = (obs.M^2)*sigma_cal*(obs.P-1)*ones(obs.P,1);
            % Calculate SIR
            SIR = (mean(num)/mean(abs(den)));
            % Calculate SNR
            SNR = ((sigma_cal)/obs.SEFD_elem);
            
            disp(['Calculated SIR: ', num2str(SIR)]);
        end
    end
end

