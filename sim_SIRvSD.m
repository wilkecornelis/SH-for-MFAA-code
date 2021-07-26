clear
close all

%% Set defaults
c = 2.99792e8;              % speed of light in m/s
f_0 = 845e6;                % Operating frequency in Hz
lambda_0 = c/f_0 ;          % Wavelength in m

d = lambda_0/2;             % Inter-element spacing in m
N = 336^2;                  % Number of elements in station

Dir = 4;                    % Directivity of a single receiveing element
Tsys = 39;                  % System temperature in K
tau_0 = 1;                  % Integration time in seconds
B = 10e6;                   % integration bandwidth in Hz

%% Start sim

% Stellenbosch
lat = -33.982;
long = 18.910;

% MANTIS with 1x1 tiles
M = 1;
MANTIS_M1 = Station(d,N,M,Dir,Tsys,tau_0,B,lat,long);
setF(MANTIS_M1,0.845e9);

% upper and lower declination limits for Haslam (degrees)
decH = 40;
decL = -90;
% upper and lower limits galactic latitude limits for Haslam (degrees)
bH = 10;
bL = -10;

% lower flux limits to investigate
flux_lims = logspace(-1.3,0,30);

%% 
% load local time used in this analysis
load('data/t_SD.mat');

r_HPBW = sin(MANTIS_M1.HPBW_stat/2);
r_FNBW = sin(MANTIS_M1.FNBW_stat/2);

% Import point source catalogue (SUMSS and NVSS combined)
load('data/srclistSHem.mat');
SHem = srclistSHem;

% Identify calibration source from PSC
srcsel_SHem = [SHem.flux]>1;
[lmn_SHem,bl,flux_SHem] = srcCat(MANTIS_M1,SHem,srcsel_SHem,t_obs);

% Apply EEP attenuation to get apparent powers
alt_app = asind(lmn_SHem(:,3));
up = alt_app>0&(bl(:,1)<deg2rad(-10)|bl(:,1)>deg2rad(10));
flux_SHem = sind(alt_app(up)).*flux_SHem(up,:);
lmn_SHem = lmn_SHem(up,:);
bl = bl(up,:);

% Call cal identify function to find flux of optimal calibration source
calIdx = cal_sel(lmn_SHem,bl,flux_SHem,r_FNBW);
cal_flux = flux_SHem(calIdx);

% loop over lower flux limits
for SD_ind = 1:length(flux_lims)
       
    % Import all sources from point source catalogue above given lower flux
    % limit
    srcsel_SHem = [SHem.flux]>flux_lims(SD_ind);
    [lmn_SHem,bl,flux_SHem] = srcCat(MANTIS_M1,SHem,srcsel_SHem,t_obs);
    alt_app = asind(lmn_SHem(:,3));
    up = alt_app>0&(bl(:,1)<deg2rad(-10)|bl(:,1)>deg2rad(10));
    flux_SHem = sind(alt_app(up)).*flux_SHem(up,:);
    lmn_SHem = lmn_SHem(up,:);

    % Import sourcelist for GP from Haslam
    [lm_Haslam,flux_Haslam] = GP_Haslam(MANTIS_M1,decH,decL,bH,bL,datenum(t_obs));
    % Select all sources
    srcsel_Haslam = flux_Haslam>0;
    lm_Haslam = lm_Haslam(srcsel_Haslam,:);
    % Calculate apparent altitude of all sources
    alt_app = acosd(sqrt(lm_Haslam(:,1).^2 + lm_Haslam(:,2).^2));
    % Apply EEP attenuation and subtract background flux
    flux_Haslam = sind(alt_app).*(flux_Haslam(srcsel_Haslam,1)-min(flux_Haslam(srcsel_Haslam,1)));

    % Combine all sources into single source list (GP included)
    lm_wGP = [lm_Haslam;lmn_SHem(:,1:2)];
    flux_wGP = [flux_Haslam;flux_SHem];
    
    % Sourcelist for point sources (PSs) only (GP excluded)
    lm_woGP = lmn_SHem(:,1:2);
    flux_woGP = flux_SHem;
    
    int_count = length(flux_SHem);

    % Find calibration source index in combined list
    calIdx_wGP = find(flux_wGP == cal_flux);
    
    % Find calibration source index in PS only list
    calIdx_woGP = find(flux_woGP == cal_flux);    
  
    % Calculate SIR for combined source list (GP included)
    [SIR,SNR,int_Nsel,sigma_cal,lm_intsel,flux_intsel] = calcSIRSNR(MANTIS_M1,lm_wGP,flux_wGP,calIdx_wGP);
    
%     save(['results/results_SIRvSD/wGP',num2str(SD_ind),'.mat'], ...
%         't_obs','sigma_cal','calIdx','SIR','SNR','int_count','int_Nsel','lm_intsel','flux_intsel');

    % Calculate SIR for PS only list (GP excluded)
    [SIR,SNR,int_Nsel,sigma_cal,lm_intsel,flux_intsel] = calcSIRSNR(MANTIS_M1,lm_woGP,flux_woGP,calIdx_woGP);
    
%     save(['results/results_SIRvSD/woGP',num2str(SD_ind),'.mat'], ...
%         't_obs','sigma_cal','calIdx','SIR','SNR','int_count','int_Nsel','lm_intsel','flux_intsel');
%     
end


   




