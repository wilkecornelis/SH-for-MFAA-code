clear 
close all

% This simulation calculates the SIR over 24 hours at 10 minute intervals for different MANTIS configs

%% Set defaults
c = 2.99792e8;              % speed of light in m/s
f_0 = 1e9;                  % center frequency in Hz
lambda_0 = c/f_0 ;          % Wavelength in m

d = lambda_0/2;             % Inter-element spacing in m
N = 336^2;                  % Number of elements in station

Dir = 4;                    % Directivity of a single receiveing element
Tsys = 39;                  % System temperature in K
tau_0 = 1;                  % Integration time in seconds
B = 10e6;                   % integration bandwidth in Hz

%% Setup sim

% SKA site Karoo
lat = -30.72;
long = 21.41;

% Define 3 MANTIS station objects with different configurations

% MANTIS with 1x1 tiles
M = 1;
MANTIS_M1 = Station(d,N,M,Dir,Tsys,tau_0,B,lat,long);
setF(MANTIS_M1,0.843e9);

% MANTIS with 2x2 tiles
M = 4;
MANTIS_M4 = Station(d,N,M,Dir,Tsys,tau_0,B,lat,long);
setF(MANTIS_M4,0.843e9);

% MANTIS with 4x4 tiles
M = 16;
MANTIS_M16 = Station(d,N,M,Dir,Tsys,tau_0,B,lat,long);
setF(MANTIS_M16,0.843e9);

% upper and lower declination limits for Haslam (degrees)
decH = 60;
decL = -90;

% upper and lower limits galactic latitude limits for Haslam region in sky model (degrees)
bH = 10;
bL = -10;

%% Loop over 24 hours to get an idea of the calibration source selection
t_obs = datetime(now,'ConvertFrom','datenum');
t_obs.Hour = t_obs.Hour-t_obs.Hour;
t_obs.Minute = t_obs.Minute-t_obs.Minute;
t_obs.Second = t_obs.Second-t_obs.Second;

% first null beamwidth of station beam
r_FNBW = sin(MANTIS_M1.FNBW_stat/2);

% Import point source catalogue (this is a pre-combined NVSS and SUMSS catalogue)
load('data/srclistSHem.mat');
SHem = srclistSHem;

% Loop over 24 hours
for t_ind = 1:144
    
    disp(['working on time index ', num2str(t_ind), ' of 144'])
    disp(' ')
    
    % set time
    t_obs.Minute = t_obs.Minute + 10;
    
    % Identify calibration source from point source catalogues 
    % first filter out sources with flux lower than 1 Jy
    srcsel_SHem = [SHem.flux]>1;
    [lmn_SHem,bl,flux_SHem] = srcCat(MANTIS_M1,SHem,srcsel_SHem,t_obs);
    
    % Apply EEP attenuation to get apparent powers
    alt_app = asind(lmn_SHem(:,3));
    up = alt_app>0&(bl(:,1)<deg2rad(-10)|bl(:,1)>deg2rad(10));
    flux_SHem = sind(alt_app(up)).^2.*flux_SHem(up,:);
    % lmn of relevant sources
    lmn_SHem = lmn_SHem(up,:);
    % galactic coordinates of relevant sources
    bl = bl(up,:);
    
    % Call cal identify function to find flux of optimal calibration source
    cal_ind = cal_sel(lmn_SHem,bl,flux_SHem,r_FNBW);
    cal_flux = flux_SHem(cal_ind);
    
    % Import all sources from point source catalogue with flux > 0.4 Jy
    srcsel_SHem = [SHem.flux]>0.4;
    [lmn_SHem,bl,flux_SHem] = srcCat(MANTIS_M1,SHem,srcsel_SHem,t_obs);
    alt_app = asind(lmn_SHem(:,3));
    up = alt_app>0&(bl(:,1)<deg2rad(-10)|bl(:,1)>deg2rad(10));
    flux_SHem = sind(alt_app(up)).^2.*flux_SHem(up,:);
    lmn_SHem = lmn_SHem(up,:);
     
    % Import sourcelist for GP from Haslam
    [lm_Haslam,flux_Haslam] = GP_Haslam(MANTIS_M1,decH,decL,bH,bL,datenum(t_obs));
    % Select all sources
    srcsel_Haslam = flux_Haslam>0;
    lm_Haslam = lm_Haslam(srcsel_Haslam,:);
    % Calculate apparent altitude of all sources
    alt_app = acosd(sqrt(lm_Haslam(:,1).^2 + lm_Haslam(:,2).^2));
    % Apply EEP attenuation and subtract background flux
    flux_Haslam = sind(alt_app).^2.*(flux_Haslam(srcsel_Haslam,1)-min(flux_Haslam(srcsel_Haslam,1)));
    
   % Combine all sources into single source list
   lm = [lm_Haslam;lmn_SHem(:,1:2)];
   flux = [flux_Haslam;flux_SHem];
   
   % Find calibration source index in combined list
   cal_ind = find(flux == cal_flux); % might need a better way to do this... 
   
   % Calculate SIR's by calling calcSIRSNR function in Station class
   
   % Calculate SIR for M = 1x1
   disp('MANTIS M1 SIR calculation outputs: ')
   disp(' ')
   [SIR,SNR,int_Nsel,sigma_cal,lm_intsel,flux_intsel] = calcSIRSNR(MANTIS_M1,lm,flux,cal_ind);
   disp(' ') 
   
%    save(['results/SIRvT/M1/SIR_t' num2str(t_ind) '.mat'], ...
%        't_obs','lm','lm_intsel','flux','flux_intsel','sigma_cal','cal_ind','SIR','SNR');
%    
   clear SIR SNR int_Nsel sigma_cal lm_intsel flux_intsel
   
   % Calculate SIR for M = 2x2
   disp('MANTIS M4 SIR calculation outputs: ')
   disp(' ')  
   [SIR,SNR,int_Nsel,sigma_cal,lm_intsel,flux_intsel] = calcSIRSNR(MANTIS_M4,lm,flux,cal_ind);
   disp(' ')  
%    save(['results/SIRvT/M4/SIR_t' num2str(t_ind) '.mat'], ...
%         't_obs','lm','lm_intsel','flux','flux_intsel','sigma_cal','cal_ind','SIR','SNR');  
%     
    clear SIR SNR int_Nsel sigma_cal lm_intsel flux_intsel
    
   % Calculate SIR for M = 4x4 
   disp('MANTIS M16 SIR calculation outputs: ')
   disp(' ')     
   [SIR,SNR,int_Nsel,sigma_cal,lm_intsel,flux_intsel] = calcSIRSNR(MANTIS_M16,lm,flux,cal_ind);
   disp(' ') 
%    save(['results/SIRvT/M16/SIR_t' num2str(t_ind) '.mat'], ...
%         't_obs','lm','lm_intsel','flux','flux_intsel','sigma_cal','cal_ind','SIR','SNR');

end    





