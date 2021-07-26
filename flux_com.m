clear all
close all

%% Set defaults
c = 2.99792e8;              % speed of light in m/s
f_0 = 845e6;                  % Operating frequency in Hz
lambda_0 = c/f_0 ;          % Wavelength in m

d = lambda_0/2;             % Inter-element spacing in m
N = 336^2;;                 % Number of elements in station

Dir = 4;                    % Directivity of a single receiveing element
Tsys = 39;                  % System temperature in K
tau_0 = 1;                  % Integration time in seconds
B = 10e6;                   % integration bandwidth in Hz

%% Start sim

% Stellenbosch
lat = -33.982;
long = 18.910;

% MANTIS with 1x1 tiles
M = 1;;
MANTIS_M1 = Station(d,N,M,Dir,Tsys,tau_0,B,lat,long);
setF(MANTIS_M1,0.845e9);


%% calculate SIR

% Current local datetime
% t_curr = datetime(now,'ConvertFrom','datenum');
% Current UTC datetime
% t_curr.Hour = t_curr.Hour+3;

% save('t_curr_GPcomp.mat','t_curr');
% 
load('t_curr_GPcomp.mat');
t_curr.Minute = t_curr.Minute;


% Get sources from Haslam for -10<b<10 and delta<-30
% for shift = 1:1
%     
%     dec_limitH = 40-shift;
%     dec_limitL = -10-shift;
% 
%     [lm_Haslam,flux_Haslam] = GP_Haslam(MANTIS_M1,dec_limitH,dec_limitL,datenum(t_curr));
% 
% 
%     % Get sources from point source cats
%     load('srclistSHem_full.mat');
%     SHem = srclistSHem;
% 
% 
%     reg_sel1 = [SHem.gal_lat]>deg2rad(-10) & [SHem.gal_lat]<deg2rad(10);
%     reg_sel2 = [SHem.delta]<deg2rad(dec_limitH)& [SHem.delta]>deg2rad(dec_limitL);
%     reg_sel = reg_sel1&reg_sel2;
%     % 
%     [lmn_SHem,bl,flux_SHem] = srcCat(MANTIS_M1,SHem,reg_sel,t_curr);
%     alt_app = asind(lmn_SHem(:,3));
%     up = alt_app>0;
% 
%     % % Select sources above horizon and above flux limit;
%     lmn_SHem = lmn_SHem(up,:);
%     flux_SHem = flux_SHem(up,:);
% 
%     flux_rat(shift) = (sum(flux_Haslam))/(sum(flux_SHem))
% end

% upper and lower declination limits for Haslam (degrees)
decH = 30;
decL = -90;

% upper and lower limits galactic latitude limits for Haslam (degrees)
bH = 10;
bL = -10;

[lm_Haslam,flux_Haslam] = GP_Haslam(MANTIS_M1,decH,decL,bH,bL,datenum(t_curr));
plotSky(MANTIS_M1,lm_Haslam,flux_Haslam,'x',1)

%% Get sources from point source cats
load('t_curr_GPcomp.mat');
load('srclistSHem.mat');
SHem = srclistSHem;
% 
reg_sel1 = [SHem.gal_lat]>deg2rad(bL) & [SHem.gal_lat]<deg2rad(bH);
reg_sel2 = [SHem.delta]<deg2rad(decH)& [SHem.delta]>deg2rad(decL);
reg_sel = reg_sel1&reg_sel2;

flux_sel = [SHem.flux]>0;

srcsel_SHem = reg_sel&flux_sel;

% load('t_curr.mat');
[lmn_SHem,bl,flux_SHem] = srcCat(MANTIS_M1,SHem,srcsel_SHem,t_curr);
alt_app = asind(lmn_SHem(:,3));
up = alt_app>0;
% 
% % % Select sources above horizon and above flux limit;
lmn_SHem = lmn_SHem(up,:);
flux_SHem = flux_SHem(up,:);
% sum(flux_SHem)
% plotSky(MANTIS_M1,lm_Haslam,flux_Haslam,'x',1)
% 
% plot(lmn_SHem(:,1),lmn_SHem(:,2),'rx');
% text(lmn_SHem(:,1),lmn_SHem(:,2),num2str(flux_SHem(:,1)));
% set(gca,'YDir', 'normal', 'XDir', 'reverse');

sum(flux_Haslam)/sum(flux_SHem)

