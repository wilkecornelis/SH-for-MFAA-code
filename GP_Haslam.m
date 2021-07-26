function [lm_sel,flux_sel] = GP_Haslam(obs,decH,decL,bH,bL,t_obs)
% this function returns a list of pixels (from the Haslam map) in the galactic plane and their associated flux
% densities. 

% inputs:

% obs           - observer object
% decH          - upper declination angle limit
% decL          - lower declination angle limit
% bH            - upper galactic latitude limit
% bL            - lower galacitc latitude limit

% outputs 

% lm_sel        - lm coordinates of selected pixels
% flux_sel       - fluxes of selected pixels 

% Constants
c = 2.99792e8;              % speed of light in m/s
kB = 1.38e-23;              % Boltzman constant

% Load original Haslam map
load data/Haslam_orig.mat

% Sample map at N*tims array resolution
% length_side = obs.d*sqrt(obs.P);
res = sin(obs.HPBW_stat);
dl = res;

% project data on (l,m)-grid
l = 0:dl:1;
l = [-fliplr(l(2:end)), l];
[lgrid, mgrid] = meshgrid(l);
[lmmap,radec] = project_radec_to_lm(ra, dec, intensity408MHz, l, l, deg2rad(obs.long), deg2rad(obs.lat), 1, t_obs);

% define list of (l,m)-pixels removing pixels too close to horizon
dist = sqrt(lgrid.^2 + mgrid.^2);
lm = [lgrid(dist(:) < 0.99), mgrid(dist(:) < 0.99)];

% apply appropriate physical conversion and corrections
% convert map to units of Kelvin (at 408 MHz)
Tmap = lmmap / 10;

% convert from temperature map to intensity map (in Jy / sr) (see TMS)
lambda408 = c / 408e6;
Imap = 1e26 * (2 * kB * Tmap / lambda408^2);
% convert to flux in Jy / pixel
dOmega = dl.^2 ./ sqrt(1 - dist.^2);
dOmega(dist >= 0.99) = NaN; % remove pixels that can potentially blow up
Smap = Imap .* dOmega;
S408MHz = Smap(dist(:) < 0.99);
lambda = c/843e6;
S = S408MHz * (lambda / lambda408).^0.55; % correction for sky noise spectrum

% Filter region from Haslam flux map
[b,gal_lon] = eq2gal(radec(:,1),radec(:,2));
% Select gal latitudes between bL and bH degrees
reg_sel1 = b>deg2rad(bL) & b<deg2rad(bH);
% Select declinations between decH and decL
reg_sel2 = radec(:,2)<deg2rad(decH)& radec(:,2)>deg2rad(decL);
% Combine conditions
reg_sel = reg_sel1&reg_sel2.';

% Return filtered GP region
lm_sel = [lm(reg_sel,1),lm(reg_sel,2)];
flux_sel = S(reg_sel,1);

end

