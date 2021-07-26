function [b,l] = eq2gal(alpha,delta)
%EQTOGAL Summary of this function goes here
%   Detailed explanation goes here
    
    % Galactic coordinate constants
    delta_0 = 27.13;
    alpha_0 = 192.86;
    l_0 = 122.93;
   
    b = asin(sin(delta.').*sin(deg2rad(delta_0)) + cos(delta.').*cos(deg2rad(delta_0)).*(cos(alpha.'-deg2rad(alpha_0))));
%     l = -1*(atan2((cos(delta.')*sin(ra_list(idx)-deg2rad(alpha_0))),(sin(dec_list(idx)).*cosd(delta_0)-cos(dec_list(idx)).*sind(delta_0).*cos(ra_list(idx)-deg2rad(alpha_0))))-deg2rad(l_0));

    l = -1*(atan2(cos(delta.').*sin(alpha.'-deg2rad(alpha_0)),sin(delta.').*cosd(delta_0)-cos(delta.').*sind(delta_0).*cos(alpha.'-deg2rad(alpha_0)))-deg2rad(l_0));
end

