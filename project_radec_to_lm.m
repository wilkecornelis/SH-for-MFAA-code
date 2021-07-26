function [lmmap,radec] = project_radec_to_lm(ra, dec, radecmap, l, m, lon, lat, epoch, tobs)

%% lmmap = project_radec_to_lm(ra, dec, radecmap, l, m, lon, lat, epoch, tobs)
%
% This function maps a skymap in (ra, dec)-coordinates to a grid in the
% local horizon system specified in (l, m)-coordinates given the
% geographical location, epoch (B1950 or J2000) of the (ra, dec)-map and
% the time of observation. Coordinate conversions are done via the ITRF
% system to include the effect of precession.
%
% Arguments
% ra       : vector of values along the right ascension axis (in rad)
% dec      : vector of values along the declination axis (in rad)
% radecmap : (ra, dec)-map (matrix) of intensity values (in K)
% l        : vector of values along the l-axis in local horizon system
% m        : vector of values along the m-axis in local horizon system
% lon      : geographical longitude of observer (in rad)
% lat      : geographical latitude of observer (in rad)
% epoch    : Boolean, true if (ra, dec)-map is in B1950, false for J2000
% tobs     : time of observation (datenum)
%
% SJW, 21 February 2015

%% define (l,m)-grid
lmmap = zeros(length(m), length(l));
[lgrid, mgrid] = meshgrid(l, m);
dist = sqrt(lgrid.^2 + mgrid.^2);
lvec = lgrid(dist(:) < 0.99);
mvec = mgrid(dist(:) < 0.99);

% lvec = lgrid(:);
% mvec = mgrid(:);

%% convert vector of (l,m)-coordinates to unit vectors in ITRF
% The (l, m, n)-system is a right-handed system with l pointing east and n
% pointing towards zenith. The ITRF-system is a right-handed coordinate
% system with the first axis pointing towards (lon, lat) = (0, 0), the
% second axis pointing towards (lon, lat) = (pi/2, 0) (in rad) and the
% third axis towards (lon, lat) = (0, pi/2). If the n axis is aligned with
% the third axis of the ITRF-system, l is aligned with the *second* axis in
% the ITRF system while the m-axis is parallel to the first axis in the
% ITRF system but points in the opposite direction.
lmn0_ITRF = [-mvec, lvec, sqrt(1 - lvec.^2 - mvec.^2)];
rotmat_lat = [cos(pi/2 - lat),  0, sin(pi/2 - lat); ...
              0,                1, 0; ...
              -sin(pi/2 - lat), 0, cos(pi/2 - lat)];
rotmat_lon = [cos(-lon),  sin(-lon), 0; ...
              -sin(-lon), cos(-lon), 0; ...
              0,          0,         1];
lmn_ITRF = (rotmat_lon * rotmat_lat * lmn0_ITRF.').';
% size(lmn_ITRF)

%% convert pointings in ITRF to pointings in (ra, dec)
% calculate JulianDay
JD = JulianDay(tobs);
% Greenwich mean sidereal time in seconds
a1 = 24110.54841;
a2 = 8640184.812866;
a3 = 0.093104;
a4 = -6.2e-6;
polcoeff = [a4, a3, a2, a1];
TU = (floor(JD) + 0.5 - 2451545) / 36525;
GMST = (JD - floor(JD) - 0.5) * 86400 * 1.002737811906 + polyval(polcoeff, TU);
% conversion to radians
GMST = (GMST / 86400) * 2 * pi; % verified
% mod(GMST * 24 / (2 * pi), 24)
% rotation matrix to convert apparent positions to ITRF positions
rotmat = [ cos(GMST), sin(GMST), 0; ...
          -sin(GMST), cos(GMST), 0; ...
                   0,         0, 1];
% inverse (transpose) to go the other way around
app_radec_pos = rotmat.' * lmn_ITRF.';

% calculate precession w.r.t. J2000
precMat = precessionMatrix(JD);
% convert apparent positions to J2000 positions
% radecCart = precMat.' * app_radec_pos;
radecCart = app_radec_pos;
% size(radecCart)
% convert J2000 to B1950 if needed
% if (epoch)
%     precMat = precessionMatrix(JulianDay(datenum(1950, 1, 1, 0, 0, 0)));
%     radecCart = precMat.' * radecCart;
% end
% find (ra, dec)-positions to interpolate
[ravec, decvec, ~] = cart2sph(radecCart(1, :).', radecCart(2, :).', radecCart(3, :).');
% plot(ravec,decvec,'x')
%% map radecmap to lmmap using interpolation
lmmap(dist < 0.99) = interp2(ra, dec, radecmap, ravec, decvec, 'cubic');

radec = [ravec,decvec];
