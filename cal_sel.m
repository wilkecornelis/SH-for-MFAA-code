function cal_sel = cal_sel(lmn,bl,flux,r_FNBW)
% This function selects an optimal calibration source from the provided PS
% list

% Identify candidate calibration sources as sources outside GP (10 degrees) and above
% 30 degrees from horizon
alt_app = asind(lmn(:,3));
calCand_sel = (bl(:,1)<deg2rad(-10)|bl(:,1)>deg2rad(10)) & alt_app>20;
lmn_calCand = lmn(calCand_sel,1:2);
flux_calCand = flux(calCand_sel,1);
[fluxCalCand_sorted,sortorder] = sort(flux_calCand,'descend');

% Distance between sources
l = meshgrid(lmn_calCand(sortorder, 1)) - meshgrid(lmn_calCand(sortorder, 1)).';
m = meshgrid(lmn_calCand(sortorder, 2)) - meshgrid(lmn_calCand(sortorder, 2)).';
lmDist = sqrt(l.^2 + m.^2);

srcSel = false(length(lmDist),length(lmDist));
compoundFlux = zeros(length(lmDist),1);

for srcIdx = 1:length(lmDist)
    % cluster sources (sources within FNBW) for this potential calibrator
    srcSel(srcIdx,:) = lmDist(srcIdx,:)<r_FNBW; 
    % compund flux for this potential calibrator
    compoundFlux(srcIdx,1) = sum(fluxCalCand_sorted(srcSel(srcIdx,:)));
end

ranked =  sort(compoundFlux,'descend');
bestCalSrc_ind = find(compoundFlux==ranked(1));

% index of optimal calibration source in provided PS list
cal_sel = find(flux==fluxCalCand_sorted(bestCalSrc_ind(1)));
    
end

