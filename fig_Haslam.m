clear
close all

load('data\Haslam_orig.mat');

imagesc((rad2deg(ra)),rad2deg(dec),log10(intensity408MHz/10));
colormap('jet')
set(gca,'YDir', 'normal', 'XDir', 'normal');

xlabel('\alpha (deg)')
ylabel('\delta (deg)')
daspect([1 1 1])

