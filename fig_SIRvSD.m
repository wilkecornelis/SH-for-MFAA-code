clear 
close all

%% Import results
flux_lims = logspace(-1.3,0,30);

for SD_ind = 1:30
    % Results without galactic plane
    load(['results/SIRvSD/woGP/results_SIRvSD',num2str(SD_ind),'.mat']);
    SIRwoGP(SD_ind) = SIR;
    % number of interfererers
    Nint(SD_ind) = int_count;
    
    % Results with galactic plane
    load(['results/SIRvSD/wGP/results_SIRvSD',num2str(SD_ind),'.mat']);
    SIRwGP(SD_ind) = SIR;
    
end
% 
figure
plot(flux_lims,10*log10(SIRwoGP))
hold on
plot(flux_lims,10*log10(SIRwGP))

grid on
xlabel('lower flux lim (Jy)')
ylabel('SIR (dB)')
legend('GP excluded','GP included')


figure
plot(flux_lims,Nint)
xlabel('lower flux lim (Jy)')
ylabel('number of sources')

grid on

