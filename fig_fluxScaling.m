clear 
close all

%% Import data
load('data/flux_scaled.mat')

%% Plot data

plot(flux_out,'-o')
xlabel('Source index');
ylabel('Flux (Jy)')
legend('SUMSS','NVSS scaled')
grid on;
hold on;
xlim([0,50])
