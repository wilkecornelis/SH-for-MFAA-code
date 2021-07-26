clear 
close all

%% Import results

for t_ind = 1:144
    load(['results/SIRvT/M1/SIR_t',num2str(t_ind),'.mat'])
    SIR_M1(t_ind) = 10*log10(SIR);
    sigma_M1(t_ind) = sigma_cal;
    t(t_ind) = t_obs;
 
    load(['results/SIRvT/M4/SIR_t',num2str(t_ind),'.mat'])
    SIR_M4(t_ind) = 10*log10(SIR);
    sigma_M4(t_ind) = sigma_cal;

    
    load(['results/SIRvT/M16/SIR_t',num2str(t_ind),'.mat'])
    SIR_M16(t_ind) = 10*log10(SIR);
    sigma_M16(t_ind) = sigma_cal;

end

figure
plot(t,SIR_M1,'-')
hold on
plot(t,SIR_M4,'-')
plot(t,SIR_M16,'-')
yyaxis right
plot(t,sigma_M1,'-')
ylabel('apparent flux (Jy)')
yyaxis left
ylabel('SIR (dB)')
xlabel('time')
grid on
pbaspect([2 1 1])
legend('M=1','M=2','M=4','\sigma_c','Location','NorthOutside')




