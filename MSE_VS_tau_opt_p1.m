% This MATLAB function was developed to generate numerical results for:
%
% Italo Atzeni and Antti Tölli, "Channel Estimation and Data Detection
%     Analysis of Massive MIMO with 1-Bit ADCs," IEEE Trans. Wireless
%     Commun. (to appear), 2021.
% -------------------------------------------------------------------------
% Description: Fig. 3 --> Plot the optimal MSE of the channel estimation
%     and the optimal transmit SNR against the pilot length with K=1 and
%     p=1_tau.
% -------------------------------------------------------------------------
% Author: Italo Atzeni
% Version: 1.0
% Last edited: 25 Jul. 2021
% -------------------------------------------------------------------------
% License: This code is licensed under the GPLv2 license. If you use this
%     code in any way for research that results in publications, please
%     cite the above article.
% -------------------------------------------------------------------------

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K=1;
rho_dB=(-20:0.001:30)';
rho=db2pow(rho_dB);
tau=(32:1:512)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MSE_blm_p1=zeros(length(tau),length(rho));
for t=1:length(tau)
    tau(t)
    
    for r=1:length(rho)
        MSE_blm_p1(t,r)=compute_MSE_blm_p1(rho(r),tau(t));
    end
end

[MSE_blm_p1_opt,r_opt]=min(MSE_blm_p1,[],2);
rho_dB_opt=rho_dB(r_opt);
rho_opt=rho(r_opt);

res=2/pi*rho_opt./sqrt(1+2*rho_opt)-Omega(rho_opt./(rho_opt+1))-1./(tau-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename=['MSE_VS_tau_opt_K=' num2str(K) '_p1.mat'];
save(['files_mat/' filename],'K','rho_dB','tau','MSE_blm_p1_opt','rho_dB_opt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on;
grid on;
box on;
plot(tau,MSE_blm_p1_opt,'k');
xlim([64,512]);
xticks(64:64:512);
xlabel('$\tau$','interpreter','latex');
ylabel('Optimal MSE','interpreter','latex');
legend({'$\textrm{MSE}_{\textrm{BLM}}$'},'Location','NorthEast','Interpreter','latex');

figure
hold on;
grid on;
box on;
plot(tau,rho_dB_opt,'k');
xlim([64,512]);
xticks(64:64:512);
xlabel('$\tau$','interpreter','latex');
ylabel('Optimal $\rho$ [dB]','interpreter','latex');
legend({'$\rho^{\star}$'},'Location','NorthEast','Interpreter','latex');