% This MATLAB function was developed to generate numerical results for:
%
% Italo Atzeni and Antti Tölli, "Channel Estimation and Data Detection
%     Analysis of Massive MIMO with 1-Bit ADCs," IEEE Trans. Wireless
%     Commun. (to appear), 2021.
% -------------------------------------------------------------------------
% Description: Fig. 2 --> Plot the MSE of the channel estimation against
%     the transmit SNR with K=1 and p=1_tau.
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

rng(0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K=1;
rho_dB=(-10:40)';
rho=db2pow(rho_dB);
tau=32;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MSE_blm_p1=zeros(length(rho),1);
MSE_sls_p1=zeros(length(rho),1);
for r=1:length(rho)
    rho_dB(r)
    
    MSE_blm_p1(r)=compute_MSE_blm_p1(rho(r),tau);    
    MSE_sls_p1(r)=compute_MSE_sls_p1(rho(r),tau);
end

MSE_blm_lim_p1=ones(length(rho),1)*compute_MSE_blm_lim_rho_p1;
MSE_sls_lim_p1=ones(length(rho),1)*compute_MSE_sls_lim_rho_p1(tau);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename=['MSE_VS_rho_K=' num2str(K) '_tau=' num2str(tau) '_p1.mat'];
save(['files_mat/' filename],'K','rho_dB','tau','MSE_blm_p1','MSE_blm_lim_p1','MSE_sls_p1','MSE_sls_lim_p1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on;
grid on;
box on;
plot(rho_dB,MSE_blm_p1,'g-');
plot(rho_dB,MSE_blm_lim_p1,'g:');
plot(rho_dB,MSE_sls_p1,'k-');
plot(rho_dB,MSE_sls_lim_p1,'k:');
xlim([-10,40]);
xticks(-10:10:40);
xlabel('$\rho$ [dB]','interpreter','latex');
ylabel('MSE of the channel estimation','interpreter','latex');
title(['$\tau=$ ' int2str(tau) ', $\mathbf{p} = \mathbf{1}_{\tau}$'],'interpreter','latex');
legend({'$\textrm{MSE}_{\textrm{BLM}}$','$\lim_{\rho \to \infty} \textrm{MSE}_{\textrm{BLM}}$','$\textrm{MSE}_{\textrm{SLS}}$', ...
    '$\lim_{\rho \to \infty} \textrm{MSE}_{\textrm{SLS}}$'},'Location','NorthWest','Interpreter','latex');