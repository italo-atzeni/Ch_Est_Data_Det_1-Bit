% This MATLAB function was developed to generate numerical results for:
%
% Italo Atzeni and Antti Tölli, "Channel Estimation and Data Detection
%     Analysis of Massive MIMO with 1-Bit ADCs," IEEE Trans. Wireless
%     Commun. (to appear), 2021.
% -------------------------------------------------------------------------
% Description: Fig. 7 --> Plot the normalized variance of the estimated
%     symbols against the pilot length with p=1_tau.
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
M=128;
rho_dB=10;
rho=db2pow(rho_dB);
tau=(1:128)';
L=16;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=qammod(0:L-1,L)/sqrt(10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s_hat_v_norm_p1=zeros(length(tau),L);
for t=1:length(tau)
    tau(t)
    
    for l=1:L
        % normalized variance of the estimated symbols
        s_hat_v_norm_p1(t,l)=compute_Vl_norm_p1(M,rho,tau(t),s(l));
    end
end

s_hat_v_norm_lim_p1=zeros(length(tau),L);
for l=1:L
    % normalied variance of the estimated symbols (lim_{tau->infty})
    s_hat_v_norm_lim_p1(:,l)=ones(length(tau),1)*compute_Vl_norm_lim_tau_p1(M,rho,s(l));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename=['var_VS_tau_M=' num2str(M) '_rho=' num2str(rho_dB) 'dB_p1.mat'];
save(['files_mat/' filename],'M','rho_dB','tau','s_hat_v_norm_p1','s_hat_v_norm_lim_p1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on;
grid on;
box on;
plot(tau,s_hat_v_norm_p1(:,14),'k-');
plot(tau,s_hat_v_norm_p1(:,10),'b-');
plot(tau,s_hat_v_norm_p1(:,9),'r-');
plot(tau,s_hat_v_norm_lim_p1(:,14),'k:');
plot(tau,s_hat_v_norm_lim_p1(:,10),'b:');
plot(tau,s_hat_v_norm_lim_p1(:,9),'r:');
xlim([0,128]);
xlabel('$\tau$','interpreter','latex');
ylabel('$V_{\ell}/|E_{\ell}|^{2}$','interpreter','latex');
legend({'$(1+j)/\sqrt{10}$','$(3+j)/\sqrt{10}$','$(3+j3)/\sqrt{10}$'},'Location','NorthEast','Interpreter','latex');
title(['$M=$ ' int2str(M) ', $\rho=$ ' int2str(rho_dB) ' dB'],'interpreter','latex');