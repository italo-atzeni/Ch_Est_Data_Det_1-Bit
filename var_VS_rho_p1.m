% This MATLAB function was developed to generate numerical results for:
%
% Italo Atzeni and Antti Tölli, "Channel Estimation and Data Detection
%     Analysis of Massive MIMO with 1-Bit ADCs," IEEE Trans. Wireless
%     Commun. (to appear), 2021.
% -------------------------------------------------------------------------
% Description: Fig. 6 --> Plot the normalized variance of the estimated
%     symbols against the transmit SNR with p=1_tau.
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
rho_dB=(-10:40)';
rho=db2pow(rho_dB);
tau=32;
L=16;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=qammod(0:L-1,L)/sqrt(10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s_hat_v_norm_p1=zeros(length(rho),L);
for r=1:length(rho)
    rho_dB(r)
    
    for l=1:L
        % normalized variance of the estimated symbols
        s_hat_v_norm_p1(r,l)=compute_Vl_norm_p1(M,rho(r),tau,s(l));
    end
end

s_hat_v_norm_lim_p1=zeros(length(rho),L);
for l=1:L
    % normalied variance of the estimated symbols (lim_{rho->infty})
    s_hat_v_norm_lim_p1(:,l)=ones(length(rho),1)*compute_Vl_norm_lim_rho_p1(M,s(l));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename=['var_VS_rho_M=' num2str(M) '_tau=' num2str(tau) '_p1.mat'];
save(['files_mat/' filename],'M','rho_dB','tau','s_hat_v_norm_p1','s_hat_v_norm_lim_p1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
semilogy(rho_dB,s_hat_v_norm_p1(:,14),'k-');
hold on;
grid on;
box on;
semilogy(rho_dB,s_hat_v_norm_p1(:,10),'b-');
semilogy(rho_dB,s_hat_v_norm_p1(:,9),'r-');
semilogy(rho_dB,s_hat_v_norm_lim_p1(:,14),'k:');
semilogy(rho_dB,s_hat_v_norm_lim_p1(:,10),'b:');
semilogy(rho_dB,s_hat_v_norm_lim_p1(:,9),'r:');
xlim([-10,40]);
xlabel('$\rho$ [dB]','interpreter','latex');
ylabel('$V_{\ell}/|E_{\ell}|^{2}$','interpreter','latex');
legend({'$(1+j)/\sqrt{10}$','$(3+j)/\sqrt{10}$','$(3+j3)/\sqrt{10}$'},'Location','NorthEast','Interpreter','latex');
title(['$M=$ ' int2str(M) ', $\tau=$ ' int2str(tau)],'interpreter','latex');