% This MATLAB function was developed to generate numerical results for:
%
% Italo Atzeni and Antti Tölli, "Channel Estimation and Data Detection
%     Analysis of Massive MIMO with 1-Bit ADCs," IEEE Trans. Wireless
%     Commun. (to appear), 2021.
% -------------------------------------------------------------------------
% Description: Fig. 6 --> Plot the normalized variance of the estimated
%     symbols against the transmit SNR.
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

p=compute_p_star(tau);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s_hat_m=zeros(length(rho),L);
s_hat_v=zeros(length(rho),L);
for r=1:length(rho)
    rho_dB(r)
    
    Delta=compute_delta(K,rho(r),tau,p,K);
    for l=1:L
        % expected value of the estimated symbols
        s_hat_m(r,l)=compute_El(M,rho(r),tau,Delta,p,s(l));
        
        % variance of the estimated symbols
        s_hat_v(r,l)=compute_Vl(M,rho(r),tau,Delta,s_hat_m(r,l));
    end
end
s_hat_v_norm=s_hat_v./abs(s_hat_m).^2;

Delta_bar=compute_delta_bar(K,tau,p,K);
s_hat_m_lim=zeros(1,L);
s_hat_v_lim=zeros(1,L);
for l=1:L
    % expected value of the estimated symbols (lim_{rho->infty})
    s_hat_m_lim(l)=compute_El_lim_rho(M,tau,Delta_bar,p,s(l));

    % variance of the estimated symbols (lim_{rho->infty})
    s_hat_v_lim(l)=compute_Vl_lim_rho(M,tau,Delta_bar,s_hat_m_lim(l));
end
s_hat_v_norm_lim=ones(length(rho),1)*s_hat_v_lim./abs(s_hat_m_lim).^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename=['var_VS_rho_M=' num2str(M) '_tau=' num2str(tau) '.mat'];
save(['files_mat/' filename],'M','rho_dB','tau','s_hat_v_norm','s_hat_v_norm_lim');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
semilogy(rho_dB,s_hat_v_norm(:,14),'k-');
hold on;
grid on;
box on;
semilogy(rho_dB,s_hat_v_norm(:,10),'b-');
semilogy(rho_dB,s_hat_v_norm(:,9),'r-');
semilogy(rho_dB,s_hat_v_norm_lim(:,14),'k:');
semilogy(rho_dB,s_hat_v_norm_lim(:,10),'b:');
semilogy(rho_dB,s_hat_v_norm_lim(:,9),'r:');
xlim([-10,40]);
xlabel('$\rho$ [dB]','interpreter','latex');
ylabel('$V_{\ell}/|E_{\ell}|^{2}$','interpreter','latex');
legend({'$(1+j)/\sqrt{10}$','$(3+j)/\sqrt{10}$','$(3+j3)/\sqrt{10}$'},'Location','NorthEast','Interpreter','latex');
title(['$M=$ ' int2str(M) ', $\tau=$ ' int2str(tau)],'interpreter','latex');