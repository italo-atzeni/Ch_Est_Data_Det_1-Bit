% This MATLAB function was developed to generate numerical results for:
%
% Italo Atzeni and Antti Tölli, "Channel Estimation and Data Detection
%     Analysis of Massive MIMO with 1-Bit ADCs," IEEE Trans. Wireless
%     Commun. (to appear), 2021.
% -------------------------------------------------------------------------
% Description: Fig. 7 --> Plot the normalized variance of the estimated
%     symbols against the pilot length.
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
s_hat_m=zeros(length(tau),L);
s_hat_v=zeros(length(tau),L);
for t=1:length(tau)
    tau(t)
    
    if tau(t)==1
        p=1;
    else
        p=compute_p_star(tau(t));
    end
    
    Delta=compute_delta(K,rho,tau(t),p,K);
    
    for l=1:L
        % expected value of the estimated symbols
        s_hat_m(t,l)=compute_El(M,rho,tau(t),Delta,p,s(l));
        
        % variance of the estimated symbols
        s_hat_v(t,l)=compute_Vl(M,rho,tau(t),Delta,s_hat_m(t,l));
    end
end
s_hat_v_norm=s_hat_v./abs(s_hat_m).^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename=['var_VS_tau_M=' num2str(M) '_rho=' num2str(rho_dB) 'dB.mat'];
save(['files_mat/' filename],'M','rho_dB','tau','s_hat_v_norm');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on;
grid on;
box on;
plot(tau,s_hat_v_norm(:,14),'k-');
plot(tau,s_hat_v_norm(:,10),'b-');
plot(tau,s_hat_v_norm(:,9),'r-');
xlim([0,128]);
xlabel('$\tau$','interpreter','latex');
ylabel('$V_{\ell}/|E_{\ell}|^{2}$','interpreter','latex');
legend({'$(1+j)/\sqrt{10}$','$(3+j)/\sqrt{10}$','$(3+j3)/\sqrt{10}$'},'Location','NorthEast','Interpreter','latex');
title(['$M=$ ' int2str(M) ', $\rho=$ ' int2str(rho_dB) ' dB'],'interpreter','latex');