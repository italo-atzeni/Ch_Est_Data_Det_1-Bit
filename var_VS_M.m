% This MATLAB function was developed to generate numerical results for:
%
% Italo Atzeni and Antti Tölli, "Channel Estimation and Data Detection
%     Analysis of Massive MIMO with 1-Bit ADCs," IEEE Trans. Wireless
%     Commun. (to appear), 2021.
% -------------------------------------------------------------------------
% Description: Fig. 5 --> Plot the normalized variance of the estimated
%     symbols against the number of BS antennas.
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
M=(32:256)';
rho_dB=10;
rho=db2pow(rho_dB);
tau=32;
L=16;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p=compute_p_star(tau);

s=qammod(0:L-1,L)/sqrt(10);

Delta=compute_delta(K,rho,tau,p,K);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s_hat_m=zeros(length(M),L);
s_hat_v=zeros(length(M),L);
for m=1:length(M)
    M(m)
    
    for l=1:L
        % expected value of the estimated symbols
        s_hat_m(m,l)=compute_El(M(m),rho,tau,Delta,p,s(l));
        
        % variance of the estimated symbols
        s_hat_v(m,l)=compute_Vl(M(m),rho,tau,Delta,s_hat_m(m,l));
    end
end
s_hat_v_norm=s_hat_v./abs(s_hat_m).^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename=['var_VS_M_rho=' num2str(rho_dB) 'dB_tau=' num2str(tau) '.mat'];
save(['files_mat/' filename],'M','rho_dB','tau','s_hat_v_norm');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on;
grid on;
box on;
plot(M,s_hat_v_norm(:,14),'k-');
plot(M,s_hat_v_norm(:,10),'b-');
plot(M,s_hat_v_norm(:,9),'r-');
xlim([32,256]);
xlabel('$M$','interpreter','latex');
ylabel('$V_{\ell}/|E_{\ell}|^{2}$','interpreter','latex');
legend({'$(1+j)/\sqrt{10}$','$(3+j)/\sqrt{10}$','$(3+j3)/\sqrt{10}$'},'Location','NorthEast','Interpreter','latex');
title(['$\rho=$ ' int2str(rho_dB) ' dB, $\tau=$ ' int2str(tau)],'interpreter','latex');