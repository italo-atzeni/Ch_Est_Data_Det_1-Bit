% This MATLAB function was developed to generate numerical results for:
%
% Italo Atzeni and Antti Tölli, "Channel Estimation and Data Detection
%     Analysis of Massive MIMO with 1-Bit ADCs," IEEE Trans. Wireless
%     Commun. (to appear), 2021.
% -------------------------------------------------------------------------
% Description: Fig. 5 --> Plot the normalized variance of the estimated
%     symbols against the number of BS antennas with p=1_tau.
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
D=dftmtx(tau);
p=D(:,2);

s=qammod(0:L-1,L)/sqrt(10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s_hat_v_norm_p1=zeros(length(M),L);
for m=1:length(M)
    M(m)
    
    for l=1:L
        % normalized variance of the estimated symbols
        s_hat_v_norm_p1(m,l)=compute_Vl_norm_p1(M(m),rho,tau,s(l));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename=['var_VS_M_rho=' num2str(rho_dB) 'dB_tau=' num2str(tau) '_p1.mat'];
save(['files_mat/' filename],'M','rho_dB','tau','s_hat_v_norm_p1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on;
grid on;
box on;
plot(M,s_hat_v_norm_p1(:,14),'k-');
plot(M,s_hat_v_norm_p1(:,10),'b-');
plot(M,s_hat_v_norm_p1(:,9),'r-');
xlim([32,256]);
xlabel('$M$','interpreter','latex');
ylabel('$V_{\ell}/|E_{\ell}|^{2}$','interpreter','latex');
legend({'$(1+j)/\sqrt{10}$','$(3+j)/\sqrt{10}$','$(3+j3)/\sqrt{10}$'},'Location','NorthEast','Interpreter','latex');
title(['$\rho=$ ' int2str(rho_dB) ' dB, $\tau=$ ' int2str(tau)],'interpreter','latex');