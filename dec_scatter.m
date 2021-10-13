% This MATLAB function was developed to generate numerical results for:
%
% Italo Atzeni and Antti Tölli, "Channel Estimation and Data Detection
%     Analysis of Massive MIMO with 1-Bit ADCs," IEEE Trans. Wireless
%     Commun. (to appear), 2021.
% -------------------------------------------------------------------------
% Description: Fig. 4 --> Plot the estimated symbols with the MRC receiver
%     and the corresponding expected values.
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

rng(0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=128;
K=1;
rho_dB=0;
rho=db2pow(rho_dB);
tau=32;
L=16;
I=10^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p=compute_p_star(tau);

h=(randn(M,I)+1i*randn(M,I))/sqrt(2);
Zp=(randn(M,tau,I)+1i*randn(M,tau,I))/sqrt(2);
z=(randn(M,I)+1i*randn(M,I))/sqrt(2);

s=qammod(0:L-1,L)/sqrt(10);

Delta=compute_delta(K,rho,tau,p,K);
Psi_prime=compute_Psi_prime(K,rho,tau,Delta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s_hat=zeros(I,L);
for i=1:I
    % channel estimation
    Yp=sqrt(rho)*h(:,i)*p'+Zp(:,:,i);
    Rp=sqrt((rho+1)/2)*(sign(real(Yp))+1i*sign(imag(Yp)));
    h_hat=sqrt(Psi_prime)*Rp*p;
    
    % MRC
    v=h_hat;
    
    % data transmission
    for l=1:L
        y=sqrt(rho)*h(:,i)*s(l)+z(:,i);
        r=sqrt((rho+1)/2)*(sign(real(y))+1i*sign(imag(y)));
        s_hat(i,l)=v'*r;
    end
end
s_hat_real=real(s_hat(:));
s_hat_imag=imag(s_hat(:));
s_hat_m_mc=mean(s_hat,1);
s_hat_v_mc=var(s_hat,1);

s_hat_m=zeros(1,L);
s_hat_v=zeros(1,L);
for l=1:L
    % expected value of the estimated symbols
    s_hat_m(l)=compute_El(M,rho,tau,Delta,p,s(l));
    
    % variance of the estimated symbols
    s_hat_v(l)=compute_Vl(M,rho,tau,Delta,s_hat_m(l));
end
s_hat_m_real=real(s_hat_m);
s_hat_m_imag=imag(s_hat_m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename=['dec_scatter_M=' num2str(M) '_rho=' num2str(rho_dB) 'dB_tau=' num2str(tau) '.mat'];
save(['files_mat/' filename],'M','rho_dB','tau','s_hat_real','s_hat_imag','s_hat_m_real','s_hat_m_imag');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on;
grid on;
box on;
axis equal;
scatter(10^6,10^6,'MarkerEdgeColor','b');
scatter(10^6,10^6,50,'filled','MarkerEdgeColor','m','MarkerFaceColor','m');
scatter(s_hat_real,s_hat_imag,'MarkerEdgeColor','b');
scatter(s_hat_m_real,s_hat_m_imag,50,'filled','MarkerEdgeColor','m','MarkerFaceColor','m');
xlim([-1.5*real(s_hat_m(12)),1.5*real(s_hat_m(12))]);
ylim([-1.5*real(s_hat_m(12)),1.5*real(s_hat_m(12))]);
xlabel('I');
ylabel('Q');
legend({'$\hat{s}_{\ell}$','$E_{\ell}$'},'Location','SouthEast','Interpreter','latex');
title(['$M=$ ' int2str(M) ', $\rho=$ ' int2str(rho_dB) ' dB, $\tau=$ ' int2str(tau)],'interpreter','latex');