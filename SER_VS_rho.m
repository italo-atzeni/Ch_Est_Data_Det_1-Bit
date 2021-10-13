% This MATLAB function was developed to generate numerical results for:
%
% Italo Atzeni and Antti Tölli, "Channel Estimation and Data Detection
%     Analysis of Massive MIMO with 1-Bit ADCs," IEEE Trans. Wireless
%     Commun. (to appear), 2021.
% -------------------------------------------------------------------------
% Description: Fig. 8 --> Plot the SER against the transmit SNR.
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
M=128;
K=1;
rho_dB=(-10:40)';
rho=db2pow(rho_dB);
tau=32;
L=16;
I=10^4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p=compute_p_star(tau);
p_p1=ones(tau,1);

s=qammod(0:L-1,L)/sqrt(10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SER=zeros(length(rho),1);
SER_p1=zeros(length(rho),1);
for rr=1:length(rho)
    tic
    rho_dB(rr)

    Delta=compute_delta(K,rho(rr),tau,p,K);
    Psi_prime=compute_Psi_prime(K,rho(rr),tau,Delta);
    Delta_p1=tau*(tau-1)*Omega(rho(rr)/(rho(rr)+1));
    Psi_prime_p1=compute_Psi_prime(K,rho(rr),tau,Delta_p1);
    
    rng(0);
    s_hat=zeros(I,L);
    s_hat_1=zeros(I,L);
    parfor i=1:I
        h=(randn(M,1)+1i*randn(M,1))/sqrt(2);
        
        % channel estimation
        Zp=(randn(M,tau)+1i*randn(M,tau))/sqrt(2);
        Yp=sqrt(rho(rr))*h*p'+Zp;
        Yp_p1=sqrt(rho(rr))*h*p_p1'+Zp;
        Rp=sqrt((rho(rr)+1)/2)*(sign(real(Yp))+1i*sign(imag(Yp)));
        Rp_p1=sqrt((rho(rr)+1)/2)*(sign(real(Yp_p1))+1i*sign(imag(Yp_p1)));
        h_hat=sqrt(Psi_prime)*Rp*p;
        h_hat_p1=sqrt(Psi_prime_p1)*Rp_p1*p_p1;
        
        % MRC
        v=h_hat;
        v_p1=h_hat_p1;
        
        % data transmission
        z=(randn(M,1)+1i*randn(M,1))/sqrt(2);
        for l=1:L
            y=sqrt(rho(rr))*h*s(l)+z;
            r=sqrt((rho(rr)+1)/2)*(sign(real(y))+1i*sign(imag(y)));
            s_hat(i,l)=v'*r;
            s_hat_1(i,l)=v_p1'*r;
        end
    end
    
    % expected value of the estimated symbols
    s_hat_m=zeros(1,L);
    s_hat_m_p1=zeros(1,L);
    for l=1:L
        s_hat_m(l)=compute_El(M,rho(rr),tau,Delta,p,s(l));
        s_hat_m_p1(l)=compute_El(M,rho(rr),tau,Delta_p1,p_p1,s(l));
    end
    
    % SER
    d=zeros(I*L,L);
    d_p1=zeros(I*L,L);
    ind=kron(1:L,ones(1,I)).';
    for l=1:L
        d(:,l)=abs(s_hat(:)-s_hat_m(l)).^2;
        d_p1(:,l)=abs(s_hat_1(:)-s_hat_m_p1(l)).^2;
    end
    [~,ind_hat]=min(d,[],2);
    [~,ind_hat_p1]=min(d_p1,[],2);
    SER(rr)=sum(ind_hat~=ind(:))/(I*L);
    SER_p1(rr)=sum(ind_hat_p1~=ind(:))/(I*L);
    toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename=['SER_VS_rho_M=' num2str(M) '_tau=' num2str(tau) '.mat'];
save(['files_mat/' filename],'M','rho_dB','tau','SER','SER_p1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
semilogy(rho_dB,SER,'k-');
box on;
hold on;
grid on;
semilogy(rho_dB,SER_p1,'k--');
xlim([-10,40]);
xlabel('$\rho$ [dB]','interpreter','latex');
ylabel('$\textrm{SER}$','interpreter','latex');
title(['$M=$ ' int2str(M) ', $\tau=$ ' int2str(tau)],'interpreter','latex');
legend({'$\mathbf{p} = \mathbf{d}_{2}$','$\mathbf{p} = \mathbf{1}_{\tau}$'},'Location','NorthEast','Interpreter','latex');