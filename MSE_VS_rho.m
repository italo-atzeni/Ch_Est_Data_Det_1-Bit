% This MATLAB function was developed to generate numerical results for:
%
% Italo Atzeni and Antti Tölli, "Channel Estimation and Data Detection
%     Analysis of Massive MIMO with 1-Bit ADCs," IEEE Trans. Wireless
%     Commun. (to appear), 2021.
% -------------------------------------------------------------------------
% Description: Fig. 2 --> Plot the MSE of the channel estimation against
%     the transmit SNR.
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
K=4;
rho_dB=(-10:40)';
rho=db2pow(rho_dB);
tau=32;

mc=0;
M=K;
I=10^6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if K==1
    P=compute_p_star(tau);
else
    D=dftmtx(tau);
    P=D(:,1:K);
end
P_tilde=kron(P,eye(M));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MSE_blm=zeros(length(rho),1);
MSE_blm_mc=zeros(length(rho),1);
MSE_sls=zeros(length(rho),1);
MSE_sls_mc=zeros(length(rho),1);
MSE_sls_prime=zeros(length(rho),1);
MSE_sls_prime_mc=zeros(length(rho),1);
for r=1:length(rho)
    rho_dB(r)
    
    delta=zeros(K,1);
    for k=1:K
        delta(k)=compute_delta(K,rho(r),tau,P,k);
    end
    Delta=1/K*sum(delta);
    
    MSE_blm(r)=compute_MSE_blm(K,rho(r),tau,delta);
    MSE_sls(r)=compute_MSE_sls(K,rho(r),tau,Delta);
    MSE_sls_prime(r)=compute_MSE_sls_prime(K,rho(r),tau,Delta);

    if mc==1
        Phi=compute_Phi(K,rho(r),tau,P);
        Sigmap=(rho(r)*K+1)*kron(Phi,eye(M));
        Psi=compute_Psi(K,rho(r),tau);
        Psi_prime=compute_Psi_prime(K,rho(r),tau,Delta);

        err_blm=zeros(I,1);
        err_sls=zeros(I,1);
        err_sls_prime=zeros(I,1);
        for i=1:I
            H=(randn(M,K)+1i*randn(M,K))/sqrt(2);
            Zp=(randn(M,tau)+1i*randn(M,tau))/sqrt(2);

            Yp=sqrt(rho(r))*H*P'+Zp;
            Rp=sqrt((rho(r)*K+1)/2)*(sign(real(Yp))+1i*sign(imag(Yp)));

            H_hat_blm=sqrt(2/pi*rho(r))*reshape(P_tilde.'*(Sigmap\Rp(:)),[M,K]);
            H_hat_sls=sqrt(Psi)*Rp*P;
            H_hat_sls_prime=sqrt(Psi_prime)*Rp*P;
            
            err_blm(i)=1/(M*K)*norm(H_hat_blm-H,'fro')^2;
            err_sls(i)=1/(M*K)*norm(H_hat_sls-H,'fro')^2;
            err_sls_prime(i)=1/(M*K)*norm(H_hat_sls_prime-H,'fro')^2;
            
        end
        MSE_blm_mc(r)=mean(err_blm);
        MSE_sls_mc(r)=mean(err_sls);
        MSE_sls_prime_mc(r)=mean(err_sls_prime);
    end
end
delta_bar=zeros(K,1);
for k=1:K
    delta_bar(k)=compute_delta_bar(K,tau,P,k);
end
Delta_bar=1/K*sum(delta_bar);

MSE_blm_lim=ones(length(rho),1)*compute_MSE_blm_lim_rho(K,tau,delta_bar);
MSE_sls_lim=ones(length(rho),1)*compute_MSE_sls_lim_rho(K,tau,Delta_bar);
MSE_sls_prime_lim=ones(length(rho),1)*compute_MSE_sls_prime_lim_rho(K,tau,Delta_bar);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename=['MSE_VS_rho_K=' num2str(K) '_tau=' num2str(tau) '.mat'];
save(['files_mat/' filename],'K','rho_dB','tau','MSE_blm','MSE_blm_mc','MSE_blm_lim','MSE_sls','MSE_sls_mc','MSE_sls_lim','MSE_sls_prime','MSE_sls_prime_mc', ...
    'MSE_sls_prime_lim');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on;
grid on;
box on;
if mc==0
    plot(rho_dB,MSE_blm,'gx');
    plot(rho_dB,MSE_blm_lim,'g:');
    plot(rho_dB,MSE_sls,'k-');
    plot(rho_dB,MSE_sls_lim,'k:');
    plot(rho_dB,MSE_sls_prime,'r-');
    plot(rho_dB,MSE_sls_prime_lim,'r:');
    legend({'$\textrm{MSE}_{\textrm{BLM}}$','$\lim_{\rho \to \infty} \textrm{MSE}_{\textrm{BLM}}$','$\textrm{MSE}_{\textrm{SLS}}$', ...
        '$\lim_{\rho \to \infty} \textrm{MSE}_{\textrm{SLS}}$','$\textrm{MSE}_{\textrm{SLS}}^{\prime}$', ...
        '$\lim_{\rho \to \infty} \textrm{MSE}_{\textrm{SLS}}^{\prime}$'},'Location','NorthEast','Interpreter','latex');
elseif mc==1
    plot(rho_dB,MSE_blm,'gx');
    plot(rho_dB,MSE_blm_mc,'go');
    plot(rho_dB,MSE_blm_lim,'g:');
    plot(rho_dB,MSE_sls,'k-');
    plot(rho_dB,MSE_sls_mc,'k+');
    plot(rho_dB,MSE_sls_lim,'k:');
    plot(rho_dB,MSE_sls_prime,'r-');
    plot(rho_dB,MSE_sls_prime_mc,'r+');
    plot(rho_dB,MSE_sls_prime_lim,'r:');
    legend({'$\textrm{MSE}_{\textrm{BLM}}$','$\textrm{MSE}_{\textrm{BLM}}$ (MC)','$\lim_{\rho \to \infty} \textrm{MSE}_{\textrm{BLM}}$', ...
        '$\textrm{MSE}_{\textrm{SLS}}$','$\textrm{MSE}_{\textrm{SLS}}$ (MC)','$\lim_{\rho \to \infty} \textrm{MSE}_{\textrm{SLS}}$', ...
        '$\textrm{MSE}_{\textrm{SLS}}^{\prime}$','$\textrm{MSE}_{\textrm{SLS}}^{\prime}$ (MC)', ...
        '$\lim_{\rho \to \infty} \textrm{MSE}_{\textrm{SLS}}^{\prime}$'},'Location','NorthEast','Interpreter','latex');
end
xlim([-10,40]);
xticks(-10:10:40);
xlabel('$\rho$ [dB]','interpreter','latex');
ylabel('MSE of the channel estimation','interpreter','latex');
title(['$K=$ ' int2str(K) ', $\tau=$ ' int2str(tau)],'interpreter','latex');