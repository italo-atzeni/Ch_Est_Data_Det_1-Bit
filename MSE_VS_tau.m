% This MATLAB function was developed to generate numerical results for:
%
% Italo Atzeni and Antti Tölli, "Channel Estimation and Data Detection
%     Analysis of Massive MIMO with 1-Bit ADCs," IEEE Trans. Wireless
%     Commun. (to appear), 2021.
% -------------------------------------------------------------------------
% Description: Fig. 3 --> Plot the MSE of the channel estimation against
%     the pilot length.
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
rho_dB=10;
rho=db2pow(rho_dB);
K_equal=1;
if K_equal==1
    K_0=4;
    tau=(K_0:4:64)';
    K=K_0*ones(length(tau),1);
elseif K_equal==0
    c=4;
    tau=(c:c:128)';
    K=tau/c;
end

mc=0;
M=K;
I=10^6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MSE_blm=zeros(length(tau),1);
MSE_blm_mc=zeros(length(tau),1);
MSE_sls=zeros(length(tau),1);
MSE_sls_mc=zeros(length(tau),1);
MSE_sls_prime=zeros(length(tau),1);
MSE_sls_prime_mc=zeros(length(tau),1);
for t=1:length(tau)
    tau(t)
    
    if tau(t)==1
        P=1;
    else
        if K(t)==1
            P=compute_p_star(tau(t));
        else
            D=dftmtx(tau(t));
            P=D(:,1:K(t));
        end
    end
    P_bar=kron(P,eye(M));
    
    delta=zeros(K(t),1);
    for k=1:K(t)
        delta(k)=compute_delta(K(t),rho,tau(t),P,k);
    end
    Delta=1/K(t)*sum(delta);
    
    MSE_blm(t)=compute_MSE_blm(K(t),rho,tau(t),delta);
    MSE_sls(t)=compute_MSE_sls(K(t),rho,tau(t),Delta);
    MSE_sls_prime(t)=compute_MSE_sls_prime(K(t),rho,tau(t),Delta);
    
    if mc==1
        Phi=compute_Phi(K(t),rho,tau(t),P);
        Sigmap=(rho*K(t)+1)*kron(Phi,eye(M));
        Psi=compute_Psi(K(t),rho,tau(t));
        Psi_prime=compute_Psi_prime(K(t),rho,tau(t),Delta);
        
        err_blm=zeros(I,1);
        err_sls=zeros(I,1);
        err_sls_prime=zeros(I,1);
        for i=1:I
            H=(randn(M,K(t))+1i*randn(M,K(t)))/sqrt(2);
            Zp=(randn(M,tau(t))+1i*randn(M,tau(t)))/sqrt(2);

            Yp=sqrt(rho)*H*P'+Zp;
            Rp=sqrt((rho*K(t)+1)/2)*(sign(real(Yp))+1i*sign(imag(Yp)));

            H_hat_blm=sqrt(2/pi*rho)*reshape(P_bar.'*(Sigmap\Rp(:)),[M,K(t)]);
            H_hat_sls=sqrt(Psi)*Rp*P;
            H_hat_sls_prime=sqrt(Psi_prime)*Rp*P;
            
            err_blm(i)=1/(M*K(t))*norm(H_hat_blm-H,'fro')^2;
            err_sls(i)=1/(M*K(t))*norm(H_hat_sls-H,'fro')^2;
            err_sls_prime(i)=1/(M*K(t))*norm(H_hat_sls_prime-H,'fro')^2;
        end
        MSE_blm_mc(t)=mean(err_blm);
        MSE_sls_mc(t)=mean(err_sls);
        MSE_sls_prime_mc(t)=mean(err_sls_prime);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if K_equal==1
    filename=['MSE_VS_tau_K=' num2str(K_0) '_rho=' num2str(rho_dB) 'dB.mat'];
elseif K_equal==0
    filename=['MSE_VS_tau_c=' num2str(c) '_rho=' num2str(rho_dB) 'dB.mat'];
end
save(['files_mat/' filename],'K','rho_dB','tau','MSE_blm','MSE_blm_mc','MSE_sls','MSE_sls_mc','MSE_sls_prime','MSE_sls_prime_mc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on;
grid on;
box on;
if mc==0
    plot(tau,MSE_blm,'gx');
    plot(tau,MSE_sls,'k-');
    plot(tau,MSE_sls_prime,'r-');
    legend({'$\textrm{MSE}_{\textrm{BLM}}$','$\textrm{MSE}_{\textrm{SLS}}$','$\textrm{MSE}_{\textrm{SLS}}^{\prime}$'},'Location','NorthEast','Interpreter','latex');
elseif mc==1    
    plot(tau,MSE_blm,'gx');
    plot(tau,MSE_blm_mc,'go');
    plot(tau,MSE_sls,'k-');
    plot(tau,MSE_sls_mc,'k+');
    plot(tau,MSE_sls_prime,'r-');
    plot(tau,MSE_sls_prime_mc,'r+');
    legend({'$\textrm{MSE}_{\textrm{BLM}}$','$\textrm{MSE}_{\textrm{BLM}}$ (MC)','$\textrm{MSE}_{\textrm{SLS}}$','$\textrm{MSE}_{\textrm{SLS}}$ (MC)', ...
        '$\textrm{MSE}_{\textrm{SLS}}^{\prime}$','$\textrm{MSE}_{\textrm{SLS}}^{\prime}$ (MC)'},'Location','NorthEast','Interpreter','latex');
end
xlim([0,128]);
xticks(0:16:128);
xlabel('$\tau$','interpreter','latex');
ylabel('MSE of the channel estimation','interpreter','latex');
if K_equal==1
    title(['$K=$ ' int2str(K_0) ', $\rho=$ ' int2str(rho_dB) ' dB'],'interpreter','latex');
elseif K_equal==0
    title(['$c=$ ' int2str(c) ', $\rho=$ ' int2str(rho_dB) ' dB'],'interpreter','latex');
end