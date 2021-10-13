% This MATLAB function was developed to generate numerical results for:
%
% Italo Atzeni and Antti Tölli, "Channel Estimation and Data Detection
%     Analysis of Massive MIMO with 1-Bit ADCs," IEEE Trans. Wireless
%     Commun. (to appear), 2021.
% -------------------------------------------------------------------------
% Description: Eq. (48) --> Compute the normalized variance of the
%     estimated symbol with p=1_tau.
% -------------------------------------------------------------------------
% Author: Italo Atzeni
% Version: 1.0
% Last edited: 25 Jul. 2021
% -------------------------------------------------------------------------
% License: This code is licensed under the GPLv2 license. If you use this
%     code in any way for research that results in publications, please
%     cite the above article.
% -------------------------------------------------------------------------

function s_hat_v_norm_p1=compute_Vl_norm_p1(M,rho,tau,s)

s_hat_v_norm_p1=1/M*(1+(tau-1)*Omega(rho/(rho+1)))/tau*(Omega(rho*real(s)/sqrt((rho+1)*(rho*abs(s)^2+1)))^2+Omega(rho*imag(s)/sqrt((rho+1)*(rho*abs(s)^2+1)))^2)^(-1)-1/M;