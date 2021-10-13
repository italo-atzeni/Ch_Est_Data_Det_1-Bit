% This MATLAB function was developed to generate numerical results for:
%
% Italo Atzeni and Antti Tölli, "Channel Estimation and Data Detection
%     Analysis of Massive MIMO with 1-Bit ADCs," IEEE Trans. Wireless
%     Commun. (to appear), 2021.
% -------------------------------------------------------------------------
% Description: Eq. (36) --> Compute the MSE of the scaled LS estimator with
%     p=1_tau.
% -------------------------------------------------------------------------
% Author: Italo Atzeni
% Version: 1.0
% Last edited: 25 Jul. 2021
% -------------------------------------------------------------------------
% License: This code is licensed under the GPLv2 license. If you use this
%     code in any way for research that results in publications, please
%     cite the above article.
% -------------------------------------------------------------------------

function MSE_sls_p1=compute_MSE_sls_p1(rho,tau)

MSE_sls_p1=1-2/pi*rho*tau*(2/pi*rho*(tau-1)+rho+1)^(-2)*(4/pi*rho*(tau-1)+(rho+1)*(1-(tau-1)*Omega(rho/(rho+1))));