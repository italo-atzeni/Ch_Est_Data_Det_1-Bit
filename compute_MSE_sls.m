% This MATLAB function was developed to generate numerical results for:
%
% Italo Atzeni and Antti Tölli, "Channel Estimation and Data Detection
%     Analysis of Massive MIMO with 1-Bit ADCs," IEEE Trans. Wireless
%     Commun. (to appear), 2021.
% -------------------------------------------------------------------------
% Description: Eq. (20) --> Compute the MSE of the scaled LS estimator.
% -------------------------------------------------------------------------
% Author: Italo Atzeni
% Version: 1.0
% Last edited: 25 Jul. 2021
% -------------------------------------------------------------------------
% License: This code is licensed under the GPLv2 license. If you use this
%     code in any way for research that results in publications, please
%     cite the above article.
% -------------------------------------------------------------------------

function MSE_sls=compute_MSE_sls(K,rho,tau,Delta)

MSE_sls=1-2/pi*rho*(2/pi*rho*(tau-K)+rho*K+1)^(-2)*(4/pi*rho*tau*(tau-K)+(rho*K+1)*(tau-Delta));