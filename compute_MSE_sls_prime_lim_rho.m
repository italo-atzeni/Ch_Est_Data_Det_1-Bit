% This MATLAB function was developed to generate numerical results for:
%
% Italo Atzeni and Antti Tölli, "Channel Estimation and Data Detection
%     Analysis of Massive MIMO with 1-Bit ADCs," IEEE Trans. Wireless
%     Commun. (to appear), 2021.
% -------------------------------------------------------------------------
% Description: Eq. (28) --> Compute the MSE of the optimal scaled LS
%     estimator in the limit of rho->infty.
% -------------------------------------------------------------------------
% Author: Italo Atzeni
% Version: 1.0
% Last edited: 25 Jul. 2021
% -------------------------------------------------------------------------
% License: This code is licensed under the GPLv2 license. If you use this
%     code in any way for research that results in publications, please
%     cite the above article.
% -------------------------------------------------------------------------

function MSE_sls_prime_lim_rho=compute_MSE_sls_prime_lim_rho(K,tau,Delta_bar)

MSE_sls_prime_lim_rho=1-2/pi*tau^2/(K*(tau+Delta_bar));