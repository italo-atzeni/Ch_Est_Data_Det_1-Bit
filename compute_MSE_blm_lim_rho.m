% This MATLAB function was developed to generate numerical results for:
%
% Italo Atzeni and Antti Tölli, "Channel Estimation and Data Detection
%     Analysis of Massive MIMO with 1-Bit ADCs," IEEE Trans. Wireless
%     Commun. (to appear), 2021.
% -------------------------------------------------------------------------
% Description: Eq. (26) --> Compute the MSE of the BLM estimator in the
%     limit of rho->infty.
% -------------------------------------------------------------------------
% Author: Italo Atzeni
% Version: 1.0
% Last edited: 25 Jul. 2021
% -------------------------------------------------------------------------
% License: This code is licensed under the GPLv2 license. If you use this
%     code in any way for research that results in publications, please
%     cite the above article.
% -------------------------------------------------------------------------

function MSEp_blm_lim_rho=compute_MSE_blm_lim_rho(K,tau,delta_bar)

MSEp_blm_lim_rho=1-2/pi*tau^2/K^2*sum(1./(tau+delta_bar));