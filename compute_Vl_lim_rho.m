% This MATLAB function was developed to generate numerical results for:
%
% Italo Atzeni and Antti Tölli, "Channel Estimation and Data Detection
%     Analysis of Massive MIMO with 1-Bit ADCs," IEEE Trans. Wireless
%     Commun. (to appear), 2021.
% -------------------------------------------------------------------------
% Description: Eq. (46) --> Compute the variance of the estimated symbol
%     (normalized by the transmit SNR) in the limit of rho->infty.
% -------------------------------------------------------------------------
% Author: Italo Atzeni
% Version: 1.0
% Last edited: 25 Jul. 2021
% -------------------------------------------------------------------------
% License: This code is licensed under the GPLv2 license. If you use this
%     code in any way for research that results in publications, please
%     cite the above article.
% -------------------------------------------------------------------------

function s_hat_v_lim_rho=compute_Vl_lim_rho(M,tau,Delta_bar,s_hat_m_lim_rho)

s_hat_v_lim_rho=2/pi*M*tau^2/(tau+Delta_bar)-1/M*abs(s_hat_m_lim_rho)^2;