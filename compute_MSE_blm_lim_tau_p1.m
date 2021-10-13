% This MATLAB function was developed to generate numerical results for:
%
% Italo Atzeni and Antti Tölli, "Channel Estimation and Data Detection
%     Analysis of Massive MIMO with 1-Bit ADCs," IEEE Trans. Wireless
%     Commun. (to appear), 2021.
% -------------------------------------------------------------------------
% Description: Eq. (35) --> Compute the MSE of the BLM estimator in the
%     limit of tau->infty with p=1_tau.
% -------------------------------------------------------------------------
% Author: Italo Atzeni
% Version: 1.0
% Last edited: 25 Jul. 2021
% -------------------------------------------------------------------------
% License: This code is licensed under the GPLv2 license. If you use this
%     code in any way for research that results in publications, please
%     cite the above article.
% -------------------------------------------------------------------------

function MSE_blm_lim_tau_p1=compute_MSE_blm_lim_tau_p1(rho)

MSE_blm_lim_tau_p1=1-2/pi*rho/((rho+1)*Omega(rho/(rho+1)));