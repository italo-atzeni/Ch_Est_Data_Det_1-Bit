% This MATLAB function was developed to generate numerical results for:
%
% Italo Atzeni and Antti Tölli, "Channel Estimation and Data Detection
%     Analysis of Massive MIMO with 1-Bit ADCs," IEEE Trans. Wireless
%     Commun. (to appear), 2021.
% -------------------------------------------------------------------------
% Description: Eq. (43) --> Compute the expected value of the estimated
%     symbol.
% -------------------------------------------------------------------------
% Author: Italo Atzeni
% Version: 1.0
% Last edited: 25 Jul. 2021
% -------------------------------------------------------------------------
% License: This code is licensed under the GPLv2 license. If you use this
%     code in any way for research that results in publications, please
%     cite the above article.
% -------------------------------------------------------------------------

function s_hat_m=compute_El(M,rho,tau,Delta,p,s)

temp=zeros(tau,1);
for u=1:tau
    temp(u)=Omega(rho*real(p(u)*s)/sqrt((rho+1)*(rho*abs(s)^2+1)))+1i*Omega(rho*imag(p(u)*s)/sqrt((rho+1)*(rho*abs(s)^2+1)));
end
s_hat_m=sqrt(2/pi*rho)*M*tau/(tau+Delta)*p'*temp;