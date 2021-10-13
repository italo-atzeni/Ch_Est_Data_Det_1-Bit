% This MATLAB function was developed to generate numerical results for:
%
% Italo Atzeni and Antti Tölli, "Channel Estimation and Data Detection
%     Analysis of Massive MIMO with 1-Bit ADCs," IEEE Trans. Wireless
%     Commun. (to appear), 2021.
% -------------------------------------------------------------------------
% Description: Eq. (68).
% -------------------------------------------------------------------------
% Author: Italo Atzeni
% Version: 1.0
% Last edited: 25 Jul. 2021
% -------------------------------------------------------------------------
% License: This code is licensed under the GPLv2 license. If you use this
%     code in any way for research that results in publications, please
%     cite the above article.
% -------------------------------------------------------------------------

function Phi=compute_Phi(K,rho,tau,P)

Phi=zeros(tau,tau);
for u=1:tau
    for v=1:tau
        if u==v
            Phi(u,v)=1;
        else
            temp=zeros(K,1);
            for i=1:K
                temp(i)=P(u,i)*P(v,i)';
            end
            Phi(u,v)=Omega(rho*sum(real(temp))/(rho*K+1))-1i*Omega(rho*sum(imag(temp))/(rho*K+1));
        end
    end
end