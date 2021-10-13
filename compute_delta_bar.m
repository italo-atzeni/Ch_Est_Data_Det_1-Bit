% This MATLAB function was developed to generate numerical results for:
%
% Italo Atzeni and Antti Tölli, "Channel Estimation and Data Detection
%     Analysis of Massive MIMO with 1-Bit ADCs," IEEE Trans. Wireless
%     Commun. (to appear), 2021.
% -------------------------------------------------------------------------
% Description: Eq. (29).
% -------------------------------------------------------------------------
% Author: Italo Atzeni
% Version: 1.0
% Last edited: 25 Jul. 2021
% -------------------------------------------------------------------------
% License: This code is licensed under the GPLv2 license. If you use this
%     code in any way for research that results in publications, please
%     cite the above article.
% -------------------------------------------------------------------------

function delta_bar=compute_delta_bar(K,tau,P,k)

A=zeros(tau,tau);
for u=1:tau
    for v=1:tau
        if u~=v
            temp=zeros(K,1);
            for i=1:K
                temp(i)=P(u,i)*P(v,i)';
            end
            A(u,v)=real(P(u,k)'*P(v,k))*Omega(sum(real(temp))/K)-imag(P(u,k)'*P(v,k))*Omega(sum(imag(temp))/K);
        end
    end
end
delta_bar=real(sum(sum(A)));