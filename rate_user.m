function [Ri] = rate_user(H, T, P)
% Function to obtain the user rates based on the SNR, the normalized precoder and 
% the channel matrix. 
% The rate is obtained following the well-known expression of the capacity under
% gaussian assumption Rate = log2(1 + SNR) bits

    K  = size(H,1);  % Number of receivers
    Ri = zeros(1,K); % Initializing rate
    
    for i = 1:K
        Ri(i) = log2(1 + P*abs(H(i,:)*T(:,i)).^2/(1 + P*sum(abs(H(i,:)*T(:,[1:i-1,i+1:end])).^2)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Antonio Bazco-Nogueras  
% Date: 2022/05/11
% Contact: antonio.bazco@imdea.org
% License: This file can be distributed, remixed, adapted, and other work can be
% built upon it, as long as appropiate credit is included for the original creation. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%