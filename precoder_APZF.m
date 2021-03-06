function [W] = precoder_APZF(He_M, Mn, mu, P, pow_control)
% precoder_APZF: 
% Function that computes the AP-ZF precoding matrix. AP-ZF algorithm was developed 
% in the paper "Degrees of freedom of the network MIMO channel with distributed CSI,"
% P. D. Kerret and D. Gesbert, IEEE Trans. Inf. Theory, vol. 58, no. 11, pp. 6806–6824, Nov. 2012. 
%
% In this algorithm, there is a set of passive TXs, for which the precodin vectors
% do not depend on local CSIT and are simple precoders known by the other TXs (the so
% called active TXs). In turn, the active TXs use their local CSIT aiming at removing
% the interference generated by the transmission of the passive TXs. 
%
% INPUTS: He_M: Channel matrix estimates at each TX
%         Mn:   Num. antennas at each TX
%         mu:   back-off power factor. 0 < mu <= 1, allows to adapt to possible power
%               outages in case of instantaneous power constraints. (see paper)
%         P:    Power (or SNR) or the setting 
%         pow_control: Different power normalizations available. 
%                Choices are:
%
%                   - 'per_Antenna': Instantaneous power norm. per antenna 
%                   - 'per_TX': Instantaneous power norm. per TX 
%                   - 'average_perTX': Average power norm. per antenna      
%                   - 'average_perAntenna': Average power norm. per TX 
%
% OUTPUT: Vector of size Num.TX-antennas x Num.RXs 
%

    K  = size(He_M,1); % Number of RXs
    Mt = size(He_M,2); % Total number of transmit antennas
        
    W = zeros(Mt,K);    % Initializing Precoding matrix 
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Precoder at dumb TXs
    W_m = 1/sqrt(2*K*max(Mn(2:end)))*ones(Mt-Mn(1),K); % Just ones 

    if strcmp(pow_control, 'per_Antenna')
        normalization = min(max(1,log10(P)/2),sqrt(2));            
    elseif strcmp(pow_control, 'per_TX')
        normalization = min(max(1,log10(P)/2),sqrt(2));
    elseif strcmp(pow_control, 'average_perTX')
        normalization = min(max(1,log10(P)/2),sqrt(2))*sqrt(Mn(1));
    elseif strcmp(pow_control, 'average_perAntenna')
        normalization = min(max(1,log10(P)/2),sqrt(2))/2;        
        if K == 2 && Mn(1) ==1
            normalization = max(1,sqrt(log10(P)/2));
        end
    else 
        error('Unknown power control')
    end

    % Building AP-ZF precoder at passive TXs, i.e., TX 2... to TX M 
    W(Mn(1)+1:end, :) = mu*W_m/normalization; % Coef. used by TX m  
         
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Precoder at smart TXs (in the current state of the code, smart TXs = TX 1)
    He = squeeze(He_M(:,:,1)); % Channel estimate at TX 1
        
    % Assuming perfect estimate of other TXs' precoder (since it does not depend on
    % local CSIT at other TX's
    W_atTX1 = W;
        
    % Reducing error 
        W1_preNorm  = zeros(Mn(1),K); % Init precoder at TX 1  
        idx_ant_noTX1 = Mn(1)+1:Mt;     % Index of antennas for TXs 2 ... M
        idx_ant_TX1   = 1:Mn(1);        % Index of antennas for TX 1

        for i = 1:K   % For each user            
            H_i_bar         = He([1:i-1,i+1:end],:);   % Channel matrix for all the RXs except RX i
            h_i_bar_TX1     = H_i_bar(:, idx_ant_TX1);   % Coefs of TXs except TX 1
            h_i_bar_not_TX1 = H_i_bar(:, idx_ant_noTX1); % Coefs of TX 1

            h_inv = h_i_bar_TX1'/(h_i_bar_TX1*h_i_bar_TX1'+ 1/sqrt(P)*eye(K-1))*h_i_bar_not_TX1; 

            W1_preNorm(:,i) = -h_inv*W_atTX1(idx_ant_noTX1,i);    % AP-ZF Precoder at TX 1  
        end  

    %% Precoder normalization per TX  
    if strcmp(pow_control, 'per_Antenna')
        norm_perAntenna = sqrt(sum(abs(W1_preNorm).^2,2));
        normalization = max([norm_perAntenna; 1]);
        W(1:Mn(1), :) = W1_preNorm/normalization;     
        
    elseif strcmp(pow_control, 'per_TX')
        norm_perTX = norm(W1_preNorm);
        W(1:Mn(1), :)  = W1_preNorm/max([norm_perTX; 1]);  

    elseif strcmp(pow_control, 'average_perTX')
        W(1:Mn(1), :) = W1_preNorm;
        
    elseif strcmp(pow_control, 'average_perAntenna')
        W1 = W1_preNorm;
        W(1:Mn(1), :) = W1;
        
    else 
        error('Unknown power control')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Antonio Bazco-Nogueras  
% Date: 2022/05/11
% Contact: antonio.bazco@imdea.org
% License: This file can be distributed, remixed, adapted, and other work can be
% built upon it, as long as appropiate credit is included for the original creation. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%