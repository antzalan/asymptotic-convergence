function V = precoder_NAIVE(He_M, Mn, P, pow_cntrl)
% Precoder that applies a naive precoder. 
% In this case, each TX computes the standard centralized precoder assuming that 
% all other TXs *share* the same CSI as itself. However, in reality, each TX has 
% a different estimate
%
% INPUTS: He_M: Channel matrix estimates at each TX
%         Mn:   Num. antennas at each TX
%         P:    Power (or SNR) or the setting 
%         pow_cntrl: Different power normalizations available. 
%                Choices are:
%
%                   - 'per_Antenna': Instantaneous power norm. per antenna 
%                   - 'per_TX': Instantaneous power norm. per TX 
%                   - 'average_perTX': Average power norm. per antenna      
%                   - 'average_perAntenna': Average power norm. per TX 
%
% OUTPUT: Vector of size Num.TX-antennas x Num.RXs 
%
    K    = size(He_M,1); % Number of RXs
    Mt   = size(He_M,2); % Total number of transmit antennas
    M    = size(He_M,3); % Number of TXs
    
    V = zeros(Mt,K);     % Initializng Precoding matrix
    
    n_ant_prev = 0;      % Initializing index of transmit antenna
    
    for m = 1:M % At TX m        
        He_m = He_M(:,:,m); % Channel estimate at TX m

        V_m = precoder_CENTRALIZED(He_m, Mn, P, pow_cntrl); % Naive precoder at TX m

        % Antenna indices for TX m
        vec_ant_TX_m = n_ant_prev + 1 : n_ant_prev + Mn(m); 
        n_ant_prev = n_ant_prev + Mn(m); % Updating antennas idx for next TX
                
        %%% Effective precoder: Each TX transmit the precoder obtained from its own
        %%% estimate (important). Hence, inconsistencies may appear
        V(vec_ant_TX_m,:) = V_m(vec_ant_TX_m,:);

    end     
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Antonio Bazco-Nogueras  
% Date: 2022/05/11
% Contact: antonio.bazco@imdea.org
% License: This file can be distributed, remixed, adapted, and other work can be
% built upon it, as long as appropiate credit is included for the original creation. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%