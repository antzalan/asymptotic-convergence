function V_norm = precoder_CENTRALIZED(He, Mn, P, pow_control)
% precoder_CENTRALIZED:
% 
% Standard centralized precoder. The same CSI is shared by all the TXs, and the 
% decision is taken in acentralized manner. 
%
% INPUTS: He_M: Channel matrix estimates at each TX
%         Mn:   Num. antennas at each TX
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
    K  = size(He,1); % Number of RXs
    Mt = size(He,2); % Total number of transmit antennas
    M  = length(Mn); % Number of TXs
        
    V  = zeros(Mt,K);  % Initializing Precoding matrix
    
    %% Precoder for each RX
    for i = 1:K % For RX i
        He_bar_i = He([1:i-1,i+1:end],:); % Other RXs' channels        
        h_i      = He(i,:); % RX i's channel

        % Orthogonal projection matrix on the null space of other RXs with regularization
        P_oirt_i =  eye(Mt)- He_bar_i'/(He_bar_i*He_bar_i'+ K/P)*He_bar_i; 

        % Projection of matched filter onto the null space
        v_proj = P_oirt_i*h_i';

        % unit-norm vector normalization:   
        v = v_proj/norm(v_proj); 

        V(:,i) = v; % Vector for user i
    end

    %% Precoder normalization 
    if strcmp(pow_control, 'per_Antenna')
        norm_perAntenna = sqrt(sum(abs(V).^2,2));
        V_norm          = V/max([norm_perAntenna; 1]);    

    elseif strcmp(pow_control, 'per_TX')
        n_ant_prev = 0;          % Initializing index of transmit antenna
        pow_TX     = zeros(1,M); % Initializing power at each TX
        for m = 1:M  % Compute power at TX m  
            pow_TX(m) = norm(V(n_ant_prev + 1 : n_ant_prev + Mn(m),  :));
            n_ant_prev = n_ant_prev + Mn(m);% Updating index to get next TX's antennas
        end
        max_pow_TX = max([pow_TX, 1]); % Highest power comsumption
        V_norm     = V/max_pow_TX;     % Normalization

    elseif strcmp(pow_control, 'average_perTX')
        V_norm = V*sqrt(Mt/(K*max(Mn)));         
    elseif strcmp(pow_control, 'average_perAntenna')
        V_norm = V*sqrt(Mt/K);         
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
