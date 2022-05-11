function V_tot = precoder_SINGLE(H, Mn, P, pow_control)
% Precoder that transmits only from the best TX (TX 1 in the current state of the code). 
% In this case, only TX 1 transmits. The number of simultanoeusly server users
% is reduced to the dimensionality of the transmitter.
%
% INPUTS: H:    Channel matrix estimates at each TX
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
    He = squeeze(H(:,1:Mn(1),1));  % Channel estimate at TX m for antennas of TX 1
    
    K = size(H,1);   % Number of RXs
    Mt = size(H,2);  % Total number of transmit antennas
    M1 = Mn(1);      % Num. antennas at TX 1
        
    V = zeros(M1,K); % Initializng Precoding matrix
    
    n_RXs = K;  % Auxiliar variable for the number of served RXs
    
    %% Precoder for different RXs

    if M1 < K % Serving only as many users as possible
        He = He(1:M1,:); % submatrix of users of interest
        n_RXs = M1;      % Serving only the maximum possible amount (num. of TX antennas)
    end

    for i = 1:n_RXs % For RX i
        He_bar_i = He([1:i-1,i+1:end],:); % Other RXs' channels
        h_i      = He(i,:); % RX i's channel

        % Orthogonal projection matrix on the null space of other RXs with regularization
        P_oirt_i =  eye(M1)- He_bar_i'/(He_bar_i*He_bar_i'+ n_RXs/P)*He_bar_i;

        % Projection of matched filter onto the null space
        v_proj = P_oirt_i*h_i';

        % unit-norm vector normalization:   
        v = v_proj/norm(v_proj); 

        V(:,i) = v; % Vector for user i
    end
   
    %% Precoder normalization per antenna or TX    
    if strcmp(pow_control,'per_Antenna')
        V_norm = V/norm(V); 
    elseif strcmp(pow_control,'per_TX')
        V_norm = V/norm(V); 
    elseif strcmp(pow_control, 'average_perTX')
        V_norm = V/sqrt(n_RXs); 
    elseif strcmp(pow_control, 'average_perAntenna')
        V_norm = V*sqrt(Mn(1)/n_RXs);    
    else
         error('Unknown power control')
    end
    
    % Transmitting nothing to other RXs
    V_tot = [V_norm; zeros(Mt-M1,K)]; 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Antonio Bazco-Nogueras  
% Date: 2022/05/11
% Contact: antonio.bazco@imdea.org
% License: This file can be distributed, remixed, adapted, and other work can be
% built upon it, as long as appropiate credit is included for the original creation. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%