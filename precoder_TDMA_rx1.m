function V = precoder_TDMA_rx1(He_M, Mn, beam, pow_control)
% Precoder that applies TDMA (serves only one user simultaneously. 
% This code computes the precoder for RX 1, while the other users' precoder 
% follows from simmetry.
%
% INPUTS: He_M: Channel matrix estimates at each TX
%         Mn:   Num. antennas at each TX
%         beam: Whether the transmitted signal is beamfored towards the intended 
%               user or not. 
%         pow_control: Different power normalizations available. 
%                Choices are:
%
%                   - 'per_Antenna': Instantaneous power norm. per antenna 
%                   - 'per_TX': Instantaneous power norm. per TX 
%                   - 'average_perTX': Average power norm. per antenna      
%                   - 'average_perAntenna': Average power norm. per TX 
%
% OUTPUT: Vector of size 1 x sum(Mn) with the precoder for RX 1 from all TXs
%

    Mt   = size(He_M,2); % Total number of transmit antennas
    M    = size(He_M,3); % Number of TXs

    V = zeros(Mt,1);     % Initializng Precoding matrix

    n_ant_prev = 0;      % Initializing index of transmit antenna
    
    if beam == 1 %Apply beamforming towards intended user
        for m = 1:M % At TX m:

            He_m     = He_M(:,:,m); % Channel estimate at TX m
            He_m_RX1 = He_m(1,:);   % Channel estimate for RX 1

            V_m = He_m_RX1'; % Precoder is matched filter (mean(|V_m|^2) = num. antennas)
            
            %%% Now we apply the power normalization
            if strcmp(pow_control, 'per_Antenna') 
                V_m = V_m/max(abs(He_m_RX1));
            elseif strcmp(pow_control, 'average_perAntenna')
                V_m = V_m;
            elseif strcmp(pow_control, 'per_TX') 
                %%% We obtain the power required at each TX (*based on the estimation
                %%% of TX m*, very impportant). And then we normalize to ensure that
                %%% all TXs satisfy the constraint. 
                n_ant_prev_mm = 0;   % Initializing index of transmit antenna
                pow_TX = zeros(1,M); % Initializing power at each TX
                for mm = 1:M % Compute power at TX mm based on CSIT of TX m 
                    pow_TX(mm) = norm(V_m(n_ant_prev_mm + 1 : n_ant_prev_mm + Mn(mm),  :));
                    n_ant_prev_mm = n_ant_prev_mm + Mn(mm); % Updating index to get next TX's antennas
                end
                V_m = V_m/max(pow_TX);
            elseif strcmp(pow_control, 'average_perTX')
                V_m = V_m/sqrt(max(Mn)); % Normalization must be the same at each TX. Power limited by TX with more antennas. 
            else
                error('Unknown power control')
            end       
            
        %%%%   %%%%  %%%%   %%%%   %%%%   %%%%   %%%%   %%%%   %%%%   %%%%
        %%% Effective precoder: Each TX transmit the precoder obtained from its own
        %%% estimate (important). Hence, inconsistencies may appear
            vec_ant_TX_m = n_ant_prev + 1 : n_ant_prev + Mn(m); % Antenna indices at TX m
            V(vec_ant_TX_m) = V_m(vec_ant_TX_m); % Assigning final precoder
            n_ant_prev = n_ant_prev + Mn(m); % Updating previous antennas
        end

    else % No matched filter precoding, just ones        
        if strcmp(pow_control, 'per_Antenna') || strcmp(pow_control, 'average_perAntenna')
            V = ones(Mt,1);
        elseif strcmp(pow_control, 'per_TX') || strcmp(pow_control, 'average_perTX')
            n_ant_prev_mm = 0;
            for mm = 1:M %At TX m
                precoder_m = ones(Mn(mm),1)/sqrt(Mn(mm)); % Precoder at TX m
                vec_ant_TX_m = n_ant_prev_mm + 1 : n_ant_prev_mm + Mn(mm);  % Antenna indices at TX m
                V(vec_ant_TX_m) = precoder_m;
                n_ant_prev_mm = n_ant_prev_mm + Mn(mm); % Updating num. of antennas already managed
            end 
        else
            error('Unknown power control')
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Antonio Bazco-Nogueras  
% Date: 2022/05/11
% Contact: antonio.bazco@imdea.org
% License: This file can be distributed, remixed, adapted, and other work can be
% built upon it, as long as appropiate credit is included for the original creation. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%