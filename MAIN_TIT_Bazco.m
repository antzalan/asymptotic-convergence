%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Code for the publication: 
%     "Asymptotically Achieving Centralized Rate on the Decentralized 
%       Network MISO Channel", IEEE TRANSACTIONS ON INFORMATION THEORY, VOL. 68, 
%       NO. 1, JANUARY 2022, A. Bazco-Nogueras, P. de Kerret, D. Gesbert, and N.
%       Gresset.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this script, we can find the code for simulating the scenario anlyzed in the
% above work. The parameters follow the same notation and there exist some variables
% that allow for several settings and configurations. 
% Any doubt or question, please contact: antonio.bazco@imdea.org
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
clear variables, rng(1);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% System size configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    M = 2;      % Number of TXs
    K = 2;      % Number of RXs

    Mn = [1, 1*ones(1,M-1)];  % Number of antennas at each TX
    Kn = ones(1,K);           % Number of antennas at each RX

    Mt = sum(Mn);  % Total number of TX antennas
    Kt = sum(Kn);  % Total number of RX antennas
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% System power configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    P_nom_dB = [linspace(5,80,9)];   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CSIT accuracy configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Definition of the vector alpha values.
    %%% Now, we have two different levels of accuracy, alpha_1 and alpha_2.
    %%% The first TX has accuracay alpha_1, the other TXs have accuracy alpha_2
    %%% (Code can be modified to allow other topologies). 

    alpha_1   = 1;    % Scaling of CSIT accuracy 
    alpha_2   = 0.6;  % Scaling of CSIT accuracy 
    
    alphas = [alpha_1*ones(1,1),  alpha_2*ones(1, M-1)];

    %%% Create matrix with CSI accuracy for each channel coef. at each TX
    csi_TXs = shiftdim(repmat(alphas,[Mt, 1,Kt]),2); 

    %%%%%   %%%%%%%%%%   %%%%%%%%%%   %%%%%%%%%%   %%%%%%%%%%   %%%%%
    
    %%% Possible values of mu (back-off factor for power. see paper)
    mu_v_apzf = linspace(0.1, 1, 20); 
    mu_v_hdzf = linspace(0.8, 1, 10); 
    mu_v_cdzf = linspace(0.7, 1, 31);    
    mu        = 1; % In case only one value is selected
    
    search_apzf_v = 1;  %% Boolean. 0 if mu is fixed for AP-ZF. 
                        %% 1 if search over possible values of mu in mu_v_apzf is applied
    search_cdzf_v = 1;  %% Boolean. 0 if mu is fixed for CD-ZF. 
                        %% 1 if search over possible values of mu in mu_v_cdzf is applied

    %%%%%   %%%%%%%%%%   %%%%%%%%%%   %%%%%%%%%%   %%%%%%%%%%   %%%%%

    %%% QUANTIZATION OF CSIT FOR PROPOSED ALGORITHM                   
    q_v       = [linspace(1e-3, .9, 1e2)]; % Quantization step 
    
    search_q_v    = 1;  %% Boolean. 0 if quantization step is fixed. 
                        %% 1 if search over possible values of quantization step is applied

    %%% Select if channel quantization or precoder quantization. (It is Better to 
    %%% choose precoder quantization, better performance, reduced dimensionality. see paper)
    %%% quant_channel == 1 => the quantization is done at the input (channel matrix)
    %%% quant_channel == 0 => the quantization is done at the output (precoder vector)
    quant_channel = 0; 

    % Values for quantization step. If quantization step is fixed, the quantization
    % step is given by q_step = d*sqrt(P_nom^(-alpha_2/k))
    d = .7;  k = 2;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Type of power normalization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Uncomment the wanted case.
%%%  'per_Antenna':        Instantaneous power norm. per antenna
%%%  'per_TX':             Instantaneous power norm. per TX
%%%  'average_perAntenna': Average power norm. per antenna
%%%  'average_perTX':      Average power norm. per TX    

    pow_control = 'per_Antenna';
%     pow_control = 'per_TX';
%     pow_control = 'average_perTX';    
%     pow_control = 'average_perAntenna';
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Montecarlo configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    Num_channels_iter   = 1*1e3;   % Iterations over channel
    
    period_display = Num_channels_iter/10;  % Display computation time every "period_display"
                                            % iterations of codebook  

    %%% Name of the file where data is saved
    file_name = ['test_titBazco',date];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialization of vector and matrix variables
    N_P     = length(P_nom_dB);        % Number of different SNR to test
    P_nom_v = db2pow(P_nom_dB);        % SNR in linear scale 
    logP_v  = log2(db2pow(P_nom_dB));  % SNR in logarithm-2 scale
      
    R_tdma_1        = zeros(N_P,  Num_channels_iter, K); 
    R_tdma_2        = zeros(N_P,  Num_channels_iter, K); 
    R_naive         = zeros(N_P,  Num_channels_iter, K); 
    R_sing          = zeros(N_P,  Num_channels_iter, K); 
    R_centr         = zeros(N_P,  Num_channels_iter, K); 
    R_perf          = zeros(N_P,  Num_channels_iter, K); 
    
    P_tdma          = zeros(N_P, Num_channels_iter, Mt);    
    P_naiv          = zeros(N_P, Num_channels_iter, Mt, K);    
    P_sing          = zeros(N_P, Num_channels_iter, Mt, K);           
    P_cent          = zeros(N_P, Num_channels_iter, Mt, K);    
    P_perf          = zeros(N_P, Num_channels_iter, Mt, K);    

    if search_q_v == 1
        R_cdzf    = zeros(N_P,  Num_channels_iter, length(q_v),  length(mu_v_cdzf), K); 
        agree_v   = zeros(N_P,  Num_channels_iter, length(q_v),  length(mu_v_cdzf), Mt-Mn(1), Kt); 
        P_cdzf    = zeros(N_P,  Num_channels_iter, length(q_v),  length(mu_v_cdzf), Mt, K);    
    else 
        R_cdzf    = zeros(N_P,  Num_channels_iter, K);
        P_cdzf    = zeros(N_P, Num_channels_iter, Mt, K); 
        agree_v   = zeros(N_P, Num_channels_iter, Mt-Mn(1), Kt);
    end
    
    if search_cdzf_v == 1
        R_hdzf    = zeros(N_P,  Num_channels_iter, length(mu_v_hdzf), K ); 
        P_hdzf    = zeros(N_P,  Num_channels_iter, length(mu_v_hdzf), Mt, K);  
    else
        R_hdzf    = zeros(N_P,  Num_channels_iter, K );         
        P_hdzf    = zeros(N_P, Num_channels_iter, Mt, K);    
    end
    
    if search_apzf_v == 1
        R_apzf    = zeros(N_P,  Num_channels_iter, length(mu_v_apzf), Kt ); 
        P_apzf    = zeros(N_P,  Num_channels_iter, length(mu_v_apzf), Mt, K);    
    else
        R_apzf    = zeros(N_P,  Num_channels_iter, K );
        P_apzf    = zeros(N_P, Num_channels_iter, Mt, K);    
    end    

    t_sum = 0;  % Initializing variable for time
              
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Loop for the power values
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     for n_p = 1:length(P_nom_dB)       
        
        P_nom = db2pow(P_nom_dB(n_p));  % SNR in linear scale 
         
        %%% Displaying some information in the command window
        disp(' ');  disp('--------------');        
        disp(['Power num: ',num2str(n_p), ' of ', num2str(N_P) ,'.// equal to: ',num2str(P_nom_dB(n_p)),' dB = ', num2str(P_nom)])
        disp('--------------');      disp(' ')  ;      

        %%% The following three lines are related to CSI model and equation (13) in the paper
        sigma2_TXs = P_nom.^(-csi_TXs); % Variance of the CSI error obtaiend from alpha values
        z_M        = sqrt(sigma2_TXs);  % Standard deviation of CSI error for each TX
        zop_M      = sqrt(1 - sigma2_TXs); % Power normalization of correct channel value to obtain the estimation.
        
        %%% Creating the channel and random gaussian additive error values for all iterations.
        H_tot     = 1/sqrt(2)*(randn(Kt, Mt, Num_channels_iter)    + 1i*randn(Kt ,Mt, Num_channels_iter   )); % channel
        error_tot = 1/sqrt(2)*(randn(Kt, Mt, Num_channels_iter, M) + 1i*randn(Kt, Mt, Num_channels_iter, M)); % error, one value per TX

        %%% In case the quantization step is fixed, we take this value
        %%% (i.e., in case search_q = false)
        q_step = d*sqrt(P_nom^(-alpha_2/k));        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Loop for the channel realizations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        for n_h = 1:Num_channels_iter       
        tic,
            %%% Generate the channel matrix and the CSIT matrices
            H     = squeeze(H_tot(:,:,n_h));
            error = squeeze(error_tot(:,:,n_h,:));
            H_M   = zop_M.*repmat(H,[1,1,M]) + z_M.*error;  % A different CSI estimate for each TX 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Precoding Schemes (From Most Naive to Perfect)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   

            %%%%%%%%%%%%%%%%%%%%%%%%%
            %% TDMA Precoding (beamforming) Uses matched filter precoder
            %%%%%%%%%%%%%%%%%%%%%%%%%      
            %%% We only compute the rate for RX 1, as the performance is the same 
            %%% for all RXs due to simmetry. 
            prec_rx1_tdma_1 = precoder_TDMA_rx1(H_M, Mn, 1, pow_control);
            Rate_rx1        = log2(1 + P_nom*abs(H(1,:)*prec_rx1_tdma_1).^2); 
            R_tdma_1(n_p, n_h, 1) = Rate_rx1;        % Saving the rate for this experiment
            P_tdma(n_p, n_h, :)   = prec_rx1_tdma_1; % Saving the precoder for this experiment
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %% TDMA Precoding (no beamforming) % Commented because previous case performs better
            %%%%%%%%%%%%%%%%%%%%%%%%%
%             v_rx1_tdma_2 = precoder_TDMA_rx1(H_M, Mn, 0, inst_pow_control);
%             R_rx1  = log2(1 + P_nom*abs(H(1,:)*v_rx1_tdma_2).^2); % Symmetry for other users.
%             R_tdma_2(n_p, n_h, 1) = R_rx1;
%             P_tdma2(n_p, n_h, :,:) = v_rx1_tdma_2;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %% Naive I Precoding
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Each TX believes that the other TXs have access to the same CSI as itself
            prec_naive = precoder_NAIVE(H_M, Mn, P_nom, pow_control);
            R_naive(n_p, n_h,:)   = rate_user(H, prec_naive, P_nom); % Saving the rates for this experiment
            P_naiv(n_p, n_h, :,:) = prec_naive; % Saving the precoders for this experiment
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %% Single TX Precoding
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Only the best TX (TX 1) transmits
            V_sing = precoder_SINGLE(H_M, Mn, P_nom, pow_control);
            R_sing(n_p, n_h,:)    = rate_user(H, V_sing, P_nom); % Saving the rates for this experiment           
            P_sing(n_p, n_h, :,:) = V_sing; % Saving the precoders for this experiment 
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %% AP-ZF 
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Scheme from work: "Degrees of freedom of the network MIMO channel with distributed CSI,"
            %%% P. D. Kerret and D. Gesbert, IEEE Trans. Inf. Theory, vol. 58, no. 11, pp. 6806–6824, Nov. 2012. 

            %%% If search_apzf_v == 1, we evaluate the algorithm for a set of
            %%% different values of mu_apzf (power back-off factor). 
            %%% Otherwise, this value is fixed. 
            if search_apzf_v == 1
                for i_mu = 1:length(mu_v_apzf) 
                    mu_apzf = mu_v_apzf(i_mu);

                    V_apzf = precoder_APZF(H_M, Mn, mu_apzf, P_nom, pow_control);                    
                    R_apzf(n_p, n_h, i_mu, :)   = rate_user(H, V_apzf, P_nom);  % Saving the rates for this experiment
                    P_apzf(n_p, n_h, i_mu, :,:) = V_apzf; % Saving the precoders for this experiment     
                end
            else
                V_apzf = precoder_APZF(H_M, Mn, mu, P_nom, pow_control);
                R_apzf(n_p, n_h,:)  = rate_user(H, V_apzf, P_nom);  % Saving the rates for this experiment
                P_apzf(n_p, n_h, :,:) = V_apzf; % Saving the precoders for this experiment     
            end
                      
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %% CD-ZF
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %%% If search_q_v == 1, we evaluate the algorithm for a set of
            %%% different values of mu_apzf. Otherwise, this value is fixed.           
            if search_q_v == 1
                for i_q = 1:length(q_v) % for each quantization step considered
                    q_step = q_v(i_q);
                    %%% Since the quantization applied is quantizing the closest value towards 0, the
                    %%% quantization applied to the channel matrix reduces the mean power of the
                    %%% coefficients proportionally to the following expression:
                    quant_pw_reduct_ch = (1-cdf('Exponential', q_step, 1)); 

                    for i_mu = 1:length(mu_v_cdzf) 
                        mu_cdzf = mu_v_cdzf(i_mu);

                        [V_cdzf, agree] = precoder_CDZF(H_M, Mn, q_step, mu_cdzf, quant_channel, P_nom, pow_control, quant_pw_reduct_ch);
                        R_cdzf(n_p, n_h, i_q, i_mu, :)    = rate_user(H, V_cdzf, P_nom); % Saving the rates for this experiment
                        P_cdzf(n_p, n_h, i_q, i_mu, :, :) = V_cdzf;  % Saving the precoders for this experiment     
                    end
                end
            else
                quant_pw_reduct_ch = (1-cdf('Exponential', q_step, 1)); 
                [V_cdzf, agree] = precoder_CDZF(H_M, Mn, q_step, mu, quant_channel, P_nom, pow_control, quant_pw_reduct_ch);
                R_cdzf(n_p, n_h, :)    = rate_user(H, V_cdzf, P_nom); % Saving the rates for this experiment
                P_cdzf(n_p, n_h, :, :) = V_cdzf;  % Saving the precoders for this experiment                      
            end

            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %% CD-ZF Hierarchical
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %%% In this case, the CSIT is hierarchical, meaning that the more
            %%% accurate TXs have access to the CSIT of the less accurate TXs. This
            %%% is valid for example in cases where the the same information is
            %%% sent with multiple resolution when the link are heterogeneous. 


            %%% If search_cdzf_v == 1, we evaluate the algorithm for a set of
            %%% different values of mu_hdzf (power back-off factor). 
            %%% Otherwise, this value is fixed. 
            if search_cdzf_v == 1
                for i_mu = 1:length(mu_v_hdzf) 
                    mu_hdzf = mu_v_hdzf(i_mu);

                    V_cdzf_h2 = precoder_CDZF_HIER(H_M, Mn, mu_hdzf, P_nom, pow_control);
                    R_hdzf(n_p, n_h, i_mu, :)   = rate_user(H, V_cdzf_h2, P_nom); % Saving the rates for this experiment
                    P_hdzf(n_p, n_h, i_mu, :,:) = V_cdzf_h2; % Saving the precoders for this experiment        
                end
            else
                V_cdzf_h2 = precoder_CDZF_HIER(H_M, Mn, mu, P_nom, pow_control);
                R_hdzf(n_p, n_h,:)  = rate_user(H, V_cdzf_h2, P_nom); % Saving the rates for this experiment           
                P_hdzf(n_p, n_h, :,:) = V_cdzf_h2; % Saving the precoders for this experiment                       
            end                         

            %%%%%%%%%%%%%%%%%%%%%%%%%
            %% Centralized CSI Precoding
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Standard precoding where all TXs share the same CSI. The CSI used
            %%% is the one of the best TX
            V_centr = precoder_CENTRALIZED(squeeze(H_M(:,:,1)), Mn, P_nom, pow_control);
            R_centr(n_p, n_h,:) = rate_user(H, V_centr, P_nom); % Saving the rates for this experiment      
            P_cent(n_p, n_h, :,:) = V_centr; % Saving the precoders for this experiment              
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %% Perfect CSI Precoding
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Case where the CSI is perfect (actual channel used as CSI)
            V_perf = precoder_CENTRALIZED(H, Mn, P_nom, pow_control);
            R_perf(n_p, n_h,:)  = rate_user(H, V_perf, P_nom); % Saving the rates for this experiment  
            P_perf(n_p, n_h, :,:) = V_perf; % Saving the precoders for this experiment                      
                
            
            %%%%%%%%%%%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%
            %%  Computing time of simulation... 
            ttt = toc;
            if mod(n_h, period_display) == 0 % time to display the passed time
                disp(['time last ',num2str(period_display) ,' channel uses: ', num2str(period_display*ttt)])
            end   
            t_sum = t_sum + ttt;
            
        end
    end        
     
    % Displaying total time passed for the simulation
    disp(' ');  
    disp('--------------');       
    disp('--------------  --------------  --------------  --------------');       
    disp(['Total Time: ',num2str(t_sum), ' segs// or ', num2str(t_sum/60) ,' mins'])
    disp('--------------  --------------  --------------  --------------');       
    disp('--------------');
    disp(' ')  ; 
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Averaging the saved values to obtain the metrics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Expected rate per user
    ExpRate_rxi_tdma_1  = squeeze(mean(R_tdma_1,2));
    ExpRate_rxi_tdma_2  = squeeze(mean(R_tdma_2,2));
    ExpRate_rxi_naive   = squeeze(mean(R_naive, 2));
    ExpRate_rxi_sing    = squeeze(mean(R_sing,  2));
    ExpRate_rxi_apzf    = squeeze(mean(R_apzf,  2));
    ExpRate_rxi_cdzf    = squeeze(mean(R_cdzf,  2));
    ExpRate_rxi_hdzf    = squeeze(mean(R_hdzf,  2));
    ExpRate_rxi_centr   = squeeze(mean(R_centr, 2));
    ExpRate_rxi_perf    = squeeze(mean(R_perf,  2));
    
    %% Expected rate mean
    ExpRate_tdma_1  = squeeze(mean(ExpRate_rxi_tdma_1, 2));
    ExpRate_tdma_2  = squeeze(mean(ExpRate_rxi_tdma_2, 2));
    ExpRate_naive   = squeeze(mean(ExpRate_rxi_naive,  2));
    ExpRate_sing    = squeeze(mean(ExpRate_rxi_sing,   2));
    ExpRate_centr   = squeeze(mean(ExpRate_rxi_centr,  2));
    ExpRate_perf    = squeeze(mean(ExpRate_rxi_perf,   2));
    
    % Expected rate mean for APZF 
    if search_apzf_v == 1                
        ExpRate_apzf_t = squeeze(mean(squeeze(ExpRate_rxi_apzf),  3));
        [max_mu_q_apzf, pos_max_mu_q_apzf] = max(ExpRate_apzf_t,[],[2],'linear');
        [I1_apzf,I2_apzf] = ind2sub([N_P, length(mu_v_apzf)], squeeze(pos_max_mu_q_apzf));

        ExpRate_apzf = max_mu_q_apzf;        
    else
        ExpRate_apzf    = squeeze(mean(ExpRate_rxi_apzf,  2));    
    end
    
    % Expected rate mean for CDZF
    if search_q_v == 1
        ExpRate_cdzf_t = squeeze(mean(squeeze(ExpRate_rxi_cdzf),  4));
        [max_mu_q_cdzf, pos_max_mu_q_cdzf] = max(ExpRate_cdzf_t,[],[2,3],'linear');
        [~,I2,I3,~] = ind2sub([N_P,length(q_v),length(mu_v_cdzf)], squeeze(pos_max_mu_q_cdzf));

        ExpRate_cdzf = max_mu_q_cdzf;    
    else
        ExpRate_cdzf    = squeeze(mean(ExpRate_rxi_cdzf,  2));
    end
    
    % Expected rate mean for CDZF Hierarchical
    if search_cdzf_v == 1
        ExpRate_hdzf_t = squeeze(mean(squeeze(ExpRate_rxi_hdzf),  3));
        [max_mu_q_hdzf, pos_max_mu_q_hdzf] = max(ExpRate_hdzf_t,[],[2],'linear');
        [I1_hdzf,I2_hdzf] = ind2sub([N_P, length(mu_v_hdzf)], squeeze(pos_max_mu_q_hdzf));

        ExpRate_hdzf = max_mu_q_hdzf;    
    else
        ExpRate_hdzf    = squeeze(mean(ExpRate_rxi_hdzf,  2));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %% Plotting
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%% Plotting expected sum rate
        figure, hold on, grid on
%         plot(P_nom_dB, ExpRate_perf, 'k-')        
        plot(P_nom_dB, ExpRate_centr,   'k--')   
        plot(P_nom_dB, ExpRate_hdzf, '-v', 'Color', '#0072BD')      
        plot(P_nom_dB, ExpRate_cdzf,    '-o', 'Color', '#0072BD')      
        plot(P_nom_dB, ExpRate_apzf,    '-s', 'Color', '#4DBEEE')      
        plot(P_nom_dB, ExpRate_naive,   '-.', 'Color', '#D95319') 
        plot(P_nom_dB, ExpRate_sing,    '-.', 'Color', '#A2142F') 
%         plot(P_nom_dB, ExpRate_tdma_1,  '-', 'Color', '#77AC30')     

        legend(... % 'Perfect CSIT', ...
            'Centralized CSIT', ...
            'CD-ZF Hierar', ...
            'CD-ZF', ...
            'AP-ZF', ...
            'Naive ZF', ...
            'Single', ... % 'TDMA'
            'location','north west')     
            
        xlabel('P [dB]')
        ylabel('Rate [bits/s/Hz]')
        title('Average Rate')      
               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
%% Saving Data...
%     time_print = clock;
%     time_str = [num2str(time_print(1)),'_',num2str(time_print(2)),'_',num2str(time_print(3)),'_Time=_',num2str(time_print(4)),'h',num2str(time_print(5))];file_name = sprintf(time_str);
%     file_name = ['simulation_k=3_nh=500_nc=100_',date];
%     save(['saved_files/' file_name '.mat'])    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%------------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Antonio Bazco-Nogueras  
% Date: 2022/05/11
% Contact: antonio.bazco@imdea.org
% License: This file can be distributed, remixed, adapted, and other work can be
% built upon it, as long as appropiate credit is included for the original creation. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%