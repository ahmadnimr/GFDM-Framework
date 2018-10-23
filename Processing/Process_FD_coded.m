function [an, meas] = Process_FD_coded(params,rchan, N0,Nrun, Ncp,ChCode)
% Frequency domain processing 
% parmas: Precoder parameters,  see GenParameters.m
% rchan: Channel parameters, see Channel_functions.GenFadingChannel
% N0 : Noise variance (scalar)
% Nrun: number of realizations of channel (Noise and fading)
% Ncp: OFDM CP
% ChCode: code parameters, see Data_functions/GenerateCode
an.SNR = zeros(params.K,1);
an.SER = zeros(params.K,1);
% Used libraries
f_data = Data_functions();
f_channel = Channel_functions();
b_func = Block_functions();
meas = CreateMeasurments(); % measurement
N = params.N;
Nt = 2*params.N; % OFDM FFT number of subcarriers
xfu = zeros(Nt,1);
yf = zeros(params.N,1);
hf = zeros(params.N,1);

for nr=1:Nrun % different realization
    %% Encoder
    %-------------------------------------------------
    % Gen encoded data for one block
    info_bits = randi([0 1],ChCode.Nr_infobits,1);
    [ qam_data, encoded_bits] = f_data.Encoder(ChCode, info_bits);
    Tx.InfoBits =  info_bits; % info bits
    Tx.EncodedBits = encoded_bits; % encoded bits
    Tx.D = reshape(qam_data, params.K, params.M); % QAM symbols
    %-------------------------------------------------
    
    %% FD modulator (Precoder)
    %-------------------------------------------------
    xf = Precode(Tx.D, params); % 
    %-------------------------------------------------
    %% OFDM TX
    %-------------------------------------------------
    % Allocation to OFDM no DC
    xfu((1:N/2)+1) =  xf(1:N/2);
    xfu(end+1-(1:N/2)) =  xf(end+1-(1:N/2));
    % OFDM signal: Es per symbols is preserved
    %-------------------------------------------------
    xt = sqrt(Nt)*ifft(xfu);
    x = b_func.AddCpCs( xt, Nt, Ncp, 0); % add CP
    %% Channel
    %-------------------------------------------------
    [he, yr ] = f_channel.ApplyChannel( rchan, x, Ncp);
    yr = yr+f_channel.GenRandomNoise([Nt+Ncp,1],N0);
    %-------------------------------------------------
    
    %% OFDM Receiver
    %-------------------------------------------------
    yr =  b_func.RemoveCpCs( yr, Nt, Ncp); % remove CP
    yrf = 1/sqrt(Nt)*fft(yr);
    % Channel estimation
    he = he(Ncp+(1:Nt),:);
    hfr  = fft(he);
    % Deallocation
    yf(1:N/2) =  yrf((1:N/2)+1);
    yf(end+1-(1:N/2)) =  yrf(end+1-(1:N/2));
    hf(1:N/2) =  hfr((1:N/2)+1);
    hf(end+1-(1:N/2)) =  hfr(end+1-(1:N/2));
    %-------------------------------------------------
    %% Equalization and inverse precoding
    %-------------------------------------------------
    % Compute Equalization params depending on the Precoding
    eq = GenChEq(params,hf, 1, N0);
    % For anylitical measure
    an.SNR  = an.SNR + eq.SNR ;
    an.SER  = an.SER + eq.SER ;
    
    % Equalization and Inverse precoding
    Rx.D = Equalization(yf, eq, params); % QAM SISO constructed
    %-------------------------------------------------
    %% Decoder
    %-------------------------------------------------
    %QAM demod to compute SER
    Tx.Di = f_data.QAM_demod(Tx.D, params.Mc);
    Rx.Di = f_data.QAM_demod(Rx.D, params.Mc);
    %-------------------------------------------------
    ChCode.info_bits = info_bits;
    % Soft decoding
    [ bits_e, encoded_bits_e ] = f_data.Decoder( ChCode, Rx.D(:), eq.SNR(:));
    %% Simulation measurment
    Rx.InfoBits =  bits_e;
    Rx.EncodedBits = encoded_bits_e;
    meas =  meas.Measure(meas,Tx, Rx);
    
end % end realization loop
%% Averaging over number of realizations
meas = meas.Avg(meas);
an.SNR  = an.SNR/Nrun;
an.SER  = an.SER/Nrun;
end
