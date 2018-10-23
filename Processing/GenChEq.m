function eq = GenChEq(params,hf, Es, N0)
% Channel equalization and precoding
% eq.SNR: scalar SNR per symbols
% eq.AG: gain per symbols
% eq.Heq: channel equlization
u_func = Utility_functions();
rx_func = Precoding_Rx_functions();
data_func = Data_functions();
Vh = u_func.V(hf,params.K,params.M);
if strcmp(params.rxType, 'pass')% OFDM
    H_eq = ones(params.K,params.M);
else % others
    H_eq = rx_func.CH_diag(Vh, params.rxType, Es, N0);
    eq.H_eq = H_eq;
end

switch params.Type
    case 'OFDM'
        [SNR, Ag] = rx_func.OFDM_SNR(Vh, Es, N0);
    case 'GFDM-Chirp'
        [SNR, Ag] =  rx_func.GFDM_SNR(Vh, H_eq, params.Atx, Es, N0);
        Ag = repmat(Ag, 1,params.M);
        SNR = repmat(SNR, 1,params.M);
    case 'SFFT'
        [SNR, Ag] =  rx_func.F3_off_SNR(Vh, H_eq, params.wtx, params.wrx, 1/params.K*Es, N0);
        Ag = repmat(Ag, params.K,params.M);
        SNR = repmat(SNR, params.K,params.M);
    case {'DFT-S-OFDM', 'SC'}
        [SNR, Ag] =  rx_func.DFT_S_FDM_SNR(Vh, H_eq, params.wtx, params.wrx, Es, N0);
        Ag = repmat(Ag, 1,params.M);
        SNR = repmat(SNR, 1,params.M);
    case '2D-DFT'
        [SNR, Ag] =  rx_func.F2_off_SNR(Vh, H_eq, params.wtx, params.wrx, params.K*Es, N0);
        Ag = repmat(Ag, params.K,params.M);
        SNR = repmat(SNR, params.K,params.M);
    otherwise
        eq.Ag = 0;
        eq.SNR = 0;
        eq.SER = 0;
end
eq.Ag = Ag;
eq.SNR = SNR;
eq.SER = data_func.QAM_SER_an(params.Mc, eq.SNR);
end