% u_func =  Utility_functions;
% K = 16; M = 4;
% N = K*M;
% U_KM1 = u_func.U(K,M);
% U_MK2 = u_func.U(M,K)';
% U_MK3 = u_func.U(M,K);
% %U_MK3 = eye(N);
% PI_KM = u_func.PI(K,M);
% PI_MK = u_func.PI(M,K);
% Ltx = exp(1j*2*pi*rand(N,1));
% wtx = 1/sqrt(M)*reshape(Ltx, K,M);
% Atx = 1/K*fft(wtx);
% wrx = 1./wtx;
% Af = PI_KM*U_MK3*diag(Ltx)*U_MK2*PI_MK*U_KM1*PI_KM;
% hf = randn(N,1)+1j*randn(N,1);
% N0 = 0.01;
% Es = 1;
% H = repmat(hf,1,N).*Af;
% pre_func = Precoding_Rx_functions();
% rxType = 'MMSE';
% Vh = u_func.V(hf, K,M);
% H_eq = pre_func.CH_diag(Vh, rxType, Es, N0);
% WH = pre_func.General(H, rxType, Es, N0);
% [SNR, LG, sig, I] = pre_func.General_SNR(WH,H, Es, N0);
% LG = reshape(LG, K,M);
% SNR = reshape(SNR, K,M);
% [SNR1, LG1] = pre_func.GFDM_SNR(Vh, H_eq, Atx, Es, N0);
%% DFT-S-OFDM
% u_func =  Utility_functions;
% K = 1; M = 64;
% N = K*M;
% U_KM1 = u_func.U(K,M);
% U_MK2 = u_func.U(M,K)';
% U_MK2 = eye(N);
% U_MK3 = u_func.U(M,K);
% U_MK3 = eye(N);
% PI_KM = u_func.PI(K,M);
% PI_MK = u_func.PI(M,K);
% Ltx = randn(N,1).*exp(1j*2*pi*rand(N,1));
% wtx = 1/sqrt(M)*reshape(Ltx, K,M);
% Atx = 1/K*fft(wtx);
% wrx = 1./wtx;
% Af = PI_KM*U_MK3*diag(Ltx)*U_MK2*PI_MK*U_KM1*PI_KM;
% hf = randn(N,1)+1j*randn(N,1);
% N0 = 0.01;
% Es = 1;
% H = repmat(hf,1,N).*Af;
% pre_func = Precoding_Rx_functions();
% rxType = 'MMSE';
% Vh = u_func.V(hf, K,M);
% H_eq = pre_func.CH_diag(Vh, rxType, Es, N0);
% WH = pre_func.General(H, rxType, Es, N0);
% [SNR, LG, sig, I] = pre_func.General_SNR(WH,H, Es, N0);
% LG = reshape(LG, K,M);
% SNR = reshape(SNR, K,M);
% [SNR1, LG1] = pre_func.DFT_S_FDM_SNR(Vh, H_eq, wtx,wrx, Es, N0);


 %% OTFS
u_func =  Utility_functions;
K = 16; M = 4;
N = K*M;
U_KM1 = u_func.U(K,M);
U_MK2 = u_func.U(M,K)';
%U_MK3 = u_func.U(M,K);
U_MK3 = eye(N);
PI_KM = u_func.PI(K,M);
PI_MK = u_func.PI(M,K);
Ltx = rand(N,1).*exp(1j*2*pi*rand(N,1));
wtx = 1/sqrt(M)*reshape(Ltx, K,M);
Atx = 1/K*fft(wtx);
wrx = 1./wtx;
Af = PI_KM*U_MK3*diag(Ltx)*U_MK2*PI_MK*U_KM1*PI_KM;
hf = randn(N,1)+1j*randn(N,1);
N0 = 0.01;
Es = 1;
H = repmat(hf,1,N).*Af;
pre_func = Precoding_Rx_functions();
rxType = 'ZF';
Vh = u_func.V(hf, K,M);
H_eq = pre_func.CH_diag(Vh, rxType, Es, N0);
WH = pre_func.General(H, rxType, Es, N0);
[SNR, LG, sig, I] = pre_func.General_SNR(WH,H, Es, N0);
LG = reshape(LG, K,M);
SNR = reshape(SNR, K,M);
[SNR1, LG1] = pre_func.F3_off_SNR(Vh, H_eq, wtx,wrx, Es, N0);
%%

 %% 2D-dft
% u_func =  Utility_functions;
% K = 16; M = 4;
% N = K*M;
% U_KM1 = u_func.U(K,M);
% %U_MK2 = u_func.U(M,K)';
% U_MK2 = eye(N);
% U_MK3 = u_func.U(M,K);
% %U_MK3 = eye(N);
% PI_KM = u_func.PI(K,M);
% PI_MK = u_func.PI(M,K);
% Ltx = exp(1j*2*pi*rand(N,1));
% wtx = 1/sqrt(M)*reshape(Ltx, K,M);
% Atx = 1/K*fft(wtx);
% wrx = 1./wtx;
% Af = PI_KM*U_MK3*diag(Ltx)*U_MK2*PI_MK*U_KM1*PI_KM;
% hf = randn(N,1)+1j*randn(N,1);
% N0 = 0.01;
% Es = 1;
% H = repmat(hf,1,N).*Af;
% pre_func = Precoding_Rx_functions();
% rxType = 'ZF';
% Vh = u_func.V(hf, K,M);
% H_eq = pre_func.CH_diag(Vh, rxType, Es, N0);
% WH = pre_func.General(H, rxType, Es, N0);
% [SNR, LG, sig, I] = pre_func.General_SNR(WH,H, Es, N0);
% LG = reshape(LG, K,M);
% SNR = reshape(SNR, K,M);
% [SNR1, LG1] = pre_func.F3_off_SNR(Vh, H_eq, wtx,wrx, Es, N0);


% func.General_SNR = @General_SNR;
%
% % Vh_eq = CH_diag(Vh,wtx, rxType, Es, N0, Mon)
% func.CH_diag = @CH_diag;
%
% %SNR, LG] = GFDM_SNR(Vh, H_eq, wtx,wrx, Es, N0)
% func.GFDM_SNR = @GFDM_SNR;
%
% %[SNR, LG] = OFDM_SNR(Vh, Es, N0)
% func.OFDM_SNR = @OFDM_SNR;
%
% % [SNR, LG] = SC_SNR(Vh, Es, N0)
% func.SC_SNR = @SC_SNR;
% %D = Demod(y,Wrx, K_set, dom)
%
% %[SNR, LG] = DFT_S_FDM_SNR(Vh, H_eq, wrx, Es, N0, isChZF)
% func.DFT_S_FDM_SNR = @DFT_S_FDM_SNR;
% %[SNR, LG] = F3_off_SNR(Vh, H_eq, wtx, wrx, Es, N0, isChZF)
% func.F3_off_SNR = @F3_off_SNR;
% %[SNR, LG] = F2_off_SNR(Vh, H_eq, wtx, wrx, Es, N0, isChZF)
% func.F2_off_SNR = @F2_off_SNR;