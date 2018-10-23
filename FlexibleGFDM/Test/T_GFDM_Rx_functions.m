% Creat functions
clear all
rng(16445);
u_func = Utility_functions();
w_func = Window_functions();
rx_func = GFDM_Rx_functions();

% parameters
K = 16; M = 32;
N = K*M;
%g = randn(K*M,1)+1j*randn(K*M,1);
gf = zeros(N,1);
gf(1:2*M-4) = 1/2*(randn(2*M-4,1)+1j*randn(2*M-4,1));
gf(M/2+(1:M-2)) = 1;
g = ifft(gf);
g = g/norm(g);
% orthogonal 
%     g = w_func.Wtx2g(exp(1j*2*pi*rand(K,M)), 'TD');
%     g = g/norm(g);
W = w_func.g2Wtx(g, K,M, 'FD');
gf = fft(g);
% FD window 
Af = w_func.Amtx(gf,K,M, 'FD');
Wtx_FD =  w_func.g2Wtx(gf,K,M,'FD');
Wrx_FD = w_func.Wtx2Wrx(Wtx_FD, 'ZF', 1, 1);
L = 45;
h = 1; 1/sqrt(L*2)*(randn(L,1)+1j*randn(L,1));
hf = fft(h,N);


% channel
Es = 1;
N0 = 1;
H = diag(hf)*Af;
Vh = u_func.V(hf,K,M);

% General MMSE 
WH = rx_func.General(H, 'MMSE', Es, N*N0);
[SNR_ref_MMSE,LG_ref_MMSE, sig_MMSE, I_MMSE] = rx_func.General_SNR(WH,H, Es, N*N0);

SNR_ref_MMSE = reshape(SNR_ref_MMSE, K,M);
LG_ref_MMSE = reshape(LG_ref_MMSE, K,M);

WH = rx_func.General(H, 'ZF', Es, N*N0);
[SNR_zf,LG_ref_MMSE, sig_MMSE, I_MMSE] = rx_func.General_SNR(WH,H, Es, N*N0);
SNR_zf = reshape(SNR_zf, K,M);

% Full MMSE channel + ZF demodulator 
H_eq =  rx_func.CH_MMSE_exact(Vh,Wtx_FD, Es, N*N0, M);
[SNR_f, LG_f, A_f, sig_f, ICI_f, ISI_f] = rx_func.GFDM_SNR(Vh, H_eq, Wtx_FD,Wrx_FD, Es, N*N0);
SNR_f = repmat(SNR_f, 1,M); % Equal SNR per subsymbol
LG_f = repmat(LG_f, 1,M);

% diagonal MMSE channel + ZF demodulator 
Vh_eq  = rx_func.CH_diag(Vh,Wtx_FD, 'ZF', Es, N*N0, M);
[SNR_d, LG_d, A_d, sig_d, ICI_d, ISI_d] = rx_func.GFDM_SNR(Vh, Vh_eq, Wtx_FD,Wrx_FD, Es, N*N0);
SNR_d = repmat(SNR_d, 1,M);repmat(
LG_d = repmat(LG_d, 1,M);


% diagonal ZF channel + MMSE demodulator 
Vh_eq  = rx_func.CH_diag(Vh,Wtx_FD, 'ZF', Es, N*N0, M);
Vh_abs = u_func.abs2(Vh_eq);
A = sum(Vh_abs);

Omeg = repmat(A,K,1);

Wrx_FD = w_func.Wtx2Wrx(Wtx_FD, 'MMSE', Es, N0*Omeg);
[SNR_zm, LG_zm, A_zm, sig_zm, ICI_zm, ISI_zm] = rx_func.GFDM_SNR(Vh, Vh_eq, Wtx_FD,Wrx_FD, Es, N*N0);
SNR_zm = repmat(SNR_zm, 1,M);
LG_zm = repmat(LG_zm, 1,M);



res.SNR_MMSE_full_vs_Ref = norm(SNR_ref_MMSE(:)- SNR_f(:));
res.LG_full_vs_Ref = norm(LG_ref_MMSE(:)- LG_f(:));
res.SNR_MMSE_full_vs_approx = norm(SNR_ref_MMSE(:)- SNR_d(:));
res.SNR_MMSE_full_vs_approx_ZFMMSE = norm(SNR_ref_MMSE(:)- SNR_zm(:));
res.LG_full_vs_approx = norm(LG_ref_MMSE(:)- LG_d(:));
res