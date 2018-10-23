K = 8; M = 13;
N = K*M;
g = randn(K*M,1)+1j*randn(K*M,1);
g = g/norm(g);
% Generate random data symbols
D = randn(K,M)+1j*randn(K,M);
w_func = Window_functions();
hf = randn(N,1)+1j*randn(N,1);
Es = 1;
N0 = .01;
res_rand = Test_Equations_Receiver(g,hf,K,M, Es, N0);
g = w_func.Wtx2g(exp(1j*2*pi*rand(K,M)), 'TD');
g = g/norm(g);
res_orth = Test_Equations_Receiver(g,hf,K,M, Es, N0);
gf = zeros(N,1);
gf(1:M) = 1;
g = ifft(gf);
g = g/norm(g);
res_orth2 = Test_Equations_Receiver(g,hf,K,M, Es, N0);

function res2 = Test_Equations_Receiver(g,hf, K,M,Es, N0)
w_func = Window_functions();
rx_func = GFDM_Rx_functions();
gm = rand(K,1);
F = dftmtx(K);
gfm = fft(gm);
%-------------------------------------------------------
N = K*M;
gf = fft(g);
Af = w_func.Amtx(gf,K,M, 'FD');
Wtx_FD =  w_func.g2Wtx(gf,K,M,'FD');
Wrx_FD = w_func.Wtx2Wrx(Wtx_FD, 'ZF', 1, 1);

%MMSE test
% MMSE both
H = diag(hf)*Af;
tic
WH = rx_func.General(H, 'MMSE', Es, N*N0);
[SNR_ref_MMSE,LG_ref_MMSE, sig_MMSE, I_MMSE] = rx_func.General_SNR(WH,H, Es, N*N0);
toc
SNR_ref_MMSE = reshape(SNR_ref_MMSE, K,M);
LG_ref_MMSE = reshape(LG_ref_MMSE, K,M);

Vh = V(hf,K,M);
% MMSE channel + ZF 
tic
H_eq =  rx_func.CH_MMSE_exact(Vh,Wtx_FD, Es, N*N0, M);
[SNR_f, LG, A, sig, ICI] = rx_func.GFDM_SNR(Vh, H_eq, Wtx_FD,Wrx_FD, Es, N*N0);
SNR_f = repmat(SNR_f, 1,M);
toc 
% ZF test
WH = rx_func.General(H, 'ZF', Es, N*N0);
SNR_ref_ZF = rx_func.General_SNR(WH,H, Es, N*N0);
SNR_ref_ZF = reshape(SNR_ref_ZF, K,M);
% ZF channel
tic
Vh_eq  = rx_func.CH_diag(Vh,Wtx_FD, 'ZF', Es, N*N0, M);
SNR_d  = rx_func.GFDM_SNR(Vh, Vh_eq, Wtx_FD,Wrx_FD, Es, N*N0);
SNR_d = repmat(SNR_d, 1,M);
toc
% MMSE approx

Vh_eq = rx_func.CH_diag(Vh,Wtx_FD, 'MMSE', Es, N*N0, M);
[SNR_d_mmse, LG_d]= rx_func.GFDM_SNR(Vh, Vh_eq, Wtx_FD,Wrx_FD, Es, N*N0);
SNR_d_mmse = repmat(SNR_d_mmse, 1,M);

res2.SNR_MMSE_full_vs_Ref = norm(SNR_ref_MMSE(:)- SNR_f(:));
res2.SNR_ZF_diag_vs_Ref = norm(SNR_ref_ZF(:)- SNR_d(:));
res2.SNR_MMSE_full_vs_approx = norm(SNR_ref_MMSE(:)- SNR_d_mmse(:));
end
function V_PQ = V(x,P,Q)
% V = unvect_Q_P(x)^T
V_PQ = transpose(reshape(x, Q,P));
end