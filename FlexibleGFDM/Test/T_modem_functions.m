% Create  functions 
m_func = Modem_functions();
w_func = Window_functions();
%% Modem parameters
K = 128; M = 16;
N = K*M;
g = randn(K*M,1)+1j*randn(K*M,1);
g = g/norm(g);
% Generate random data symbols
D = randn(K,M)+1j*randn(K,M);
K_set = 1:K;
gf = fft(g);
%% Create Tx widows 
Wtx_TD =  w_func.g2Wtx(g,K,M,'TD');
Wtx_FD =  w_func.g2Wtx(gf,K,M,'FD');
%% Create Rx widows 
Wrx_TD = w_func.Wtx2Wrx(Wtx_TD, 'ZF', 1, 1);
Wrx_FD = w_func.Wtx2Wrx(Wtx_FD, 'ZF', 1, 1);
%% Rx pulses
gr = w_func.Wrx2gamma(Wrx_TD, 'TD');
gfr = w_func.Wrx2gamma(Wrx_FD, 'FD');
%% Rx FD convoutional implementation pulses 
gbr = w_func.Wrx2gamma(conj(Wrx_TD), 'TD');
gbfr = w_func.Wrx2gamma(M*K*conj(Wrx_FD), 'FD');
%% Modulator
% Basic reference modem
x_ref =  m_func.Basic_Mod(D,g,K,M,'TD');
xf_ref = m_func.Basic_Mod(D,gf,K,M,'FD');
% 4 step modulator
x4 = m_func.Mod(D,Wtx_TD,K_set,'TD');
xf4 = m_func.Mod(D,Wtx_FD,K_set,'FD');
% TD convolution 
x_TD = m_func.Conv_TD_Mod(D,g,K,M);
% FD convolution
L_set = 0:K-1;
xf_FD = m_func.Conv_FD_Mod(D,gf,K,M, L_set);
%% Demodulator
% Basic reference demod
De_ref_TD = m_func.Basic_Demod(x_ref,gr,K,M,'TD');
De_ref_FD = m_func.Basic_Demod(xf_ref,gfr,K,M, 'FD');
% 4 step demodulator
De_TD = m_func.Demod(x_ref,Wrx_TD,K_set,'TD');
De_FD = m_func.Demod(xf_ref,Wrx_FD,K_set, 'FD');
% TD convolution 
D_TD = m_func.Conv_TD_Demod(x_ref,gbr,K,M);
% FD convolution 
D_FD = m_func.Conv_FD_Demod(xf_ref,gbfr,K,M, L_set);
%% check errors
res.TD_refvsFD_ref = norm(fft(x_ref)-xf_ref);

res.TD_refvs4step = norm(x4-x_ref);
res.FD_refvs4step = norm(xf4-xf_ref);
res.TD_RefvsTD_imp = norm(x_ref-x_TD, 'fro');
res.FD_RefvsFD_imp = norm(xf_ref-xf_FD, 'fro');

res.TD_D_ref = norm(De_ref_TD-D, 'fro');
res.FD_D_ref = norm(De_ref_FD-D, 'fro');


res.TD_D_D4step = norm(De_TD-D, 'fro');
res.FD_D_D4step = norm(De_FD-D, 'fro');

res.TD_demod_imp = norm(D-D_TD, 'fro');
res.FD_demod_imp = norm(D-D_FD, 'fro');

res