K = 4; M = 2;
N = K*M;
g = randn(K*M,1)+1j*randn(K*M,1);
g = g/norm(g);
% Generate random data symbols
D = randn(K,M)+1j*randn(K,M);
res = Test_Equations_Modem(g,D,K,M)

function res = Test_Equations_Modem(g,D,K,M)
m_func = Modem_functions();
w_func = Window_functions();
u_func = Utility_functions();
K_set = 1:K;
gf = fft(g);
%% Compute Decomposition parameters
% tx Window
Wtx_TD =  w_func.g2Wtx(g,K,M,'TD');
Wtx_FD =  w_func.g2Wtx(gf,K,M,'FD');
Wrx_TD = w_func.Wtx2Wrx(Wtx_TD, 'TD', 'ZF', 1, 1);
Wrx_FD = w_func.Wtx2Wrx(Wtx_FD, 'FD', 'ZF', 1, 1);
% diagonal TD
lam2_tx_TD = u_func.Ltx2(Wtx_TD);
lam1_tx_TD = u_func.Ltx1(Wtx_TD);
% diagonal FD
lam2_tx_FD = u_func.Ltx2(Wtx_FD);
lam1_tx_FD = u_func.Ltx1(Wtx_FD);
% unitary
UKM = u_func.U(K,M);
UMK = u_func.U(M,K);
% commutation
PKM = u_func.PI(K,M);
PMK = u_func.PI(M,K);
% U1, U2
Ut = UMK'*PMK*UKM*PKM;
Uf = UKM*PKM*UMK';

% V1, v2
Vt = PMK*UKM';
Vf = PKM*UMK;
%% Matrix formulas
A = w_func.Amtx(g,K,M, 'TD');
Af = w_func.Amtx(gf,K,M, 'FD');
% time domain decomposition
A2 = 1/sqrt(K) * Vt       * lam2_tx_TD * Uf;
A1 = 1/sqrt(K) * Vt * PKM * lam1_tx_TD * Ut;
% frequency domain decomposition
Af1 = sqrt(M) * Vf       * lam1_tx_FD * Ut;
Af2 = sqrt(M) * Vf * PMK * lam2_tx_FD * Uf;

%% Modulator test
x_ref =  m_func.Basic_Mod(D,g,K,M,'TD');
xf_ref = m_func.Basic_Mod(D,gf,K,M,'FD');

x_mat = A*D(:);
xf_mat = Af*D(:);

x4 = m_func.Mod(D,Wtx_TD,K_set,'TD');
xf4 = m_func.Mod(D,Wtx_FD,K_set,'FD');

%% Test implementation
x_TD = m_func.Conv_TD_Mod(D,g,K,M);
L_set = 0:K-1;
xf_FD = m_func.Conv_FD_Mod(D,gf,K,M, L_set);

%% Test Demodulator
De_TD = m_func.Demod(x_ref,Wrx_TD,K_set,'TD');
De_FD = m_func.Demod(xf_ref,Wrx_FD,K_set, 'FD');



gr = w_func.Wrx2gamma(Wrx_TD, 'TD');
gfr = w_func.Wrx2gamma(Wrx_FD, 'FD');

De_ref_TD = m_func.Basic_Demod(x_ref,gr,K,M,'TD');
De_ref_FD = m_func.Basic_Demod(xf_ref,gfr,K,M, 'FD');

%gbr = u_func.iZ(Wrx_TD.');
gbr = w_func.Wrx2gamma(conj(Wrx_TD), 'TD');
%gbfr = u_func.iZbar(Wrx_FD);
gbfr = w_func.Wrx2gamma(M*K*conj(Wrx_FD), 'FD');
D_TD = m_func.Conv_TD_Demod(x_ref,gbr,K,M);
L_set = 0:K-1;
D_FD = m_func.Conv_FD_Demod(xf_ref,gbfr,K,M, L_set);

%% check errors
res.TD_refvsFD_ref = norm(fft(x_ref)-xf_ref);

res.TD_refvsMat = norm(x_mat-x_ref);
res.FD_refvsMat = norm(xf_mat-xf_ref);

res.TD_refvs4step = norm(x4-x_ref);
res.FD_refvs4step = norm(xf4-xf_ref);
res.AvsAf = norm(Af-fft(A), 'fro');
res.TD_AvsA1 = norm(A-A1, 'fro');
res.TD_AvsA2 = norm(A-A2, 'fro');
res.TD_A2vsA1 = norm(A2-A1, 'fro');
res.FD_AvsAf1 = norm(Af-Af1, 'fro');
res.FD_AvsAf2 = norm(Af-Af2, 'fro');
res.FD_Af1vsAf2 = norm(Af1-Af2, 'fro');

res.TD_RefvsTD_imp = norm(x_ref-x_TD, 'fro');
res.FD_RefvsFD_imp = norm(xf_ref-xf_FD, 'fro');

res.TD_D_D4step = norm(De_TD-D, 'fro');
res.FD_D_D4step = norm(De_FD-D, 'fro');

res.TD_D_ref = norm(De_ref_TD-D, 'fro');
res.FD_D_ref = norm(De_ref_FD-D, 'fro');

res.TD_demod_imp = norm(D-D_TD, 'fro');
res.FD_demod_imp = norm(D-D_FD, 'fro');

end