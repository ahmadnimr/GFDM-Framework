% create instance of the window functions
w_funcs = Window_functions();
% choose parameters
K = 12;
M = 6;
N = K*M;
g = rand(N,1);
gf = fft(g);
% set domain
dom = 'FD';
% set receiver type
rxType = 'ZF';
% Modulation matrix
A = w_funcs.Amtx(g,K,M, dom);
% Transmitter windwo
Wtx = w_funcs.g2Wtx(g,K,M,dom);
% receiver window
Wrx = w_funcs.Wtx2Wrx(Wtx, rxType, 1, 1);
% Receiver pulse
gamm = w_funcs.Wrx2gamma(Wrx, dom);
% Deodulation matrix
B = w_funcs.Amtx(gamm,K,M, dom);
% Test ZF relation
norm(B'*A-eye(N),'fro')