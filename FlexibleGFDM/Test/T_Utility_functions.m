K = 8; 
M = 8;
N = K*M;

func = Utility_functions();
g = rand(N,1);
gf = fft(g);
Zg = func.Z(g, M,K);
Zgbar = func.Zbar(gf, K, M);

zg = reshape(Zg.', N,1);
zgbar = reshape(Zgbar, N,1);

UKM = func.U(K,M);
UMK = func.U(M,K);
PMK = func.PI(M,K);
PKM = func.PI(K,M);
FN= dftmtx(N);
S = sqrt(M)*PMK*UKM*PMK';
zg_e = S*g;
Sbar = 1/sqrt(K)*UMK'*PKM'*FN;
zgbar_e = Sbar*g;

err1 = norm(zg_e-zg)
err12 = norm(zgbar_e-zgbar)
B = FN;
[V D] = eig(B);
d = diag(D);
index = find(abs(d)<1e-6);
v = V(:,index(1));
%v =  V(:,index)*rand(length(index),1);
vf = fft(v);

Zg = func.Z(v, M,K);
Zgbar = func.Zbar(vf, K, M);

zg = reshape(Zg.', N,1);
zgbar = reshape(Zgbar, N,1);
zg_e = S*v;
zgbar_e = Sbar*v;

err1 = norm(zg_e-zg)
err12 = norm(zgbar_e-zgbar)
