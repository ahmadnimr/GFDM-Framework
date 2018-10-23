K = 16; M = 8;
N = K*M;
g = randn(K*M,1)+1j*randn(K*M,1);
g = g/norm(g);
% Generate random data symbols
D = randn(K,M)+1j*randn(K,M);
w_func = Window_functions();
rx_func = GFDM_Rx_functions();
util_func = Utility_functions();
m_func = Modem_functions();

Vhf = randn(K,M)+1j*randn(K,M);
hf = util_func.iV(Vhf);
Es = 1;
N0 = .01;
K_set = 1:K;

gf = fft(g);
Wtx_FD =  w_func.g2Wtx(gf,K,M,'FD');
Wrx_FD = w_func.Wtx2Wrx(Wtx_FD, 'FD', 'ZF', 1, 1);



H_eq =  rx_func.CH_MMSE_exact(Vhf,Wtx_FD, Es, N*N0, M);
% X = util_func.V(m_func.Mod(D,Wtx_FD, K_set, 'FD'), K,M);
% % Rx
% Y = X.*Vhf+Vv;

[SNR_f, LG, A, sig, ICI] = rx_func.GFDM_SNR(Vhf, H_eq, Wtx_FD,Wrx_FD, Es, N*N0);
Cm = cell(M,1);
Vm = cell(M,1);
for m=1:M
      Cm{m} = H_eq{m}.*repmat(Vhf(:,m).',K,1);
        Cm{m} = transpose(fft(transpose(ifft(Cm{m}))));
        Cm{m} = Cm{m}.*(Wrx_FD(:,m)*Wtx_FD(:,m).');
        Cm{m} = transpose(ifft(transpose(fft(Cm{m}))));
        % Bm'* lam_heq
        Vm{m} = util_func.CircMat(1/K*fft(Wrx_FD(:,m)))*H_eq{m};
end

%Y_eq = Y;
Atild = zeros(K,M);
Etild = zeros(K,M);
Vtild = zeros(K,M);
sig_e = zeros(K,M);
ICI_e = zeros(K,M);
NR = 10000;
Imk = zeros(NR,1);
for nr =1:NR
    D = sqrt(1/2)*(rand(K,M)+1j*rand(K,M));
    % General GFDM modulation
    X = m_func.spread(D, K_set);
    X = m_func.window(Wtx_FD,X);
    X = m_func.transform(X,'FD');
    Vv = sqrt(N*N0/2)*(randn(K,M)+1j*randn(K,M));
    Y = X.*Vhf+Vv;
    DfK = transpose(fft(D.'));
    % equalization
    for m=1:M
        Y_eq(:,m) =  H_eq{m}*Y(:,m);
        
        for k=1:K
            Atild(k,m) = Cm{m}(k,k);
            Etild(k,m) = Cm{m}(k,:)*DfK(:,m)-Cm{m}(k,k)*DfK(k,m);
            Vtild(k,m) = Vm{m}(k,:)*Vv(:,m);
        end
    end
    sig_e = sig_e+ util_func.abs2(ifft(Vtild.').');
    ICI_e = ICI_e+ util_func.abs2(ifft(Etild.').');
    Imk(nr) =  Etild(k,m) ;
    De = m_func.detransform(Y_eq,'FD');
    De = m_func.window(Wrx_FD, De);
    De = m_func.despread(De, K_set, 0,1);
    Vb = Vtild + Etild;
    De1 = D;
    DfM = fft(D.'); % M-FFT
    for k =1:K
        De1(k,:) = transpose(diag(Atild(k,:))*DfM(:,k)) +   Vb(k,:);
    end
end
% estimate model;

sig_e = sig_e/NR;
ICI_e = ICI_e/NR;


