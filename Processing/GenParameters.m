function params = GenParameters(K,M,Type, RxType, Mc)
% Generate Precoding parameters
% Type = precoder type see options bellow
% RxType = {MMSE, ZF, MF}
% Mc: Modulation order% used to compute SER
params.K = K;
params.K_set = 1:K;
params.M = M;
params.N = K*M;
w_func = Window_functions();
params.Type = Type;
params.EW = true;
switch Type
    case 'OFDM'
        params.K = K*M;
        params.M = 1;
        params.K_set = 1:params.K;
        gf = zeros(K*M,1);
        gf(1) = 1;
        gf = gf/norm(gf);
        params.EM1 = false;
        params.EK2 = false;
        params.EK3 = false;
        params.EW = false;
    case 'GFDM-Chirp'
        gt = zeros(K*M,1);
        c1 = -cot(-pi/4)/2/K;
        gt(1:K) = exp(1j*2*pi*(0:K-1).^2*c1);
        gf = fft(gt);
        gf = gf/norm(gf);
        params.EM1 = true;
        params.EK2 = true;
        params.EK3 = true;
    case 'SFFT'
        params.EM1 = true;
        params.EK2 = true;
        params.EK3 = false;
        gf = zeros(K*M,1);
        gf(1:M) = ones(M,1);
        gf = gf/norm(gf);
    case 'DFT-S-OFDM'
        params.EM1 = true;
        params.EK2 = false;
        params.EK3 = false;
        gf = zeros(K*M,1);
        gf(1:M) = ones(M,1);
        gf = gf/norm(gf);
    case '2D-DFT'
        params.EM1 = true;
        params.EK2 = false;
        params.EK3 = true;
        gf = zeros(K*M,1);
        gf(1:M) = ones(M,1);
        gf = gf/norm(gf);
    case 'SC'
        params.EM1 = true;
        params.EK2 = false;
        params.EK3 = true;
        params.M = K*M;
        params.K = 1;
        params.K_set = 1;
        gf = zeros(params.N,1);
        gf(1:params.M) = ones(params.M,1);
        gf = gf/norm(gf);
        
    otherwise
        gf = zeros(K*M,1);
        gf(1:M) = ones(M,1);
        gf = gf/norm(gf);
end
params.gf = gf;
params.wtx = w_func.g2Wtx(params.gf,params.K,params.M,'FD');
params.wrx = w_func.Wtx2Wrx(params.wtx,'ZF');
params.Mc = Mc;
params.Atx = 1/K*fft(params.wtx);
params.rxType = RxType;
if strcmp(params.Type, 'SFFT')
    params.wtx = ones(params.K,params.M)*sqrt(params.K/params.M);
    params.wrx = 1./params.wtx;
elseif strcmp(params.Type, '2D-DFT')
    params.wtx = ones(params.K,params.M)/sqrt(params.N);
    params.wrx = 1./params.wtx;
end
end

