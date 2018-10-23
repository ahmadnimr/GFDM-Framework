%% Settings 
seed = 1010029; % random seed 
K = 32; M = 16;  Mc = 16;

ChType = 'AWGN';
fs = 8e6; fD =  0; % EVA
L = 24; expd = 0.5; % EXP, UNI, AWGN

SNR = 0:4:32;
N_SNR = numel(SNR);
SIM_on = true;
Nrun = 10*ones(N_SNR,1);
Nch = 2;
Pre_inds = [1 3];
RxType = 'MMSE'; %RxType = {'ZF', 'MMSE'};

%% Initialization
Pre_types = {'OFDM','GFDM-Chirp', 'SFFT', 'DFT-S-OFDM','2D-DFT', 'SC'};
% PHY params
N_types = numel(Pre_types);
params = cell(N_types,1);
for i = Pre_inds
    params{i} = GenParameters(K,M,Pre_types{i}, RxType, Mc);
    params{i}.seed = seed;
end
% (type index , snr index, Doppler index) for consistency  with plots
SNR_an = cell(N_types,N_SNR,1);
SER_an = cell(N_types,N_SNR,1);
SNR_sim = cell(N_types,N_SNR,1);
SER_sim = cell(N_types,N_SNR,1);
for i = Pre_inds
    for nsnr = 1: N_SNR
        SNR_an{i,nsnr,1} = zeros( params{i}.K, 1);
        SER_an{i,nsnr,1} = zeros(params{i}.K, 1);
        SNR_sim{i,nsnr,1} = zeros(params{i}.K, 1);
        SER_sim{i,nsnr,1} = zeros(params{i}.K, 1);
    end
end
N = K*M;
rng(seed);
Hf = GenDummyChannel(ChType, N, L, Nch, expd, fs, fD);
N0 = 10.^(-SNR/10);
%% channel loop
for nsnr = 1: N_SNR
    rng(seed*113); % for the random data
    disp(['start (snr): (' num2str(nsnr) ')' '--' num2str(toc) '\n']);
    for nch = 1:Nch
        hf = Hf(:,nch);
        for i = Pre_inds
            
            if SIM_on
                [eq, meas] = Process_SER(params{i},hf, N0(nsnr), Nrun(nsnr));
                SNR_an{i,nsnr,1} = SNR_an{i,nsnr,1} + mean(eq.SNR,2);
                SER_an{i,nsnr,1} = SER_an{i,nsnr,1} + mean(eq.SER,2);
                SNR_sim{i,nsnr,1} = SNR_sim{i,nsnr,1} + 1./mean(meas.SymbolError,2);
                SER_sim{i,nsnr,1} = SER_sim{i,nsnr,1} + mean(meas.SymSER,2);
            else
                eq = GenChEq(params{i},hf, 1, N0(nsnr));
                SNR_an{i,nsnr,1} = SNR_an{i,nsnr,1} + mean(eq.SNR,2);
                SER_an{i,nsnr,1} = SER_an{i,nsnr,1} + mean(eq.SER,2);
            end
        end
        
    end
    % average over channels
    for i =1:N_types
        SNR_an{i,nsnr,1} = SNR_an{i,nsnr,1}/Nch;
        SER_an{i,nsnr,1} = SER_an{i,nsnr,1}/Nch;
        SNR_sim{i,nsnr,1} = SNR_sim{i,nsnr,1}/Nch;
        SER_sim{i,nsnr,1} = SER_sim{i,nsnr,1}/Nch;
    end
    
end





