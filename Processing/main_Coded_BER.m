seed = 19092014; % seed for evaluation
K = 32; M = 16; Mc = 16;
N = K*M;
RxType = 'MMSE';
% coding
CodeType = 'TC-1-2';
payloadSize = 2*N/8-2;
% channel
ChType = 'EVA';
fs = 15360000 ;
Ncp = 64;
FD_options = false;
nu_M = 0:0.005:0.04; %fD =  nu_M*fs/Nt;
SNR = [12 18 24];

Pre_types = {'OFDM','GFDM-Chirp', 'SFFT', 'DFT-S-OFDM','2D-DFT', 'SC'};
Pre_inds = [1,2,4,5,3,6]; % PHY select

Nruns_SNR{1} = [100, 200, 400, 800, 1000, 1500,2000, 3000, 4000]';
Nruns_SNR{2} = [100, 200, 400, 600, 8000, 1000,1500, 2000, 3000]';
Nruns_SNR{3} = [100, 200, 200, 400, 600, 800,1000, 1500, 2000]';

Nruns_FD{3} = [8000, 6000, 5000, 4000, 2000, 1200, 800, 600, 400];
Nruns_FD{2} = [4000, 3000, 2000, 1500, 1000, 800, 600, 400, 300];
Nruns_FD{1} = [3000, 2000, 1500, 1000, 800, 600, 400, 300, 300];


%% Initialization
N_typ = numel(Pre_types);
N_SNR = numel(SNR);
N_fD = numel(nu_M);

% PHY params
N_types = numel(Pre_types);
params = cell(N_types,1);
for pre_i = Pre_inds
    params{pre_i} = GenParameters(K,M,Pre_types{pre_i}, RxType, Mc);
    params{pre_i}.seed = seed;
    Nruns{pre_i} = zeros(N_SNR, N_fD);
end
if FD_options
Nruns{1}(:,1) = Nruns_SNR{1};
Nruns{1}(:,2) = Nruns_SNR{2};
Nruns{1}(:,3) = Nruns_SNR{3};
Nruns{2} = 4*Nruns{1};
Nruns{2}(8,1) = 15000;
Nruns{2}(9,1) = 20000;
Nruns{3} = Nruns{2};
Nruns{4} = 3*Nruns{1};
Nruns{5} = Nruns{2};
Nruns{6} = Nruns{2};
else    
Nruns{1}(1,:) = Nruns_FD{2};
Nruns{1}(2,:) = Nruns_FD{3};
Nruns{1}(3,:) = Nruns_FD{2};

Nruns{2} = 4*Nruns{1};
%Nruns{2}(1,3) = 20000;
%Nruns{2}(1,2) = 15000;
Nruns{3} = Nruns{2};
Nruns{4} = 3*Nruns{1};
Nruns{5} = Nruns{2};
Nruns{6} = Nruns{2};
end

SNR_an = cell(N_typ,N_SNR,N_fD);
SER_an = cell(N_typ,N_SNR,N_fD);
SNR_sim = cell(N_typ,N_SNR,N_fD);
SER_sim = cell(N_typ,N_SNR,N_fD);

BER_sim = cell(N_typ,N_SNR,N_fD);
CBER_sim = cell(N_typ,N_SNR,N_fD);
FER_sim = cell(N_typ,N_SNR,N_fD);
UFER_sim = cell(N_typ,N_SNR,N_fD);
for fd_i = 1:N_fD
    for pre_i = 1:N_typ
        for nsnr = 1: N_SNR
            SNR_an{pre_i,nsnr,fd_i} = zeros( params{pre_i}.K, 1);
            SER_an{pre_i,nsnr,fd_i} = zeros(params{pre_i}.K, 1);
            SNR_sim{pre_i,nsnr,fd_i} = zeros(params{pre_i}.K, 1);
            SER_sim{pre_i,nsnr,fd_i} = zeros(params{pre_i}.K, 1);
            
            BER_sim{pre_i,nsnr,fd_i} = zeros(params{pre_i}.K, 1);
            CBER_sim{pre_i,nsnr,fd_i} = zeros(params{pre_i}.K, 1);            
            FER_sim{pre_i,nsnr,fd_i} = zeros(params{pre_i}.K, 1);
            UFER_sim{pre_i,nsnr,fd_i} = zeros(params{pre_i}.K, 1);
        end
    end
end

d_func = Data_functions();
ChCode = d_func.GenerateCode( CodeType, payloadSize,N,Mc);

ch_func = Channel_functions();
Nt = 2*N; % OFDM length [half allocation ]
fD =  nu_M*fs/Nt;
N0 = 10.^(-SNR/10);

tic
%% channel loop
ttemp = num2str(tic);
for fd_i = 1:N_fD
    for nsnr = 1: N_SNR
        for  pre_i = Pre_inds
            disp(['start (fd, snr, type): (' num2str(fd_i)...
                ','  num2str(nsnr) ',' num2str(pre_i) ')' '--' num2str(toc)]);
            rng(seed); % for the random data to be the same for all
            rchan = ch_func.GenFadingChannel( ChType, fD(fd_i), fs);
            %---------------------------
            [an, meas] = Process_FD_coded(params{pre_i},rchan, N0(nsnr),Nruns{pre_i}(nsnr,fd_i),Ncp,ChCode);
            %---------------------------
            % SER, SNR measure
            SNR_an{pre_i,nsnr,fd_i} =    mean(an.SNR(:));
            SER_an{pre_i,nsnr,fd_i} =   mean(an.SER(:));
            SNR_sim{pre_i,nsnr,fd_i} = 1./mean(meas.SymbolError(:));
            SER_sim{pre_i,nsnr,fd_i} =  mean(meas.SymSER(:));
            % BER, FER measure
            BER_sim{pre_i,nsnr,fd_i} = meas.BER;
            CBER_sim{pre_i,nsnr,fd_i} =  meas.CBER;
            FER_sim{pre_i,nsnr,fd_i} = meas.FER;
            UFER_sim{pre_i,nsnr,fd_i} = meas.UFER;
            
            save([ttemp 'temp']); % temp saving
        end
    end
end
if FD_options
fname = [CodeType ttemp 'FDotions'];
else
    fname = [CodeType ttemp 'SNRotions'];
end
save([fname ttemp]); % temp saving
toc