function func = GFDM_Rx_functions()
%% --------------- Memeber Variables --------------------------------------
%% --------------- Memeber functions Declaration --------------------------
% LG gain, for linear it is of size Non
% For GFDM it is of size K\times 1
% N0 noise of the linear model after N-DFT if applicale
% Mon number of active subsubsymbols with Kon = K
%SNR: total SNR
%LG: gain
% AFMd+ v, sig_v =  sig+ICI

%General Linear rceiver
%WH = General(H, rxType, Es, N0)
func.General = @General;
%[SNR, LG, sig, I]  = General_SNR(WH,H, Es, N0)
func.General_SNR = @General_SNR;

%GFDM Receiver-Channel equalization

%H_eq = CH_MMSE_exact(Vh,wtx, Es, N0, Mon)
func.CH_MMSE_exact = @CH_MMSE_exact;

% Vh_eq = CH_diag(Vh,wtx, rxType, Es, N0, Mon)
func.CH_diag = @CH_diag;

%SNR, LG, A, sig, ICI, ISI] = GFDM_SNR(Vh, H_eq, wtx,wrx, Es, N0)
func.GFDM_SNR = @GFDM_SNR;

%D = Demod(y,Wrx, K_set, dom)

%% --------------- Including of library function --------------------------
util_func = Utility_functions();
% utility
abs2 = util_func.abs2;
CircMat = util_func.CircMat;
%% --------------- Implementation -----------------------------------------
%-------------------------------
% SISO Receiver equations [FD] Full allocated subcarriers
%Ref equalizer
%H = repmat(h,1,Non).*Af_on;
    function WH = General(H, rxType, Es, N0)
        % N0 is the variance of the noise in the model after N-DFT
        [~,Non] = size(H);
        if strcmp(rxType, 'ZF')
            % ZF
            WH = (H'*H)\H';
        elseif strcmp(rxType, 'MF')
            %  MF
            WH = H';
        elseif strcmp(rxType, 'MMSE')
            WH = (H'*H+N0/Es*eye(Non))\H';
        end
    end

% Full Channel Equalization for non-orthogonal matrix, full allocation
    function H_eq = CH_MMSE_exact(Vh,wtx, Es, N0, Mon)
        [K,M] = size(wtx);
        if nargin < 5
            Mon = M;
        end
        H_eq = cell(M,1);
        
        for m=1:M
            % To compute the invers [AmAm']^{-1}
            Dm = N0/K/Mon/Es*fft(1./abs2(wtx(:,m)));
            %[ lam^{h}_m*lam^{h}_m' + 1/K*FK*[N*N0/Es/Mon *lam_m*lam^H]^{-1} *FK']^{-1};lma^{h}_m';
            H_eq{m} =  (diag(abs2(Vh(:,m)))+CircMat(Dm))\diag(conj(Vh(:,m)));
        end
    end

% Diag approxination or exact for orthogonal window
    function Vh_eq = CH_diag(Vh,wtx, rxType, Es, N0, Mon)
        [K,M] = size(wtx);
        
        if strcmp(rxType, 'ZF')
            % ZF
            Vh_eq = 1./Vh;
        elseif strcmp(rxType, 'MF')
            %  MF
            Vh_eq = conj(Vh);
        elseif strcmp(rxType, 'MMSE')
            % MMSE
            Vh_eq = zeros(K,M);
            for m=1:M
                % To compute the invers [AmAm']^{-1}
                % Pm = N0/K/Mon/Es*sum(1./abs2(wtx(:,m)));
                Pm = sum(abs2(wtx(:,m)));
                Pm = N0*K/Mon/Es/Pm;
                %[ lam^{h}_m*lam^{h}_m' + 1/K*FK*[N0/Es/Mon *lam_m*lam^H]^{-1} *FK']^{-1};lma^{h}_m';
                Vh_eq(:,m) =  conj(Vh(:,m))./(abs2(Vh(:,m)) + Pm);
            end
        end
        
    end
%
    function  [SNR, LG, sig, I] = General_SNR(WH,H, Es, N0)
        % The model
        % WH*y = WH*H d + WH*v;
        % LG is the gain matrix
        Non = size(H,2);
        C = WH*H;
        LG = diag(C);
        P = abs2(LG);
        I = abs2(C)*ones(Non,1)-P;
        sig = abs2(WH)*ones(Non,1); % sum;
        SNR = P./(I+sig*N0/Es);
        % output
        sig = sig*N0;
        I = Es*I;
    end

% SINR computation considering full equalizer
    function [SNR, LG, A, sig, ICI, ISI] = GFDM_SNR(Vh, H_eq, wtx,wrx, Es, N0, isChZF, isDemZF)
        % LG is the gain per subcarrier
        [K,M] = size(wtx);
        A = zeros(K,M);
        I =  zeros(K,1);
        ISI =  zeros(K,1);
        sig = zeros(K,1);
        EQ_FULL = iscell(H_eq);
        for m=1:M
            % CM = 1/FK * [lam_rx *FK' * diag(hm)*diag(hm_eq) 1/K*FK] *lam_tx *FK'
            if EQ_FULL
                Cm = H_eq{m}.*repmat(Vh(:,m).',K,1);
                Cm = transpose(fft(transpose(ifft(Cm))));
                Cm = Cm.*(wrx(:,m)*wtx(:,m).');
                Cm = transpose(ifft(transpose(fft(Cm))));
                % Bm'* lam_heq
                %Vm = CircMat(1/K*fft(wrx(:,m)))*H_eq{m};
                Vm = fft(repmat(wrx(:,m), 1,K).*ifft(H_eq{m}));
            else
                if isChZF
                    if isDemZF
                        Cm = eye(K);
                    else
                        Cm = CircMat(1/K*fft(wrx(:,m).*wtx(:,m)));
                    end
                else
                    Cm = CircMat(ifft(Vh(:,m).*H_eq(:,m)));
                    Cm = Cm.*(wrx(:,m)*wtx(:,m).');
                    Cm = transpose(ifft(transpose(fft(Cm))));
                end
                % Bm'* lam_heq
                Vm = CircMat(1/K*fft(wrx(:,m))).*repmat(H_eq(:,m).',K,1);
            end
            
            A(:,m) =  diag(Cm);
            I = I + abs2(Cm)*ones(K,1); % sum
            ISI =  ISI+abs2(A(:,m));
            sig = sig + abs2(Vm)*ones(K,1); % sum
        end
        Ag = A*ones(M,1); % sum
        P = abs2(Ag);
        SNR = P./(M*I-P + sig*N0/Es);
        LG  = Ag/M;
        sig = N0/M^2*sig;
        ICI = Es/M*(I -ISI);
        ISI =  Es/M*ISI-Es*P/M^2;
    end
%% --------------- END of Implementation ----------------------------------
end
