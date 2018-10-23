function func = Precoding_Rx_functions()
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

% Vh_eq = CH_diag(Vh,wtx, rxType, Es, N0, Mon)
func.CH_diag = @CH_diag;

%GFDM_SNR(Vh, H_eq, Atx, Es, N0, isChZF)
func.GFDM_SNR = @GFDM_SNR;

%[SNR, LG] = OFDM_SNR(Vh, Es, N0)
func.OFDM_SNR = @OFDM_SNR;

% [SNR, LG] = SC_SNR(Vh, Es, N0)
func.SC_SNR = @SC_SNR;
%D = Demod(y,Wrx, K_set, dom)

%[SNR, LG] = DFT_S_FDM_SNR(Vh, H_eq, wrx, Es, N0, isChZF)
func.DFT_S_FDM_SNR = @DFT_S_FDM_SNR;
%[SNR, LG] = F3_off_SNR(Vh, H_eq, wtx, wrx, Es, N0, isChZF)
func.F3_off_SNR = @F3_off_SNR;
%[SNR, LG] = F2_off_SNR(Vh, H_eq, wtx, wrx, Es, N0, isChZF)
func.F2_off_SNR = @F2_off_SNR;
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
% Diag approxination or exact for orthogonal window
    function Vh_eq = CH_diag(Vh, rxType, Es, N0)
        [K,M] = size(Vh);
        
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
                %Pm = sum(abs2(wtx(:,m)));
                Pm = N0/Es;
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
    function [SNR, LG] = GFDM_SNR(Vh, H_eq, Atx, Es, N0)
        % LG is the gain per subcarrier
        [K,M] = size(Vh);
        H = Vh.*H_eq;      
        A = zeros(K,1);
        zeta = zeros(K,1);        
           xi = zeros(K,1);
        for m=1:M
            % CM = 1/FK * [lam_rx *FK' * diag(hm)*diag(hm_eq) 1/K*FK] *lam_tx *FK'
            % Note that BM' = M*Am'
               Am = CircMat((Atx(:,m)));               
               Cm = Am'*(repmat(M*H(:,m),1,K).*Am);
               Em = (Am'.*repmat(M*H_eq(:,m).',K,1));
                A = A + diag(Cm);
                xi = xi + sum(abs2(Cm),2);
               zeta = zeta + sum(abs2(Em),2);
        end
        LG = 1/M*A;                          
        P =  Es*abs2(LG);
        I = Es*1/M*xi-P;
        sig = N0*zeta/M^2;
        SNR = P./(I+sig);                
    end
% SINR computation considering full equalizer
    function [SNR, LG] = OFDM_SNR(Vh, Es, N0)
        LG  = Vh;
        P = Es*abs2(LG);
        SNR = P./N0;
        LG  = Vh;
    end
    function [SNR, LG] = SC_SNR(Vh, H_eq, Es, N0, isChZF)
        % LG is the gain per subcarrier
        M = numel(Vh);
        if isChZF
            sig = H_eq'*H_eq;
            SNR = M^2./sig*N0/Es;
            LG  = 1;
        else
            Cm = H_eq.*Vh;
            sig = H_eq'*H_eq;
            Ag = sum(Cm);
            I = Cm'*Cm;
            P = abs2(Ag);
            SNR = P./(M*I-P + sig*N0/Es);
            LG  = Ag/M;
        end
    end

    function [SNR, LG] = DFT_S_FDM_SNR(Vh, H_eq, wtx,wrx, Es, N0)
        % LG is the gain per subcarrier
        [~,M] = size(Vh);
        Vh = Vh.*wtx;
        H_eq = H_eq.*wrx;
        H = Vh.*H_eq;
        A =  sum(H,2);
        xi =  sum(abs2(H),2);
        zeta = sum(abs2(H_eq),2);
        LG = 1/M*A;
        P =  Es*abs2(LG);
        I = Es*1/M*xi-P;
        sig = N0*zeta/M^2;
        SNR = P./(I+sig);
    end
 function [SNR, LG] = F3_off_SNR(Vh, H_eq, wtx, wrx, Es, N0)
       % effective channel  
       [~,M] = size(Vh);
       Vh = Vh.*wtx;
       H_eq = H_eq.*wrx;       
        h = H_eq(:).*Vh(:);
        A = mean(h);
        xi =  mean(abs2(h));
        zeta = mean(abs2(H_eq(:)))/M;
        LG = A;
        P =  Es*abs2(LG);
        I = Es*xi-P;
        sig = N0*zeta;
        SNR = P./(I+sig);
 end

function [SNR, LG] = F2_off_SNR(Vh, H_eq, wtx, wrx, Es, N0, isChZF)
        [~,M] = size(Vh);
       Vh = Vh.*wtx;
       H_eq = H_eq.*wrx;       
        h = H_eq(:).*Vh(:);
        A = mean(h);
        xi =  mean(abs2(h));
        zeta = mean(abs2(H_eq(:)))/M;
        LG = A;
        P =  Es*abs2(LG);
        I = Es*xi-P;
        sig = N0*zeta;
        SNR = P./(I+sig);
end
%% --------------- END of Implementation ----------------------------------
end
