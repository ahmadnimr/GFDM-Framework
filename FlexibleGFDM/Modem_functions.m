function func = Modem_functions()
%% --------------- Memeber Variables --------------------------------------
%% --------------- Memeber functions Declaration --------------------------
% structure of GFDM matrix related  functions
%dom = 'TD'; 'FD';
% Symbol mapping
%D = Map(d, K_set, M_set, K,M)
func.Map = @Map;

%d = Demap(D, K_set, M_set)
func.Demap = @Demap;
% GFDM General Modem
%x = Mod(D,Wtx, K_set, dom)
func.Mod = @Mod;
%D = Demod(y,Wrx, K_set, dom)
func.Demod = @Demod;

% Costumized Modems
% OFDM
% DFT-S-OFDM
% OTFS
%

% EM, EK enable K/M-[i]fft
%Ds = spread(D, K_set, EM, EK)
func.spread = @spread;
%D = despread(Ds, K_set, EM, EK)
func.despread = @despread;

%Xw = window(W,X)
func.window = @window;

%Vx = transform(X,dom, ET)
func.transform = @transform;

%Y = detransform(Vy,dom, ET)
func.detransform = @detransform;

%x = allocation(Vx, Ea, index)
func.allocation = @allocation;

%Vy = deallocation(y, K,M, dom, Ea, index)
func.deallocation = @deallocation;

% --------- ref modem
%x = Basic_Mod(D,g,K,M, dom)
func.Basic_Mod = @Basic_Mod;
% D = Basic_Demod(x,gr,K,M)
func.Basic_Demod = @Basic_Demod;
% -------- convolution implementation
% x = Conv_TD_Mod(D,g,K,M)
func.Conv_TD_Mod = @Conv_TD_Mod;
%D = Conv_TD_Demod(y,gbr,K,M)
func.Conv_TD_Demod = @Conv_TD_Demod;
%xf = Conv_FD_Mod(D,gf,K,M, L_set)
func.Conv_FD_Mod = @Conv_FD_Mod;
%D = Conv_FD_Demod(yf,gfbar,K,M, L_set)
func.Conv_FD_Demod = @Conv_FD_Demod;
%% --------------- Including of library function --------------------------
%% --------------- Implementation -----------------------------------------

    function D = Map(d, K_set, M_set, K,M)
        D = zeros(K,M);
        K_on = length(K_set);
        M_on = length(M_set);
        D(K_set, M_set) = reshape(d, K_on, M_on);
    end

    function d = Demap(D, K_set, M_set)
        K_on = length(K_set);
        M_on = length(M_set);
        d = reshape(D(K_set, M_set), M_on*K_on,1);
    end

    function x = Mod(D,Wtx, K_set, dom)
        % General GFDM modulation
        X = spread(D, K_set);
        X = window(Wtx,X);
        X = transform(X,dom);
        x = allocation(X);
    end
    function D = Demod(y,Wrx, K_set, dom)
        %General GFDM demodulation
        [K,M] = size(Wrx);
        D = deallocation(y,K,M,dom);
        D = detransform(D,dom);
        D = window(Wrx, D);
        D = despread(D, K_set);
    end

    function Ds = spread(D, K_set, EM, EK)
        [K,M] = size(D);
        if nargin < 2
            K_set = 1:K;
            EM = true;
            EK = true;
        elseif nargin < 3
            EM = true;
            EK = true;
        elseif nargin < 4
            EK = true;
        end
        
        % EK: enable K-ifft
        % EM: Enable M-fft
        Ds = zeros(K,M);
        % DS = 1/K F_K^H D F_M
        if EM &&  M> 1
            Ds(K_set,:) = transpose(fft(transpose(D(K_set,:))));
        else
            Ds(K_set,:) = D(K_set,:);
        end
        if EK && K>1
            Ds = ifft(Ds);
        end
        
    end

    function D = despread(Ds, K_set, EM, EK)
        % EK: enable K-fft
        % EM: Enable M-ifft
        
        [K,M] = size(Ds);
        if nargin < 2
            K_set = 1:K;
            EM = true;
            EK = true;
        elseif nargin < 3
            EM = true;
            EK = true;
        elseif nargin < 4
            EK = true;
        end
        if EK && K>1
            D = fft(Ds);
        else
            D = Ds;
        end
        
        if EM && M>1
            D(K_set,:) = transpose(ifft(transpose(D(K_set,:))));
        end
    end

    function Xw = window(W,X, Ew)
         if nargin < 3
            Ew = true;
        end
        % W: KxM
        % X: KxM
        if Ew
        Xw = W.*X;
        else
           Xw = X; 
        end
    end
    function Vx = transform(X,dom, ET)
        if nargin < 3
            ET = true;
        end
        % X: KxM
        % Vx = 1/M F_M^H X^T
        [K,M] = size(X);
        if strcmp(dom,'TD')
            if ET && M>1
                Vx = ifft(X.');
            else
                Vx = X.';
            end
        else % FD
            if ET && K>1
                % X: KxM
                % Vx = F_K X
                Vx = fft(X);
            else
                Vx = X;
            end
        end
    end
    function Y = detransform(Vy,dom, ET)
        if nargin < 3
            ET = true;
        end
        if strcmp(dom,'TD')
            % Vy: MxK
            M = size(Vy,1);
            % Y = (F_M Vy).'
            if ET && M>1
                Y = transpose(fft(Vy));
            else
                Y = transpose(Vy);
            end
        else % FD
            % Vy: KxM
            K = size(Vy,1);
            % Y = 1/K*F_K^H Vy
            if ET && K>1
                Y = ifft(Vy);
            else
                Y = Vy;
            end
        end
    end

    function x = allocation(Vx, Ea, index)
        % this allocation is still part of the modem not GFDMA
        % the final vector is x = X(:).
        if nargin < 2
            % defualt is transpose
            X = Vx.';
        elseif nargin < 3
            if Ea
                X = Vx.';
            else % no transpose
                X = Vx;
            end
        else % indexing
            X(index) = Vx;
        end
        x = X(:);
    end

    function Vy = deallocation(y, K,M, dom, Ea, index)
        
        if nargin < 5
            Ea = true;
        end
        if nargin < 6
            if Ea
                if strcmp(dom,'TD')
                    Vy = reshape(y, K,M);
                else
                    Vy = reshape(y, M, K);
                end
                Vy = Vy.';
            else
                if strcmp(dom,'TD')
                    Vy = reshape(y, M,K);
                else
                    Vy = reshape(y, K, M);
                end
            end
        else % index
            y = y(index);
            
            if strcmp(dom,'TD')
                Vy = reshape(y, M,K);
            else
                Vy = reshape(y, K, M);
            end
            
        end
    end
    function x = Basic_Mod(D,g,K,M, dom)
        % refrence equatin
        N = M*K;
        tmp = zeros(N,1);
        x  = zeros(N,1);
        TD = strcmp(dom, 'TD');
        for m=0:M-1
            for k=0:K-1
                % Note that N_cp can be included here
                for n=0:N-1
                    if TD
                        tmp(n+1) = D(k+1,m+1) * g(mod(n-m*K, N)+1)*exp(2*1j*pi*n*k/K);
                    else
                        tmp(n+1) = D(k+1,m+1) * g(mod(n-k*M, N)+1)*exp(-2*1j*pi*n*m/M);
                    end
                end
                x = x+ tmp;
            end
        end
    end
    function D = Basic_Demod(x,gr,K,M,dom)
        % refrence equatin
        N = M*K;
        D  = zeros(K,M);
        TD = strcmp(dom, 'TD');
        for m=0:M-1
            for k=0:K-1
                % Note that N_cp can be included here
                for n=0:N-1
                    if TD
                        D(k+1,m+1) = D(k+1,m+1) + x(n+1)* conj(gr(mod(n-m*K,N)+1))*exp(-2*1j*pi*n*k/K);
                    else
                        D(k+1,m+1) = D(k+1,m+1) + x(n+1)* conj(gr(mod(n-k*M,N)+1))*exp(2*1j*pi*n*m/M);
                    end
                end
            end
        end
    end

    function x = Conv_TD_Mod(D,g,K,M)
        % 1/KF_K^H*D
        if K > 1
            Df = ifft(D);
        else
            Df = D;
        end
        G = cell(M,1);
        Db = cell(M,1);
        % V_MK_g^T
        Vgt = K*reshape(g,K,M);
        p = 0:M-1;
        % build matrix G^{m}
        for m=0:M-1
            %G{m}(q,p) = vg (<p-m>_M,q);
            G{m+1} = Vgt(:, mod(p-m,M)+1);
            % \bar{D}^{m}
            Db{m+1} = repmat(Df(:,m+1),1,M);
        end
        X = zeros(K,M);
        for m = 1:M
            X = X+ G{m}.*Db{m};
        end
        x = X(:);
    end

    function D = Conv_TD_Demod(y,gbr,K,M)
        % gbr = \gamma[-n]^* % demodulator pulse
        Gamm = cell(M,1);
        Yb = cell(M,1);
        % reshape received signal
        Vy = reshape(y,K,M);
        % reshape pulse
        Vgr = reshape(gbr,K,M);
        p = 0:M-1;
        % build matrix G^{m}
        for m=0:M-1
            %Gamm{m}(q,p) = vgr (<p-m>_M,q);
            Gamm{m+1} = Vgr(:, mod(p-m,M)+1);
            % \bar{D}^{m}
            Yb{m+1} = repmat(Vy(:,m+1),1,M);
        end
        X = zeros(K,M);
        for m = 1:M
            X = X+ Gamm{m}.*Yb{m};
        end
        if K > 1
            D = fft(X);
        else
            D = X;
        end
    end

    function xf = Conv_FD_Mod(D,gf,K,M, L_set)
        % L_set the set of non zero subcarriers
        L = numel(L_set);
        % F_M*D'
        if M > 1
            Df = fft(D.');
            
        else
            Df = D.';
        end
        %  V_KM_gf^T
        Vgft = reshape(gf,M,K);
        G = cell(L,1);
        Db = cell(L,1);
        q = 0:K-1;
        % build matrix G^{m}
        for li= 1:L
            %G{l}
            G{li} = repmat(Vgft(:,L_set(li)+1), 1, K);
            % Db{l} = FM D^T Pl^T: circuilar shift
            Db{li} = Df(:, mod(q - L_set(li), K)+1);
        end
        Xf = zeros(M,K);
        for li = 1:L
            Xf = Xf+ G{li}.*Db{li};
        end
        xf = Xf(:);
    end

    function D = Conv_FD_Demod(yf,gfbar,K,M, L_set)
        % L_set the set of non zero subcarriers
        L = numel(L_set);
        % constuct V_KM^{yf}^T
        Vy = reshape(yf,M,K);
        %reshape the filter
        Vgfr = reshape(gfbar,M,K)/K;
        % sub pulse shapes
        Gamm = cell(L,1);
        % sub signal
        Yb = cell(L,1);
        q = 0:K-1;
        % build matrix G^{m}
        for li= 1:L
            %G{l}
            Gamm{li} = repmat(Vgfr(:,L_set(li)+1), 1, K);
            % Db{l} = FM D^T Pl^T: circuilar shift
            Yb{li} = Vy(:, mod(q - L_set(li), K)+1);
        end
        Xf = zeros(M,K);
        for li = 1:L
            Xf = Xf+ Gamm{li}.*Yb{li};
        end
        if M > 1
            D = ifft(Xf);
        else
            D = Xf;
        end
        D = D.';
    end
%% --------------- END of Implementation ----------------------------------
end
