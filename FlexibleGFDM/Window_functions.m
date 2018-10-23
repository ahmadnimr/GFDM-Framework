function func = Window_functions ()
%% --------------- Memeber Variables --------------------------------------
%% --------------- Memeber functions Declaration --------------------------
% Core functions
%  dom = 'TD', 'FD'
% rxType = 'ZF', 'MMSE', 'MF'
% W: K x M tx/rx window

%A = Amtx(g,K,M, dom, K_set, M_set)
func.Amtx = @Amtx;

%Wtx = g2Wtx(g,K,M,dom):
func.g2Wtx = @g2Wtx;

% wrx = Wtx2Wrx(wtx, dom, psi_d, psi_v)
% Wtx of size K x M
% psi_d: size K x M diagonal matrix represnts the power per  symbol
% psi_v:  size K x M diagonal matrix represnts the noise power per channel
func.Wtx2Wrx = @Wtx2Wrx;

%g = Wtx2g(W, dom)
func.Wtx2g = @Wtx2g;

%gamma = Wrx2gamma(W, dom)% normalization is considered
func.Wrx2gamma = @Wrx2gamma;
%alpha is the roll-off
%[g_t, g_f] = RC_filter(K, M, alpha)
func.RC_filter = @RC_filter;
%% --------------- Including of library function --------------------------
util_func = Utility_functions();
Z = util_func.Z;
Zbar = util_func.Zbar;
iZ = util_func.iZ;
iZbar = util_func.iZbar;
abs2 = util_func.abs2;
%% --------------- Implementation -----------------------------------------
% Matrix Generation
    function A = Amtx(g,K,M, dom, K_set, M_set)
        % TD
        if nargin<5
            M_set = 1:M;
            K_set = 1:K;
        elseif nargin<6
            M_set = 1:M;
        end
        N = K*M;
        Kon = numel(K_set);
        Mon = numel(M_set);
        Non = Kon*Mon;
        A = zeros(N,Non);
        n=(0:N-1)';
        TD = strcmp(dom,'TD');
        n_on = 1;
        M_set = M_set-1;
        K_set = K_set-1;
        for m=1:Mon
            for k=1:Kon
                % Note that N_cp can be included here
                if TD
                    A(n+1,n_on) = g(mod(n-M_set(m)*K, N)+1).*exp(2*1j*pi*n*K_set(k)/K);
                else
                    A(n+1,n_on) = g(mod(n-K_set(k)*M, N)+1).*exp(-2*1j*pi*n*M_set(m)/M);
                end
                n_on = n_on+1;
            end
        end
    end
% Pulse to window
    function Wtx = g2Wtx(g,K,M,dom)
        if strcmp(dom,'TD')
            %Wtx_TD = K Z_MK ^T
            Wtx = K* Z(g, M,K).';
        else
            %Wtx_FD = K Z_MK
            Wtx = K* Zbar(g, K,M);
        end
    end
% Rx window
    function wrx = Wtx2Wrx(wtx, rxType, psi_d, psi_v)
        % Wtx: of size K x M
        % psi_d, psi_v are  of size K x M
        
        % find MMSE pulse
        if strcmp(rxType, 'ZF')
            % ZF
            wrx = 1./wtx;
        elseif strcmp(rxType, 'MF')
            %  MF
            wrx = conj(wtx);
        elseif strcmp(rxType, 'MMSE')
            % MMSE
            wrx = conj(wtx)./(abs2(wtx) + psi_v./psi_d);
        end
    end

% TX Pulse from TX window
    function g = Wtx2g(W, dom)
        K = size(W,2);
        if strcmp(dom,'TD')
            g = 1/K*iZ(W.');
        else
            g = 1/K*iZbar(W);
        end
    end
% RX Pulse from RX window% Normalization
    function g = Wrx2gamma(W, dom)
        N = numel(W);
        if strcmp(dom,'TD')
            g = iZ(W');
        else
            g = 1/N*iZbar(conj(W));
        end
    end
    function [g_t, g_f] = RC_filter(K, M, alpha)
        %K: number of samples per (sub)symbol, M: number of (sub)symbols
        % alpha: roll-off factor
        % Func: increasing function from -1 to 1;
        % shift sampling shift
        N = K*M;
        %alpha;
        if mod(M,2)== 0
            shift = 0.5;
        else
            shift = 0;
        end
        f=(-1/2:1/N: (1/2-1/N))'+shift/N;
        g_f = zeros(length(f),1);
        ind = abs(f)<= (1-alpha)/(2*K);
        g_f(ind) = 1;
        ind = ((abs(f)> (1-alpha)/(2*K))& (abs(f) <= (1+alpha)/(2*K)));
        f1 = f(ind);
        g_f(ind)= 1/2*(1+sin(-pi*K/alpha*(abs(f1)-1/(2*K))));
        g_f = sqrt(N)*g_f/norm(g_f);
        g_f = fftshift(g_f);
        g_t = ifft(g_f);
    end
%% --------------- END of Implementation ----------------------------------
end

