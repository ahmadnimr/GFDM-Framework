function func = Utility_functions()
%% --------------- Memeber Variables --------------------------------------
%% --------------- Memeber functions Declaration --------------------------
%V_PQ = V(x,P,Q)
func.V = @V;
%x = iV(Vx)
func.iV = @iV;
%Z_QP =Z(x, Q,P)
func.Z = @Z;
% x = iZ(Z_QP)
func.iZ = @iZ;
% Z_PQ = Zbar(xf, P,Q)
func.Zbar = @Zbar;
%xf = iZbar(Z_PQ)
func.iZbar = @iZbar;
%lam_1 = Ltx1(W)
func.Ltx1 = @Ltx1;
%lam_2 = Ltx2(W)
func.Ltx2 = @Ltx2;
%U_PQ = U(P,Q)
func.U = @U;
%PI_QP = PI(Q,P)
func.PI = @PI;
%G = CircMat(g)
func.CircMat = @CircMat;
% Xa = abs2(X)
func.abs2 = @abs2;
%% --------------- Including of library function --------------------------
%% --------------- Implementation -----------------------------------------
    function V_PQ = V(x,P,Q)
        % V = unvect_Q_P(x)^T
        V_PQ = transpose(reshape(x, Q,P));
    end
    function x = iV(Vx)
        Vx = Vx.';
        x = Vx(:);
    end
% ZAK TD
    function Z_QP = Z(x, Q,P)
        % Zak transform TD
        Z_QP = V(x, Q,P);
        if Q>1
            Z_QP = fft(Z_QP);
        end
    end
% inverse ZAK TD
    function x = iZ(Z_QP)
        % inverse Zak transform TD
        Q = size(Z_QP,1);
        Vx = Z_QP;
        if Q>1
            Vx = ifft(Z_QP);
            
        end
        x = iV(Vx);
    end
% ZAK FD
    function Zbar_PQ = Zbar(xf, P,Q)
        % Zak transform FD
        Zbar_PQ = V(xf,P,Q);
        if P >1
            Zbar_PQ = ifft(Zbar_PQ);
        end
    end
% inverse ZAK FD
    function xf = iZbar(Z_PQ)
        % inverse Zak transform FD
        P = size(Z_PQ,1);
        Vx = Z_PQ;
        if P > 1
            Vx = fft(Z_PQ);
        end
        xf = iV(Vx);
    end
% Matrix used in the decomposition
% ---- diagonal matrix types
    function lam_2 = Ltx2(W)
        % lam_2^_tx = diag(Vect(W_tx^T))
        W = W.';
        lam_2 = diag(W(:));
    end
    function lam_1 = Ltx1(W)
        % lam_1^_tx = diag(Vect(W_tx))
        lam_1 = diag(W(:));
    end
% unitary
    function U_PQ = U(P,Q)
        U_PQ = 1/sqrt(Q)*kron(eye(P),dftmtx(Q));
    end
% Commutation matrix
    function PI_QP = PI(Q,P)
        QP = Q*P;
        PI_QP = zeros(QP,QP);
        for q=0:Q-1
            for p=0:P-1
                PI_QP(p+q*P+1,q+p*Q+1) = 1;
            end
        end
    end
% generate circuila matrix
    function Am = CircMat(gm)
        K = numel(gm);
        gm = reshape(gm,K,1);
        Am = zeros(K,K);
        q = 0:K-1;
        for k=0:K-1
            Am(:,k+1) = gm(mod(q-k, K)+1);
        end
    end

    function xa = abs2(x)
        % effecient coputation of abs(x).^2
        xa = real(x).^2+imag(x).^2;
    end
%% --------------- END of Implementation ----------------------------------
end