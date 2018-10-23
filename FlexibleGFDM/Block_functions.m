function func = Block_functions()
%% --------------- Memeber Variables --------------------------------------
% It can be used to store variable, for example constant.
%% --------------- Memeber functions Declaration --------------------------
% reference to functions
% Xh = Filtering(X, h)
func.Filtering = @Filtering;
% X=Windowing(X, w, Ncs)
func.Windowing = @Filtering;
% Yo = RemoveCpCs( Yr, N, n0)
func.RemoveCpCs = @RemoveCpCs;
% X = AddCpCs( X0, N, Ncp, Ncs)
func.AddCpCs = @AddCpCs;
% Yr = BlockDemultiplexing(y, Ns, No)
func.BlockDemultiplexing = @BlockDemultiplexing;
% [xt, Ns] = BlockMultiplexing(X, No)
func.BlockMultiplexing = @BlockMultiplexing;
% S = Psd_th( G, R, Ns)
func.Psd_th = @Psd_th;
%[Sx, f] = PSD_Estimation( x, Nw)
func.PSD_Estimation = @PSD_Estimation;
%% --------------- Including of library function --------------------------
% call of other structures 
%% --------------- Implementation -----------------------------------------
% function implementation
% The block in X are transmitted with overlap and add of No samples
% No > 0, ovelapping
% No < 0, Zero padding
function [xt, Ns] = BlockMultiplexing(X, No)
[Nt, Nb] = size(X);
Ns = Nt-No;
N_sig = Ns*(Nb-1)+ Nt;
%x = X(:);
xt = zeros(N_sig,1);
xt(1:Nt) = X(:,1);
for nb=1:Nb-1
    xt((1:Nt)+nb*Ns) =  xt((1:Nt)+nb*Ns)+X(:,nb+1);
end
end
% extract blocks of length Nt 
% it is used to separate blocks before processing
% Ns block spacing
function Yr = BlockDemultiplexing(y, Ns, No)
Nt = No+Ns;
Ny = length(y);
Nb = ceil((Ny-No)/Ns);
Yr = zeros(Nt, Nb);
for nb=0:Nb-1
    Yr(:,nb+1) =  y((1:Nt)+nb*Ns);
end
end
% $$ X = X_0(<n-N_{cp}>_N,:),~ n = 0,\cdots, N+N_{cp}+N_{cs}-1  $$
function X = AddCpCs( X0, N, Ncp, Ncs)
No = Ncp+ Ncs;
Nt = N+No;
nt = 0:Nt-1;
X = X0(mod(nt-Ncp, N)+1,:);
end
% Extract blocks of length N, it start counting from n0, which can be the
% CP length itself.
function Yo = RemoveCpCs( Yr, N, n0)
Yo = Yr(n0+(1:N),:);
end
% Window over symbol 
% N_Cs length of Cs used for edges of the window.
function X = Windowing(X, w, Ncs)
Nb = size(X,2);
% Widow works only on the edges
for nb=1:Nb
    X(1:Ncs,nb) = w(1:Ncs).*X(1:Ncs,nb);
    X(end+1-(1:Ncs),nb) = w(end+1-(1:Ncs)).*X(end+1-(1:Ncs),nb);
end
end

% filter the wole block (filtering might include CP)
function Xh = Filtering(X, h)
[Nt, Nb] = size(X);
Nh = length(h)+Nt-1;
Xh = zeros(Nh,Nb);
for nb=1:Nb
    Xh(:,nb) = conv(h,X(:,nb));
end
end
% theoritical PSD of filter matrix G
function S = Psd_th( G, R, Ns)
% R is the resolution 
N = size(G,1);
if nargin<3
    Ns = N;
end
S = real(sum(abs(fft(G, N*R)).^2,2)/Ns);
end
% of signal x with window Nw
function [Sx, f] = PSD_Estimation( x, Nw)
%PSD_ESTIMATION Summary of this function goes here
%   Detailed explanation goes here
f = (0:Nw-1)'./Nw;
NT = floor(length(x)/Nw);
X = reshape(x(1:NT*Nw), Nw,NT);
Xf = abs(fft(X)).^2;
Sx = mean(Xf,2)/Nw;
end

%% --------------- END of Implementation ----------------------------------
end
