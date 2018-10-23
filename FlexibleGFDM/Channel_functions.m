function func = Channel_functions()
%% --------------- Memeber Variables --------------------------------------
% It can be used to store variable, for example constant.
%% --------------- Memeber functions Declaration --------------------------
% reference to functions
%gain power profile, L number of taps, NR number of realizations
%H = DummyChannel(gain, L, NR)
func. DummyChannel = @DummyChannel;
%v = GenRandomNoise(siz, 0)
func.GenRandomNoise = @GenRandomNoise;
% fs sampling frequency 
% fD Doppler frequency
%rchan = GenFadingChannel( ChType, fD, fs)
func.GenFadingChannel = @GenFadingChannel;
% He estimatated channel after removing CP.
%[ He, Y ] = ApplyChannel( rchan, X, Ncp)
func.ApplyChannel = @ApplyChannel;
%% --------------- Including of library function --------------------------
% call of other structures 
%% --------------- Implementation -----------------------------------------
% function implementation
    function H = DummyChannel(gain, NR)
	 L = numel(gain);
        H = repmat(gain, 1, NR).*GenRandomNoise([L, NR], 1);
    end
function v = GenRandomNoise(siz, N0)
v = sqrt(N0/2) * (randn(siz)+1j*randn(siz));
end

function [ He, Y ] = ApplyChannel( rchan, X, Ncp)
release(rchan);
rchan.Seed = rchan.Seed+1; % change realization
[Ns, NB] = size(X);
D = zeros(Ns,NB);
% Estimate the channel appling a pulse
D(Ncp+1,:) = 1;
He = zeros(size(D));
for nb=1:NB
    He(:,nb) = step(rchan, D(:,nb));
end
% reset the channel to first state which correspond to the estimation
reset(rchan);
y = step(rchan, X(:));
Y = reshape(y,Ns, NB);
end

% ChType = EPA, EVA, TGn, FLAT
function rchan = GenFadingChannel( ChType, fD, fs)
%GENFADINGCHANNEL Summary of this function goes here
%   Detailed explanation goes here
EPAT = [0, 30, 70, 90, 110, 190, 410]*1E-9;
EVAT = [0, 30, 150, 310, 370, 710, 1090, 1730, 2510]*1E-9;

EPAP = [0, -1, -2, -3, -8, -17.2, -20.8];
EVAP = [0, -1.5, -1.4, -3.6, -0.6, -9.1, -7.0, -12.0 -16.9];

tau = [0 10 20 30 40 50 60 70 80]*1e-9;   % Path delays, in seconds

% Average path gains of cluster 1, in dB
pdb1 = [0 -5.4 -10.8 -16.2 -21.7 -inf -inf -inf -inf];
% Average path gains of cluster 2, in dB
pdb2 = [-inf -inf -3.2 -6.3 -9.4 -12.5 -15.6 -18.7 -21.8];
% Total average path gains for both clusters, in dB
pdb = 10*log10(10.^(pdb1/10) + 10.^(pdb2/10));

switch ChType
    case 'EPA'
        pathDelays = EPAT;
        avgPathGains = EPAP;
    case 'EVA'
        pathDelays = EVAT;
        avgPathGains = EVAP;
    case 'TGn'
         pathDelays = tau;
        avgPathGains = pdb;
    otherwise
        pathDelays = 0;
        avgPathGains = 1;
end
rchan = comm.RayleighChannel('SampleRate',fs, ...
    'PathDelays',pathDelays, ...
    'AveragePathGains',avgPathGains, ...
    'MaximumDopplerShift',fD,...
    'RandomStream','mt19937ar with seed', ...
    'Seed',22);
end
%% --------------- END of Implementation ----------------------------------
end
