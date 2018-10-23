function H = GenDummyChannel(ChType, N, L, Nch, expd, fs, fd)
% Generate dummy channel
ch_func = Channel_functions();
if strcmp(ChType, 'AWGN')
    H = ones(N,Nch);
elseif strcmp(ChType, 'UNI')
    gain = ones(L,1)/sqrt(L);
    H = fft(ch_func.DummyChannel(gain, Nch), N);
elseif strcmp(ChType, 'EXP')
    gain = exp(-expd*(0:L-1)');
    gain = gain/norm(gain);
    H = fft(ch_func.DummyChannel(gain, Nch),N);
else %% use channel models
    H = zeros(N,Nch);
    x = zeros(2*N,1);
    x(1) = 1;
    rchan = ch_func.GenFadingChannel( ChType, fd, fs);
    for nch = 1:Nch
        [ h, ~ ] = ch_func.ApplyChannel( rchan, x, 0);
        hf = fft(h);
        H(1:N/2,nch) = hf(1+(1:N/2));
        H(end+1-(1:N/2),nch) = hf(end+1-(1:N/2));
    end
end
end