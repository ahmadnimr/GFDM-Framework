function [eq, meas] = Process_SER(params,hf, N0,Nrun)
% To Test with one fading channel realization and
% Nrun noise realizations
% N0 noise variance (scalar)
% hf: channel gain
f_data = Data_functions();
f_channel = Channel_functions();
N = params.N;

meas = CreateMeasurments();

eq = GenChEq(params,hf, 1, N0);

for nr =1:Nrun
    % gen data
    Tx.Di = reshape(f_data.QAM_gen(params.Mc, params.N), params.K, params.M);
    Tx.D = f_data.QAM_mod(Tx.Di, params.Mc);
    
    % modulator
    xf = Precode(Tx.D, params);
    % apply channel
    vf = f_channel.GenRandomNoise([N,1],N0);
    yf = hf.*xf + vf;
    % demodulator
    Rx.D = Equalization(yf, eq, params);
    Rx.Di = f_data.QAM_demod(Rx.D, params.Mc);
    
    Tx.InfoBits = 0;
    Tx.EncodedBits = 0;
    Rx.InfoBits = 0;
    Rx.EncodedBits = 0;
    meas =  meas.Measure(meas,Tx, Rx);
end

meas = meas.Avg(meas);

end



