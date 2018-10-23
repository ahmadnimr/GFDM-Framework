%% Measurmet structure
% all functions related to measurments
function measure = CreateMeasurments()
measure.NB = 0;
measure.SymSER = 0;
measure.SymbolError = 0;
measure.CodedERR = 0;
measure.CodedLength = 0;
measure.FrameERR = 0;
measure.UncodedFrameERR = 0;
measure.FrameCount = 0;
measure.UnCodedError = 0;
measure.UnCodedLength = 0;

%% functions
measure.Measure = @Measure;
measure.Avg = @Avg;
%% Librarary
util_func = Utility_functions();
abs2 = util_func.abs2;

    function measure = Measure(measure,Tx, Rx)
        measure.NB = measure.NB + 1;
        measure.SymSER = measure.SymSER + (Tx.Di~= Rx.Di);
        measure.SymbolError = measure.SymbolError+ abs2(Tx.D-Rx.D);
        codedError_frame  = sum(Tx.InfoBits~=Rx.InfoBits);      
        codedError  = sum(codedError_frame);
        measure.CodedERR = measure.CodedERR + codedError;
        measure.CodedLength = measure.CodedLength + numel(Rx.InfoBits);
        N_frame = size(Tx.InfoBits, 2);
        measure.FrameERR = measure.FrameERR + sum(codedError_frame > 0);
        
        
        measure.FrameCount = measure.FrameCount + N_frame;
        uncodedError = sum(sum(Tx.EncodedBits~=Rx.EncodedBits));
        UncodedError_frame  = sum(Tx.EncodedBits~=Rx.EncodedBits);
        measure.UnCodedError = measure.UnCodedError + uncodedError;
        measure.UnCodedLength = measure.UnCodedLength + numel(Rx.EncodedBits);
          measure.UncodedFrameERR = measure.UncodedFrameERR + sum(UncodedError_frame > 0);
    end
    function meas = Avg(meas)
        meas.SymSER = meas.SymSER/meas.NB;
        meas.SymbolError = meas.SymbolError/meas.NB;
        meas.CBER = meas.CodedERR/meas.CodedLength;
        meas.BER = meas.UnCodedError/meas.UnCodedLength;
        meas.FER = meas.FrameERR/meas.FrameCount;
        meas.UFER = meas.UncodedFrameERR/meas.FrameCount;
    end
end
