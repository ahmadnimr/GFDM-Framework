function func = Data_functions()
%% --------------- Memeber Variables --------------------------------------
% It can be used to store variable, for example constant.
%% --------------- Memeber functions Declaration --------------------------
% reference to functions
%[ser, PeI] = QAM_SER(Mc, SNR_l)
func.QAM_SER_an = @QAM_SER_an;
%d_i = QAM_gen(Mc, Non)
func.QAM_gen = @QAM_gen;
%d_qam = QAM_mod(Mc, d_i)
func.QAM_mod = @QAM_mod;
%d_i = QAM_demod(Mc, d_qam)
func.QAM_demod = @QAM_demod;
%ser = Ser(di, di_e)
func.Ser = @Ser;
%diff = Nmse(dqam, dqam_e)
func.Nmse = @Nmse;

%GenerateCode( CodeType,Payload_size,Non,Mc, seed )
func.GenerateCode = @GenerateCode;
%[ qam_data, info_bits, encoded_bits] = GenerateEncodedData( ChCode)
func.GenerateEncodedData = @GenerateEncodedData;
%[ r, y ] = CRC( x, g )
func.CRC = @CRC;
%[ bits_e, encoded_bits_e ] = Decoder( ChCode, data_e, Esno, h)
func.Decoder = @Decoder;
%[ qam_data, encoded_bits] = Encoder( ChCode, info_bits)
func.Encoder = @Encoder;


%% --------------- Including of library function --------------------------
% call of other structures
util_func = Utility_functions();
abs2 = util_func.abs2;
%% --------------- Implementation -----------------------------------------
% function implementation
    function d_i = QAM_gen(Mc, Non)
        d_i = randi([0, Mc-1], Non,1);
    end
    function d_qam = QAM_mod(d_i, Mc)
        %Generate qam data
        d_qam = qammod(d_i, Mc, 'UnitAveragePower',true);
    end
    function d_i = QAM_demod(d_qam, Mc)
        %Generate qam data
        d_i = qamdemod(d_qam, Mc, 'UnitAveragePower',true);
    end
    function [ser, PeI] = QAM_SER_an(M, SNR_l)
        % y = d+ v.
        % SNR_l = E{|d|^2}/E{|v|^2}
        % error probability per componenet
        PeI = 2*(sqrt(M)-1)/sqrt(M)* Q(sqrt(3/(M-1)*SNR_l));
        ser = 1-(1-PeI).^2;
    end
    function ser = Ser(di, di_e)
        % number of errors
        ser = (di ~= di_e);
        ser = mean(ser(:));
    end

    function diff = Nmse(dqam, dqam_e)
        % square of
        diff =  abs2(dqam-dqam_e);
        diff = mean(diff(:));
    end
    function out=Q(x)
        out=erfc(x/sqrt(2))/2;
    end
% Mc: QAM Modulation length 
% active data
% payload size in Bytes
% Code type = 'CC-1-2', 'CC-2-3', 'CC-3-4', 'CC-5-6' : Convolution
% Code type = 'TC-1-2', 'TC-2-3', 'TC-3-4', 'TC-5-6' : Turbo
% all parameters related are geneaten internally 
% random interleavers are used
% punctures are used
% sdeed is used to for random interleaver 
function ChCode = GenerateCode( CodeType,Payload_size,Non,Mc, seed )
if nargin < 5
    seed = 199;
end
% Nf : number of frames 
bps = log2(Mc);
if Mc == 4
    ChCode.Constellation  = CreateConstellation( 'QAM', Mc);
else
    ChCode.Constellation = CreateConstellation( 'QAM', Mc,'gray');
end
ChCode.Constellation = reshape(ChCode.Constellation,numel(ChCode.Constellation),1);
ChCode.Non = Non;
ChCode.bps = bps;
ChCode.g           = [1 0 1 1 0 1 1; 1 1 1 1 0 0 1]; %133 171
ChCode.nsc_flag    = 1;
ChCode.Nr_infobits = Payload_size*8 ;

% Turbo LTE
ChCode.Turbo_g = [1,0,1,1; 1,1,0,1]; %13 15 
ChCode.Turbo_nsc_flag = 0;
ChCode.Turbo_pun_pattern = [1; 1; 0; 0]; % rate 1/2;
ChCode.Turbo_tail_pattern = ones(4,3);
ChCode.Turbo_interleaver = CreateLTEInterleaver(ChCode.Nr_infobits);

%hSCR      = comm.Scrambler(2, [1 0 0 0 1 0 0 1],[1 0 1 1 1 0 1]);
%hDSCR     = comm.Descrambler(2, [1 0 0 0 1 0 0 1],[1 0 1 1 1 0 1]);
%scr_seq   = step(hSCR, zeros(Nr_infobits,1));
%scr_seq   = scr_seq';
ChCode.CodeType = CodeType;
rate = CodeType(end-3+1:end);
alg = ChCode.CodeType(1:2);
switch alg
    case 'CC'
        ChCode.Alg  = 0; % convolutional 
    case 'TC'
        ChCode.Alg  = 1; % turbo 
      otherwise
        warning('Unexpected Coding type.')
end
switch rate
    case '1-2'
        Nbit_per_block  = Non*bps/2;
    case '2-3'
        Nbit_per_block  = Non*bps/3*2;
    case '3-4'
        Nbit_per_block  = Non*bps/4*3;
    case '5-6'
        Nbit_per_block = Non*bps/6*5;
    otherwise
        warning('Unexpected Coding Rate.')
end
% 6: is the constraint number. I.e. number of shift registers.
cnst = size(ChCode.g,2)-1;
Nb  = ceil((ChCode.Nr_infobits+cnst)/Nbit_per_block);
ChCode.Nb = Nb;
% to fix numeric precesion error
ChCode.Nr_padbits = Nb*Nbit_per_block - ChCode.Nr_infobits;
Nr_unpuncbits     = (ChCode.Nr_infobits + ChCode.Nr_padbits)*2;
switch rate
    case '1-2'
        ChCode.punc_flg   = ones(1,Nr_unpuncbits);
    case '2-3'
        ChCode.punc_flg   = ones(1,Nr_unpuncbits);
        ChCode.punc_flg(4:4:end) = 0;
    case '3-4'
        ChCode.punc_flg   = ones(1,Nr_unpuncbits);
        ChCode.punc_flg(3:6:end) = 0;
        ChCode.punc_flg(5:6:end) = 0;
    case '5-6'
        ChCode.punc_flg   = ones(1,Nr_unpuncbits);
        ChCode.punc_flg(4:10:end) = 0;
        ChCode.punc_flg(5:10:end) = 0;
        ChCode.punc_flg(8:10:end) = 0;
        ChCode.punc_flg(9:10:end) = 0;
    otherwise
        warning('Unexpected Coding Scheme.')
end
Ncbpsym           = Non*bps;
ChCode.Ncbpsym = Ncbpsym;
rng(seed); % to keep the interleaver 
ChCode.interleaver1 = randperm(Ncbpsym);
ChCode.interleaver2 = randperm(Ncbpsym);
ChCode.decoder_type = 2;
ChCode.demod_type = 2;
ChCode.info_bits = zeros(ChCode.Nr_infobits,1);
ChCode.number_itr = 20;
 % the decoder type
%              = 0 For linear-log-MAP algorithm, i.e. correction function is a straght line.
%              = 1 For max-log-MAP algorithm (i.e. max*(x,y) = max(x,y) ), i.e. correction function = 0.
%              = 2 For Constant-log-MAP algorithm, i.e. correction function is a constant.
%              = 3 For log-MAP, correction factor from small nonuniform table and interpolation.
%              = 4 For log-MAP, correction factor uses C function calls.
end

function [ qam_data, info_bits, encoded_bits] = GenerateEncodedData( ChCode)
% CC Encoding
info_bits = randi([0 1],ChCode.Nr_infobits,1); % here we skip scrambling, coz we have randomly generated the source bit sequence.
[ qam_data, encoded_bits] = Encoder( ChCode, info_bits);
end
% ChCode: a structure contains information about CMS generated by
% 'GenerateCode.m' and frame length, etc.
% info bits
function [ qam_data, encoded_bits] = Encoder( ChCode, info_bits)
% CC Encoding
% info_bits: column vector
bps = ChCode.bps;
Nb = ChCode.Nb;
Non = ChCode.Non;
Ncbpsym = ChCode.Ncbpsym;
constellation =  ChCode.Constellation;

% zero padding
pad_bits  = randi([0 1],ChCode.Nr_padbits,1);  % the type of pad_bits does not matter here
if ChCode.Alg == 0 %CC encoder
    codeword_o   = ConvEncode(info_bits',ChCode.g,ChCode.nsc_flag);
elseif ChCode.Alg == 1 % turbo encoder
    g1 = ChCode.Turbo_g; %13 15
    g2 = g1;
    nsc_flag1 = ChCode.Turbo_nsc_flag ;
    nsc_flag2 = nsc_flag1;
    pun_pattern = ChCode.Turbo_pun_pattern; % rate 1/2;
    tail_pattern = ChCode.Turbo_tail_pattern;
    interleaver = ChCode.Turbo_interleaver;
    codeword_o = TurboEncode( info_bits', interleaver, pun_pattern, tail_pattern, g1, nsc_flag1, g2, nsc_flag2 );
end
codeword_pad = ConvEncode(pad_bits',ChCode.g,ChCode.nsc_flag); % padding bit are just encoded with covolutional code
% Puncture
encoded_bits = [codeword_o'; codeword_pad'];
encoded_bits = encoded_bits(ChCode.punc_flg == 1);

% Interleaving
encoded_mtx = reshape(encoded_bits,Ncbpsym,Nb);
encoded_mtx = encoded_mtx(ChCode.interleaver1,:);
encoded_bits = encoded_mtx(:);
% QAM Modulating
data_symb_idx = bi2de(reshape(encoded_bits,bps,length(encoded_bits)/bps)','left-msb') + 1;
data_symb_idx = reshape(data_symb_idx,Non, Nb);
qam_data      = constellation(data_symb_idx);
end
% it uses soft input i.e. LLRs of the symbols
% ChCode: a structure contains information about CMS generated by
% 'GenerateCode.m'
% data a vector of demodulated data, Esno is the Es/No ratio of each symbol
function [ bits_e, encoded_bits_e ] = Decoder( ChCode, data_e, Esno, h)

constellation = ChCode.Constellation;
Nb = ChCode.Nb;
Ncbpsym = ChCode.Ncbpsym;
% Symbol to bits LLR mapping
if nargin < 4
    LLRs_symb  = Demod2D(data_e.',constellation, Esno.');
else
    LLRs_symb  = Demod2D(data_e.',constellation, Esno.', h.');
end
%LLRs_symb = LLRs_symb -repmat(max(LLRs_symb,[],1),2^bps,1);
% convert symbols LLRs to bits LLRs

llrs_dem  = Somap(LLRs_symb, ChCode.demod_type);
encoded_bits_e = llrs_dem' > 0;
llrs_dem  = reshape(llrs_dem,Ncbpsym,Nb);
% Deinterleaving
llrs_dem(ChCode.interleaver1,:) = llrs_dem;
llrs_un_dec = zeros(1,length(ChCode.punc_flg));
llrs_un_dec(ChCode.punc_flg == 1) = llrs_dem;
cnst = size(ChCode.g,2)-1;

llrs_un_dec = llrs_un_dec(1:(ChCode.Nr_infobits+cnst)*2);
if ChCode.Alg == 0 %CC encoder
bits_e = ViterbiDecode(llrs_un_dec, ChCode.g, ChCode.nsc_flag);
elseif ChCode.Alg == 1 % turbo encoder
    g1 = ChCode.Turbo_g; %13 15
    g2 = g1; %13 15
    nsc_flag1 = ChCode.Turbo_nsc_flag ;
    nsc_flag2 = nsc_flag1;
    pun_pattern = ChCode.Turbo_pun_pattern; % rate 1/2;
    tail_pattern = ChCode.Turbo_tail_pattern;
    interleaver = ChCode.Turbo_interleaver;   
    bits_e = TurboDecode( llrs_un_dec, ChCode.info_bits.', ChCode.number_itr, ChCode.decoder_type,  interleaver, pun_pattern, tail_pattern, g1, nsc_flag1, g2, nsc_flag2 );      
end
bits_e = bits_e.';
end
function [ r, y ] = CRC( x, g )
% x: input bits index zero is MSB
% g: polynomial index zero is MSB
% r: reminder  index zero is MSB
% y: div index zero is MSB
K = length(x);
Q = length(g);

r = zeros(1,Q);
y = zeros(1,K); % div result is already in MSB first
% reverse g for indexing 
g = flip(g); 
for m=1:K
y(m) = xor(r(Q),x(m));     
for q=Q:-1:2
    r(q) = xor(r(q-1), and(g(q), y(m)));
end    
r(1) = and(g(1),y(m));
end
% to store in MSB 
r = flip(r);
end

%% --------------- END of Implementation ----------------------------------
end
