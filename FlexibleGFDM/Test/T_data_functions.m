func = Data_functions();

Mc = 16; 
Non = 100000;

di = func.QAM_gen(Mc, Non);
dqam = func.QAM_mod(Mc, di);
N0 = 10;

x = dqam + sqrt(N0/2)*(randn(Non,1) + 1j*randn(Non,1));

Es = 1;
SNRl  = Es/N0;
die = func.QAM_demod(Mc, x);
ser = func.Ser(di, die);
snr = 1/func.Nmse(dqam, x);
seran = func.QAM_SER_an(Mc, SNRl);

diff = norm(seran- ser);
disp(['SER vs SER_an: ' num2str(diff)]);
%% Test CRC
g = hexToBinaryVector('1021',16, 'LSBFirst');
xbyte = {'E1', '00', 'CA', 'FE'};
Nbit = 8 * length (xbyte);
x = zeros(Nbit, 1);
for nb = 1: length(xbyte)    
x((nb-1)*8+(1:8)) = hexToBinaryVector(xbyte{nb}, 8, 'LSBFirst');
end

[r1, y] = func.CRC(x, g);
r = binaryVectorToHex(r1, 'LSBFirst');

%% Test encoding

Mc = 4;
Non = 128;
Payload_size = 10;
CodeType = 'CC-1-2';
ChCode = GenerateCode( CodeType,Payload_size,Non,Mc );


g1 = [1,0,1,1; 1,1,0,1];
g2 = g1;
nsc_flag1 = 0;
nsc_flag2 = 0;

pun_pattern = [1; 1; 0; 0];
tail_pattern = ones(4,3);
      
K = 1024;
%code_interleaver = CreateLTEInterleaver(K) ;

data = randi([0 1],1, ChCode.Nr_infobits);

ChCode.interleaver1 = randperm(ChCode.Nr_infobits);
encoded_bits = TurboEncode( data, ChCode.interleaver1, pun_pattern, tail_pattern, g1, nsc_flag1, g2, nsc_flag2 );
bps = ChCode.bps;
Nb = ChCode.Nb;
Non = ChCode.Non;
Ncbpsym = ChCode.Ncbpsym;
constellation =  ChCode.Constellation;

data_symb_idx = bi2de(reshape(encoded_bits,bps,length(encoded_bits)/bps)','left-msb') + 1;
qam_data      = constellation(data_symb_idx);

a = 2;
EsN0 = 10;
N0 = 1/EsN0;
qam_data_e =  a*qam_data + GenRandomNoise(size(qam_data), N0);
Esno = repmat(EsN0, numel(qam_data),1);
fad = repmat(a, numel(qam_data),1);
LLRs_symb  = Demod2D(qam_data_e.',constellation, Esno.', fad.');
%LLRs_symb = LLRs_symb -repmat(max(LLRs_symb,[],1),2^bps,1);
% convert symbols LLRs to bits LLRs
llrs_dem  = Somap(LLRs_symb, 3);

[detected_data, errors, output_decoder_c ] = TurboDecode( llrs_dem, data, 14, 3,  ChCode.interleaver1, pun_pattern, tail_pattern, g1, nsc_flag1, g2, nsc_flag2 );      
%% CRC
g = [0,0,1,1,0,0,0,1];
x = [1 0 0 0 1 1 1 1 1 0 0 0 1 1 1 1];

[r, y] = CRC(x, g);

Q = length(g);
Ug = toint([1,g])
Ux = toint([x,zeros(1,Q)])
Ux1 = toint([x, zeros(1,Q)])

Ur = toint(r)
Uy = toint(y)

dm = floor(Ux/Ug)
rm = mod(Ux,Ug)
r1 = mod(Ux,Uy)
function val = toint(x)
% MSB first 
N = length(x);
x = reshape(x, 1,N);
val = sum(x.*2.^(N-1:-1:0));
end

%

