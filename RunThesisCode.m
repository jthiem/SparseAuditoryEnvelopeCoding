% Run the sparse auditory envelope analysis and reconstruction
% as described in my thesis
% "A Sparse Auditory Envelope Representation with Iterative
% Reconstruction for Audio Coding"
% Ph.D. Thesis,  Joachim Thiemann, March 2011

% This code is (c) Joachim Thiemann 2011
% I make no gurantees that this code is correct, marketable,
% fit for any purpose whatsoever, or even related to the
% thesis at all.  The archive in which this script and the
% associated files are containd represent a best-effort to
% collate all the code used during the year of my research,
% but should only be regarded as supplementary information for
% the aid of future research.

% I am releasing this code under a Creative Commons CC-BY
% (Creative Commons Attribution 3.0 Unported License)
% license, the full terms of which can be found at
% http://creativecommons.org/licenses/by/3.0/
%
% This code will be posted on my blog, 
% http://signalsprocessed.blogspot.com/ and on
% http://www.ece.mcgill.ca/~jthiem1/research.html
% where the thesis can be found as well.

%% SETUP
%
% There are several data structures that are constant.  Here the 
% structures are computed before loading the audio files.  However, 
% the filterbank optimization as described in 5.1.1 requires the
% MATLAB optimisation toolbox: here, I just load the "Ms" vector
% (filter center frequencies in "ERB scale")
load Ms
% Using this vector the filterbank can be generated
FD = FBsetup( struct('Ms',Ms,'fs',16000,'L',1800) );
% Normalisation.  The procedure used in the thesis was a bit more
% elaborate, but this is almost the same: in the thesis code
% used to run the tests, max(sum(FD.G,1)) evaluates to 1.0060.
% The simple code below of course forces it to 1.
FD.G = FD.G/max(sum(FD.G,1));

% Now precompute the envelope smooting filters and sampling rates
% see section 5.1.2 of the thesis
EF = SFmakeEnvFilters( FD );

% precompute the transmultiplexed masking patterns T_k
% WARNING: the resulting array is about 121 MB.
[ ~, TF ] = makeTransLib( FD );

% precompute the threshold in quiet (5.1.3, eq. 5.5)
[ corrFact, TQ ] = makeThreshold(FD,80);
tqs = TQ/corrFact;

% clear things no longer needed
clear Ms corrFact TQ

% This sets the only variable parameter: the impact factor r_I
r = 1.0;

%% Analysis
%
% Load the audio file
% MAKE SURE THE FILE IS MONO AND SAMPLED AT 16kHz!
in = audioread( 'Databases/TIMIT/TIMIT/TEST/DR1/FAKS0/SA1_conv.wav' );

% to allow for delay and runout, extend by 1 sec
in = [ in; zeros(8000,1) ];

% intial analysis using the Basilar Membrane model (Gammatone filters)
c = SFanalysis( in, FD );

% envelope extraction (abs), smoothing, then sampling (5.1.2)
cs = SFenvSmooth( abs(c), EF.ef );
% sample envelopes
[M,L] = size(cs);
c_sampled = zeros(M,L);
for m = 1:M
    temp = 1:EF.srate(m):L;
    c_sampled(m,temp) = cs(m,temp);
end
fprintf( '%d samples in full envelope rep.\n', sum(c_sampled(:)>0) );

% remove samples due to threshold of hearing (5.1.3)
c_tq = c_sampled;
for m=1:M
    c_tq(m,:) = c_tq(m,:).*(c_tq(m,:)>(tqs(m)));
end
fprintf( '%d samples over threshold\n', sum(c_tq(:)>0) );

% apply masking model
REP = SFmask( c_tq, r, TF );
fprintf( '%d samples after masking\n', sum(REP(:)>0) );

% The array REP is the final sparse auditory envelope representation
% Efficiently encoding the representation into a bitstream is
% left as an exercise for the reader!

% don't need these anymore
clear c cs c_sampled c_tq temp

%% RESYNTHESIS
%
% Compute fixed regions (env_f) and envelope limits (env_x)
% see section 5.2.1, eqs. 5.7 - 5.9
[ env_f, env_x ] = SFsparse2Env( REP, r, tqs, TF, EF.srate );

% prepare for iterative block algorithm (section 5.2.3)
Dt = [];
Tt = [];
Na = [];

% for speed, we do filtering in frequency domain: precompute filters
BS=2048;
M = FD.M;
Gf = fft( FD.G, 2*BS, 2 );

adv = 80;
temp = hann(8*adv+1);
temp = temp(2:2:4*adv)';
wins.prev = [ ones(1,2*adv) 1-temp zeros(1,2*BS-4*adv) ];
wins.c = [ zeros(1,2*adv) temp ones(1,BS) 1-temp zeros(1,BS-6*adv) ]; 
wins.env = wins.prev+wins.c;

% do reconstruction loop
c_out = zeros(M,L);
prev = zeros(FD.M,2*BS);
F = floor((L-4*BS)/adv);
for f=0:F
    envs1 = env_f(:,(f*adv)+(1:2*BS));
    envs2 = env_x(:,(f*adv)+(1:2*BS));
    [ temp, D, T, Nout ] = ...
        BLitSparse( envs1, envs2, prev, Gf, wins, 20, 0.003 );
    c_out(:,(f*adv)+(1:adv)) = temp(:,adv+(1:adv));
    prev(:,1:(2*BS-adv)) = temp(:,(adv+1):(2*BS));
    
    if Nout>0
        fprintf('%d/%d %g (%d)\n',f,F,D(Nout)/T(Nout),Nout);
        Dt = [ Dt; D ];
        Tt = [ Tt; T ];
    end
    Na = [ Na; Nout ];
end
x_out = 2*real(SFsynthesis( c_out, FD ));
x_out = SFdelayComp( x_out, FD.L-adv );

