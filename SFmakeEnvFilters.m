function [ SAMP ] = SFmakeEnvFilters( FD )
%ENVFILTERS Make envelope smoothing filters

% (c) Joachim Thiemann 2010
% for full license details see RunThesisCode.m
% and http://creativecommons.org/licenses/by/3.0/

% shortcuts
M = FD.M;
b = FD.b;

% keep bw between 50Hz and 2kHz (experiment with this!)
bws = max(b,50);
bws = min(bws,2000);
bws = bws/FD.fs;
L = 512;

SAMP.ef = zeros(M,L);

for m = 1:M
    SAMP.ef(m,:) = ...
        firls( L-1, [ 0 bws(m)*.9 bws(m)*1.1 1 ], [ 1 1 0.1 0]);
    fprintf('.');
end
fprintf('\n');

SAMP.srate = floor(1./(2*bws));

end
