function [ corrfact, tqm ] = makeThreshold( FD, dbfs )
% Compute threshold of hearing values for envelopes
% assume I have FD

% (c) Joachim Thiemann 2010
% for full license details see RunThesisCode.m
% and http://creativecommons.org/licenses/by/3.0/

if nargin<2
    dbfs = 60;
end

% first calibrate input: ~1 kHz sine wave +/- 1
f = 984.17;
sin1k = cos(2*pi*f*(0:15999)/16000);

% filter through filterbank
c = abs(SFanalysis(sin1k,FD));
[amp,ch] = max(c(:,9000));
% value of 'amp' represents 60 dB perceived loudness.
fprintf('%f Hz in ch %d (fc = %f Hz) has amplitude %g, %f dB\n', ...
    f, ch, FD.fc(ch), amp, 20*log10(amp) );
corrfact = (10^(dbfs/20))/amp;
fprintf('Signal correction factor: %g\n', corrfact );

% now calculate threshold in quiet in dB
fck = FD.fc/1000;
tq = -6.5*exp(-0.6*(fck-3.3).^2) ...
    + 0.001*fck.^4 ...
    + 3.64*fck.^-0.8 ...
    - 80.64*exp(-4.712*fck.^0.5);
% and turn into envelope equiv
tqm = 10.^(tq/20);
fprintf('Thres. quiet at %f Hz: %g -> %f dB\n', ...
    FD.fc(ch), tqm(ch), tq(ch) );

