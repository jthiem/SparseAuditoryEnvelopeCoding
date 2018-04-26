function [ FD ] = FBsetup( FD )
%FBSETUP Setup filter data
%
%   Input is a structure that will be filled with gammatone filter
%   data.
%    Accepted fields:
%     fs: sampling frequency (if not present, assume 16k)
%     Ms: vector containgin ERB indecies for filters
%         If empty, assume Ms = 1.5, 2, 2.5, ..., 32
%     bwscale: filter bandwidths can be scaled (default 1)
%     gd: "group delay", the number of samples at which the peak
%         of the Gammatone envelope will be.  If 0, peaks will not
%         be aligned.  If empty, peak is automatic minimum.
%     L:  Length of filters. If 0 or not given, L is determined
%         automatically (and replaced).
%
%    Output fields:
%     fc: Filter center frequencies, derived from Ms
%     b:  Filter bandwitdhs, also from Ms
%     M:  length(Ms) == length(fc), number of filters
%
%     G:  fixed-length gammatone filters, MxL array
%     H:  time-reversed complex conjugate of G
%
%     Gvar: variable-length gammatone filters, M-element struct of vectors
%     Hvar: time reversed complex conjugate of Gvar
%     Gvd(m):  delay that needs to be added to filter Gvar{m} to align
%              peaks
%     Hvd(m): ditto for Hvar{m}

% (c) Joachim Thiemann 2010
% for full license details see RunThesisCode.m
% and http://creativecommons.org/licenses/by/3.0/

%%

% defaults
if nargin==0
    FD = struct();
end
if ~isfield(FD, 'Ms')
    FD.Ms = 1:.5:32;
end
if ~isfield(FD, 'fs')
    FD.fs = 16000;
end
if ~isfield(FD, 'L')
    FD.L = 0;
end

% Calculate the center frequencies
FD.fc = (1000/4.37) * (10.^(FD.Ms./21.4) - 1)';
FD.M = length(FD.Ms);

% Sanity check: Make sure largest fc is less than Nyquist
if FD.fc(FD.M) > (FD.fs/2)
    error('Highest freq filer is > Nyquist!');
end

% precompute b
FD.b = 1.019 * ERB( FD.fc );

% Peak alignment
if ~isfield(FD,'gd')
    td = 3/(2*pi*FD.b(1));
    FD.gd = ceil(td*FD.fs);
    %fprintf('Using filter peak offset at %5d samples (%f s)\n', FD.gd, td );
end

% Now build the filterbank
tl = 12/(FD.b(1)*pi);
Ltemp = ceil(tl*FD.fs);

if FD.L == 0
    tl = 12/(FD.b(1)*pi);
    FD.L = Ltemp;
    fprintf('Using max filter length %5d samples (%f s)\n', FD.L, tl );
else
    Ltemp = FD.L;
end        
FD.G = zeros(FD.M,FD.L);

for m = 1:FD.M
    temp = makeFiltFixed( m, Ltemp, FD );
    FD.G(m,:) = temp;
end
%fprintf( 'Done\n' );

end

%% Fixed length filter
function out = makeFiltFixed( m, L, FD )

if FD.gd>0
    n = (-FD.gd+1):(L-FD.gd);
    td = 3/(2*pi*FD.b(m));
else
    % traditional GT (not aligned)
    n = 1:L;
    td = 0;
end

t = n/FD.fs;
% sanity check
if (t(1)+td-1/FD.fs)>0
    fprintf('!');
end

% envelope
out = ((2*pi*FD.b(m)).^4)/6 * (t+td).^3 .* exp(-2*pi*FD.b(m)*(t+td));
out = out .* (out>0);  % half-wave rectify

% modulation
out = out .* exp(-2*pi*1i*FD.fc(m).*t);

end

%% ERB
function erb = ERB( f )
erb = 0.1079.*f + 24.7;
end
