function [ c ] = SFanalysis( in, FD )

% (c) Joachim Thiemann 2010
% for full license details see RunThesisCode.m
% and http://creativecommons.org/licenses/by/3.0/

[M,N] = size(FD.G);
Lt = length(in);
c = zeros(M,Lt);

G = FD.G;
for m = 1:M
    c(m,:) = fftfilt( G(m,:), in );
end

end
