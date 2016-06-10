function [ out ] = SFsynthesis( c, FD )

% (c) Joachim Thiemann 2010
% for full license details see RunThesisCode.m
% and http://creativecommons.org/licenses/by/3.0/

[M,Lc] = size(c);
out = zeros(1,Lc);

H = fliplr(conj(FD.G));
for m = 1:M
    out = out+fftfilt(H(m,:),c(m,:));
end

end
