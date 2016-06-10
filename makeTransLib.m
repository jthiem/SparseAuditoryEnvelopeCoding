function [ minR, maxEnv ] = makeTransLib( FD )
%MAKETRANSLIB mae a library of transmultiplexed masking patterns

% (c) Joachim Thiemann 2010
% for full license details see RunThesisCode.m
% and http://creativecommons.org/licenses/by/3.0/

[M,L] = size(FD.G);
Lt = 2*L-1;
maxEnv = zeros(M,Lt,M);
minR = zeros(M,2);

for m1=1:M
    fprintf('.');
    sig = zeros(1,Lt);
    % create transcoded pulse
    sig(1:L) = fliplr(FD.G(m1,1:L))';
    maxEnv(:,:,m1) = abs(SFanalysis( sig, FD ));
    % normalise to +1dB
    scale = max(maxEnv(m1,:,m1));
    maxEnv(:,:,m1) = (10^(1/20)/scale)*maxEnv(:,:,m1);
    % find -1dB points
    %      range = find(maxEnv(m1,:,m1)>(10^(-1/20)));
    range = find(maxEnv(m1,:,m1)>1);
    minR(m1,1) = min(range);
    minR(m1,2) = max(range);
end

fprintf('\n');
end
