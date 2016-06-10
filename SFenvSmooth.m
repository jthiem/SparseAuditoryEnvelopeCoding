function [ outenvs ] = SFenvSmooth( inenvs, envFilt )
%FBENVSMOOTHING Smooth the envelopes

% (c) Joachim Thiemann 2010
% for full license details see RunThesisCode.m
% and http://creativecommons.org/licenses/by/3.0/

[M,L] = size(inenvs);
outenvs = zeros(M,L);

for m = 1:M
    outenvs(m,:) = ...
        max(filtfilt(envFilt(m,:), 1, inenvs(m,:)),0);
end
