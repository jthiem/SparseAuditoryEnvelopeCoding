function [ out ] = SFdelayComp( in, delay )
%SFDELAYCOMP Compensate for delay to make filters noncausal
%   Tries to be smart about dimensions

% (c) Joachim Thiemann 2010
% for full license details see RunThesisCode.m
% and http://creativecommons.org/licenses/by/3.0/

[M,N] = size(in);
if N>M
    in = in.';
    flip = 1;
    [M,N] = size(in);
else
    flip = 0;
end

if M<delay
    error('delay exceeds signal length!');
end

out = [ in(delay:M,1:N); zeros(delay-1,N) ];

if flip
    out = out.';
end

end

