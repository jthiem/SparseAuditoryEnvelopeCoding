function [ env_f, env_x ] = SFsparse2Env( sampled, r, tqs, maxenv, srate )
%SFMAKEENVR based on envelope samples, compute fixed and max envs
%   sampled are the envelope samples
%   r is the correction factor
%   tqs is the threshold in quiet
%   maxenv are the transcoded pulses
%   srate is the per-channel sampling rate

% (c) Joachim Thiemann 2010
% for full license details see RunThesisCode.m
% and http://creativecommons.org/licenses/by/3.0/

[M,L] = size(sampled);
% get list of pulses
[mI,tI] = find(sampled>0);
Le = size(maxenv,2);

env_f = zeros(M,L);
env_x = zeros(M,L);

% find the peak of the transcoded env 
offsets = zeros(M,1);
for m=1:M
    [~, offsets(m)] = max(maxenv(m,:,m));
    env_x(m,:) = tqs(m);
end

K = length(mI);
for k=1:K
    m = mI(k);
    t = tI(k);
    a = sampled(m,t);
    if (t<offsets(m)) || (t>(L-offsets(m)))
        continue
    end
    
    sr2 = floor(srate(m)/2);
    trange = -sr2:sr2;
    env_f(m,trange+t) = a;

    trange = (1:Le)-offsets(m);
    env_x(:,trange+t) = max(env_x(:,trange+t),r*a*maxenv(:,:,m));
    if mod(k,1000)==0
        fprintf('.');
    end
end
fprintf('\n');
end

