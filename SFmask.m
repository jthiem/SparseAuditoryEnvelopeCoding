function [ rep ] = SFmask( rep, r, te )
%SFRMASK mask by transmultiplexer
%   in: sparse input representation
%   r:  relaxation parameter
%   te: library of transmultiplexed pulses

% (c) Joachim Thiemann 2010
% for full license details see RunThesisCode.m
% and http://creativecommons.org/licenses/by/3.0/

% find maximal offsets
M = size(rep,1);
tem = zeros(M,1);
for m=1:M
    [~,tem(m)] = max(te(m,:,m));
end

% get list of pulses
[ mI, tI, P ] = find( rep );
Lt = size(te,2);
Ls = size(rep,2);

% sort by power
[ sP, pi ] = sort(P,'descend');
mI = mI(pi);
tI = tI(pi);

K = length(sP);
for k=1:K
    % sP(k) is the masker power, mI(k) the channel and tI(k) the time index
    % first check if masker still exists
    m = mI(k);
    t = tI(k);
    p = rep(m,t);     % should be same as sP(k) unless eliminated
    if p
        % now check if within bounds for comparison
        if (t>tem(m)) && (t<(Ls-Lt+tem(m)))
           masker = r*p*te(:,:,m);
           % time bounds
           lb = t-tem(m);
           ub = lb+Lt-1;
           block = rep(:,lb:ub);
           mask = block.*(block<masker);
           block(mask>0) = 0;
           rep(:,lb:ub) = block;
           rep(m,t) = p;   % put back the masking pulse!
        end
    end
    if mod(k,1000)==0
        fprintf('.');
    end
end
fprintf('\n');
end

