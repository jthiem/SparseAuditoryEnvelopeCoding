function [ cout,Dout,Tout,n ] = BLitSparse( env_f, env_x, ch, Gf, wins, N, tau )
% Do iterative reconstruction of a short block of envelopes
% Inputs:
%   env_f   fixed portion of envelopes
%   env_x   envelope limits
%   ch      initial carrier estimate from runout of previous frame
%   Gf      gammatone filters IN FREQUENCY DOMAIN!
%   wins    the windows to handle overlap etc
%   N       maximum number of iterations to do
%   tau     early termination if error below this value

% (c) Joachim Thiemann 2010
% for full license details see RunThesisCode.m
% and http://creativecommons.org/licenses/by/3.0/

[M,ws] = size(env_f);
BS = size(Gf,2);
Hf = conj(Gf);

erange = 480:2368;

% cfix from previous ch
cfix = ch.*repmat(wins.prev,M,1);
win_c = repmat(wins.c,M,1);
% window the envelopes
envf_w = env_f.*repmat(wins.env,M,1);
envx_w = env_x.*repmat(wins.env,M,1);

Dout = zeros(N,1);
Tout = zeros(N,1);

fixed_mask = envf_w>0;
if sum(fixed_mask(:)) == 0
    cout = zeros(M,ws);
    n = 0;
    Dout = 0;
    Tout = 0;
    return
end
% loop
for n = 1:N
    % synthesize, then analyze
    xF = sum(fft(ch,BS,2).*Hf)+1e-80*randn(1,BS);
    cF = repmat(xF,M,1).*Gf;
    c = ifft(cF,BS,2);
    % apply envelopes and overlap
    c2 = cfix + c.*win_c;
    cenv = abs(c2);
    % do envelope correction
    cetemp = min(cenv,envx_w);               % handle max env
    cetemp(fixed_mask) = envf_w(fixed_mask); % handle fix env
    ch = cetemp .* c2./(cenv+1e-80);        % copy carrier
    % measure error
    em_cetemp = cetemp(:,erange).^.4;
    em_cenv = cenv(:,erange).^.4;
    Dout(n) = sum(abs(em_cetemp(:)-em_cenv(:)).^2);
    Tout(n) = sum(abs(em_cetemp(:).^2));
    if n>2 && (Dout(n)/Tout(n))<tau
        break
    end
end

cout = c2;

end
