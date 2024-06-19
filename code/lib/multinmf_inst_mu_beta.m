function [M, Q, W, H, cost] = multinmf_inst_mu_beta(V, beta, n_iter, Q, W, H, M, switch_Q, switch_W, switch_H)

% Multichannel NMF minimizing Itakura-Saito divergence through multiplicative updates
% with linear instantaneous mixing
%
% [Q, W, H, cost] = multinmf_inst_mu(V, n_iter, Q, W, H, part, switch_Q, switch_W, switch_H)
%
% Input:
%   - V: positive matrix data       (F x N x n_c)
%   - beta: the beta divergence value \in [0,2]
%   - n_iter: number of iterations
%   - init_Q: mixing matrix         (n_c x n_s)
%   - init_W: basis                 (F x K)
%   - init_H: activation coef.      (K x N)
%   - part : component indices
%   - switch_W, switch_H, switch_Q: (opt) switches (0 or 1) (def = 1)
%
% Output:
%   - Estimated Q, W and H
%   - Cost through iterations betw. data power and fitted variance.
%
%
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
%
% This code is based on
%
% A. Ozerov and C. Fevotte,
% "Multichannel nonnegative matrix factorization in convolutive mixtures for audio source separation,"
% IEEE Trans. on Audio, Speech and Lang. Proc. special issue on Signal Models and Representations
% of Musical and Environmental Sounds, vol. 18, no. 3, pp. 550-563, March 2010.
% Available: http://www.irisa.fr/metiss/ozerov/Publications/OzerovFevotte_IEEE_TASLP10.pdf
%
% The code has been adapted to the beta-divergence cost function by
% J. J. Carabias-Orti
%
% If you use this code please cite this paper
% M. Pezzoli, J. J. Carabias-Orti, M. Cobos, F. Antonacci,and A. Sarti.
% "Ray-space-based multichannel nonnegative matrix factorization for audio
% source separation". IEEE Signal Processing Letters (2021).
% doi: 10.1109/LSP.2021.3055463


if nargin < 8 || isempty(switch_Q)
    switch_Q = 1;
end

if nargin < 9 || isempty(switch_W)
    switch_W = 1;
end

if nargin < 10 || isempty(switch_H)
    switch_H = 1;
end

[F,N,n_c] = size(V);
n_s = size(Q,2);
lambda = 16;
<<<<<<< Updated upstream
% lambda = 0;
=======
>>>>>>> Stashed changes

W = reshape(W,F,size(W,2)/n_s,n_s);
H = permute(reshape(H',N,size(H,1)/n_s,n_s),[2 1 3]);

%% Definitions %%
cost    = zeros(1,n_iter);

%% Compute app. variance structure V_ap %%
P = tprod(W,[1 -1 3],H,[-1 2 3]);
MQ = tprod(Q,[-1 3],M,[1 2 -1]);
V_ap = tprod(MQ,[1 3 -1],P,[1 2 -1]);

V(V == 0) = eps;
if beta == 1
    cost(1) = sum(V(:).*log((V(:)./(V_ap(:)+eps))+eps)-V(:)+V_ap(:));
elseif beta==0
    cost(1) = sum((V(:)./(V_ap(:)+eps)) - log((V(:)./(V_ap(:)+eps))+eps)-1);
else
    cost(1) = 1/(beta*(beta-1)) * sum( V(:).^beta + ((beta-1)*V_ap(:).^beta) - (beta*V(:).*(V_ap(:).^(beta-1))) + lambda*sum(sum(Q'*Q - trace(Q'*Q))));
end

for iter = 2:n_iter
    tic
    %%% Update Q %%%
    if switch_Q
        Q_num = tprod(permute(V_ap.^(beta-2).*V,[3 1 2]),[-1 1 2],permute(M,[2 1 3]),[-1 1 3]);
        Q_num = tprod(P,[-1 -2 2],Q_num,[-1 -2 1]);
        Q_den = tprod(permute(V_ap.^(beta-1),[3 1 2]),[-1 1 2],permute(M,[2 1 3]),[-1 1 3]);
        Q_den = tprod(P,[-1 -2 2],Q_den,[-1 -2 1]);
        Q = Q.*((Q_num + lambda*Q)./(Q_den + lambda*ones(size(Q,1),size(Q,1))*Q + eps));

        Q(isnan(Q)) = 0;
        Q(~isfinite(Q)) = 0;

        MQ = tprod(Q,[-1 3],M,[1 2 -1]);
        V_ap = tprod(MQ,[1 3 -1],P,[1 2 -1]);
    end

    %%% Update W %%%
    if switch_W
        W_num = tprod(V_ap.^(beta-2).*V,[1 -1 4],H,[2 -1 3]);
        W_num = tprod(MQ,[1 -1 3],W_num,[1 2 3 -1]);
        W_den = tprod(V_ap.^(beta-1),[1 -1 4],H,[2 -1 3]);
        W_den = tprod(MQ,[1 -1 3],W_den,[1 2 3 -1]);
        W = W.*(W_num./(W_den+eps));

        W(isnan(W)) = 0;
        W(~isfinite(W)) = 0;

        P = tprod(W,[1 -1 3],H,[-1 2 3]);
        V_ap = tprod(MQ,[1 3 -1],P,[1 2 -1]);
    end

    %%% Update H %%%
    if switch_H
        H_num = tprod(V_ap.^(beta-2).*V,[1 2 -1],MQ,[1 -1 3]);
        H_num = tprod(W,[-1 1 3],H_num,[-1 2 3]);
        H_den = tprod(V_ap.^(beta-1),[1 2 -1],MQ,[1 -1 3]);
        H_den = tprod(W,[-1 1 3],H_den,[-1 2 3]);
        H = H.*(H_num./(H_den+eps));

        H(isnan(H)) = 0;
        H(~isfinite(H)) = 0;

        P = tprod(W,[1 -1 3],H,[-1 2 3]);
        V_ap = tprod(MQ,[1 3 -1],P,[1 2 -1]);
    end

    if beta == 1
        cost(iter) = sum(V(:).*log((V(:)./(V_ap(:)+eps))+eps)-V(:)+V_ap(:));
    elseif beta==0
        cost(iter) = sum((V(:)./(V_ap(:)+eps)) - log((V(:)./(V_ap(:)+eps))+eps)-1);
    else
        cost(iter) = 1/(beta*(beta-1)) * sum( V(:).^beta + ((beta-1)*V_ap(:).^beta) - (beta*V(:).*(V_ap(:).^(beta-1))) ) + lambda*sum(sum(Q'*Q - trace(Q'*Q)));
    end

    %%% Normalize %%%

%     %% Scale Q / W %%
%     if switch_Q && switch_W
%         scale = sum(Q,1);
%         Q = Q ./ repmat(scale,n_c,1);
%         W = W .* permute(repmat(scale,F,1,size(W,2)),[1 3 2]);
%     end
% 
%     %% Scale W / H %%
%     if switch_W && switch_H
%         scale = sum(W,1);
%         W = W ./ repmat(scale ,F,1,1);
%         H = H .* permute(repmat(scale,N,1,1),[2 1 3]);
%     end
    t = toc;
    fprintf('Iter: %3.0d -- Cost: %6.1f -- Time %.4f\n', iter, cost(iter), t);
end
fprintf('MU %d cost: %.4f -> %.4f\n', n_iter, cost(2), cost(iter));
