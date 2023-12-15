function Im = multinmf_recons_im(X,M,Q,W,H,part)

%
% Reconstructs source images from multi-nmf factorization (conservative)
%
% Im = multinmf_recons_im(X,Q,W,H,part)
%
%
% Input:
%   - X: truncated STFT of multichannal mixture (F x N x n_c)
%   - Q: squared filters                        (F x n_c x n_s)
%   - W: basis                                  (F x K)
%   - H: activation coef.                       (K x N)
%   - part : component indices
%
% Output:
%   - Im : reconstructed source images (F x N x n_s x n_c)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2010 Cedric Fevotte
% (cedric.fevotte -at- telecom-paristech.fr)
%
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
%
% If you use this code please cite this paper
%
% A. Ozerov and C. Fevotte,
% "Multichannel nonnegative matrix factorization in convolutive mixtures for audio source separation,"
% IEEE Trans. on Audio, Speech and Lang. Proc. special issue on Signal Models and Representations
% of Musical and Environmental Sounds, vol. 18, no. 3, pp. 550-563, March 2010.
% Available: http://www.irisa.fr/metiss/ozerov/Publications/OzerovFevotte_IEEE_TASLP10.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_s = size(Q,2);

P = tprod(W,[1 -1 3],H,[-1 2 3]);
MQ = tprod(Q,[-1 3],M,[1 2 -1]);
MQP = tprod(MQ,[1 3 4],P,[1 2 4]);

Im = permute(MQP ./ repmat(sum(MQP,4),1,1,1,n_s) .* repmat(X,1,1,1,n_s),[1 2 4 3]);
end
