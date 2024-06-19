clc; clear all; close all;

fs = 8000;                                          % Sampling frequency [Hz]
c = 340;                                            % Speed of sound [m/s]
sourceN = 2;                                        % Number of sources
signalLen = 3;                                      % Length of the source signal
winLength = 256;                                    % Window length in samples
analysisWin = hamming(winLength,'periodic');
synthesisWin = hamming(winLength,'periodic');
hop = winLength / 4;                                % Hop size (75% overlap)
nfft = 2*2^nextpow2(winLength);
fmin = 0;                                           % Minimum frequency for processing
fmax = 4000;                                        % Maximum frequency for processing


fbinsep = fs/nfft;
fbins = (0:round(nfft/2))*fs/nfft;
minbin = round(fmin/fbinsep + 1);
maxbin = round(fmax/fbinsep + 1);
kbins = 2*pi*fbins ./ c;
mybins = minbin:maxbin;
nbins = maxbin-minbin+1;

nMic = 32;                                          % Number of microphones
micAperure = 3.2;
micPos = [zeros(nMic,1), linspace(0, micAperure, nMic).' zeros(nMic,1)].';
d = micPos(2,2);                                    % Distance between two adjacent microphones

mubar = 0.06;                                       % m axis sampling interval
D = 15;                                             % Length of m axis
nubar = 4*d;                                        % nu axis sampling interval
sigma = 16*d;                                       % Gaussian window standard deviation
mu = ((0:mubar:(D-1)*mubar)-((D-1)/2*mubar))';      % [D,1] mu axis
nu = (0:nubar:micPos(2,end))';                    % [L,1] nu axis
t = 1:nMic;
tik = [ 10.^(-7:-1) 0.1:0.1:1 1.5:0.5:5 6:2:20 20:10:100 10.^(3:7)];

psi = rayspacetransformmatrix(fbins(mybins),c,d,length(micPos),mubar,D,nubar,sigma);


fLen = nfft/2+1;            % Length of frequency axis

invPsi = zeros(size(psi,2),size(psi,1),fLen);
tik_f = zeros(1,fLen);
for ff = 7:170%fLen
    for jj = 1:length(tik)
        invPsi(:,:,ff,jj) = regularizedinversematrix(psi(:,:,ff), tik(jj));
        b(ff,jj) = real(sum(sum( psi(:,:,ff)*invPsi(:,:,ff,jj) - eye(32).*(psi(:,:,ff)*invPsi(:,:,ff,jj)))));
    end
    a = b(ff,:)<3;
    c = find(a==1);
    tik_f(ff) = tik(c(1));
end

tik_f

% figure; imagesc(abs(psi(:,:,2)*invPsi(:,:,2)))