%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RS-MNMF example
% Ray-Space-Based Multichannel Nonnegative Matrix Factorization for Audio
% Source Separation.
% This script shows an example of audio source separation using the
% RS-MNMF algorithm.
%
% Copyright 2020 Mirco Pezzoli
% (mirco.pezzoli -at- polimi.it)
%
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
%
% If you use this code please cite this paper
% M. Pezzoli, J. J. Carabias-Orti, M. Cobos, F. Antonacci,and A. Sarti.
% "Ray-space-based multichannel nonnegative matrix factorization for audio
% source separation". IEEE Signal Processing Letters (2021).
% doi: 10.1109/LSP.2021.3055463

clear
close all
clc

%% Setup
fprintf('Setup...\n');
addpath(genpath('..'))

% General parameters
fs = 8000;                                          % Sampling frequency [Hz]
c = 340;                                            % Speed of sound [m/s]
sourceN = 3;                                        % Number of sources
signalLen = 3;                                      % Length of the source signal

% STFT parameters
winLength = 256;                                    % Window length in samples
analysisWin = hamming(winLength,'periodic');
synthesisWin = hamming(winLength,'periodic');
hop = winLength / 4;                                % Hop size (75% overlap)
nfft = 2*2^nextpow2(winLength);
fmin = 0;                                           % Minimum frequency for processing
fmax = 4000;                                        % Maximum frequency for processing

% Dictionary parameters
fbinsep = fs/nfft;
fbins = (0:round(nfft/2))*fs/nfft;
minbin = round(fmin/fbinsep + 1);
maxbin = round(fmax/fbinsep + 1);
kbins = 2*pi*fbins ./ c;
mybins = minbin:maxbin;
nbins = maxbin-minbin+1;

% ULA parameters
d = 0.03;                                           % Distance between two adjacent microphones
nMic = 32;                                          % Number of microphones
z = (0:d:d*(nMic-1))';                              % mic y coordinate
micPos = [zeros(nMic,1), z];                        % mic position in 2D

% Ray Space Transform parameters
mubar = 0.06;                                       % m axis sampling interval
D = 40;                                             % Length of m axis
nubar = 8*d;                                        % nu axis sampling interval
sigma = 16*d;                                       % Gaussian window standard deviation
mu = ((0:mubar:(D-1)*mubar)-((D-1)/2*mubar))';      % [D,1] mu axis
nu = (0:nubar:micPos(end,end))';                    % [L,1] nu axis
t = 1:nMic;
tik = 0.5;                                          % Regularization parameter for the psiTilde

% Point grid parameters
gridPx = [1.2 0.891 0.623 0.331];%linspace(-1,2,10);
gridPy = [1 0.783 0.442 0.086 -0.2];%linspace(0.5,5,5);
[XX, YY] = meshgrid(gridPx, gridPy);
gridPts = [XX(:), YY(:)]';

figure
scatter(micPos(:,1), micPos(:,2), 'filled')
hold on
scatter(gridPts(1,:), gridPts(2,:), 'x')
text(gridPts(1,:)+0.01, gridPts(2,:)+0.01, num2cell(1:length(gridPts)));
xlabel('x'); ylabel('y');
axis equal;
title('Setup');

% MNMF parameters
nBasisSource = 12;
nBasis = nBasisSource * sourceN;
sourceNMFidx = cell(1,sourceN);
for iSrc = 1:sourceN
    sourceNMFidx{iSrc} = (1:nBasisSource) + (iSrc-1)*nBasisSource;
end
nIter = 500;
beta = 0.9;

% Source parameters
sourceType = 'music';                                % male/female/music

% Directory where the impulse resposes are stored
rirPath = ['..' filesep 'data' filesep 'sources'];
sourceLocationFile = [rirPath, filesep, 'source_location.mat'];
sourcePath = ['..' filesep 'data' filesep 'audio' filesep, sourceType];

load(sourceLocationFile)                            % Load the location of the sources from disk
% Randomly choose N sources from the dataset
% sourceLabel = randperm(9,sourceN);
sourceLabel = [1 5 9];
sourcePosition = sourceLocation(sourceLabel, :);

%% Geometric setup
% Source location from office room experiment
figure(1)
scatter(micPos(:,1), micPos(:,2))
hold on
scatter(sourcePosition(:,1), sourcePosition(:,2), '^', 'filled')
text(sourcePosition(:,1)+0.05, sourcePosition(:,2), ...
    cellstr(strsplit(num2str(sourceLabel))))
title('Setup'), xlim([-0.1, 3.8]), grid on
legend('Mics', 'Sources')

%% RST computation
psi = rayspacetransformmatrix(fbins(mybins),c,d,length(micPos),mubar,D,nubar,sigma);
I = size(psi, 1);                       % Number of Ray space data points

dist = pdist2(micPos, gridPts');        % Source microphones distances
for ss = 1:length(gridPts)
    % Microphone signals (e.g. free-field Green's Function)
    greenF(ss, :,:) = (exp(-1i*2*pi*fbins(mybins)/c .*dist(:,ss)) ./ (4*pi*dist(:,ss))).';  
    % RST computation --> [greenFRaySpace] = [mybins, IxD, gridPts]
    greenFRaySpace(:,:,ss) = squeeze(spatialfilter(permute(greenF(ss, :,:), [2,1,3]), psi,false));
end

%% Microphone array signals
fprintf('Compute the array signal...\n');

referenceSTFT = cell(sourceN,1);
referenceSignal = referenceSTFT;
rirLength = 0.8;
rirSamp = rirLength*fs;
micSignal = zeros((signalLen+1)*fs, nMic);
sourceSignal = zeros((signalLen+1)*fs, sourceN);    % Source signal in time

for ss = 1:sourceN
    rirFile = [rirPath, filesep, 'rir_source_', ...
        num2str(sourceLabel(ss)), '.wav'];
    [rir, rirFs] = audioread(rirFile);
    sourceFile = [sourcePath, filesep, num2str(ss), '.wav'];
    [tmp, srcFs] = audioread(sourceFile);

    start = find(tmp > 0.8 * var(tmp), 1 );  % Threshold
    stop = start + srcFs * signalLen;
    tmp = tmp(start:stop-1, 1);
    tmp = resample(tmp, rirFs, srcFs);   % resampling the signal
    tmp = tmp ./max(abs(tmp(:)));


    for mm = 1:nMic
        impResp = rir(:,mm);
        impResp = impResp(1:rirSamp);               % rir shortening
        currentSignal = conv(tmp, impResp);
        currentSignal = resample(currentSignal, fs, rirFs);  % Resample signal
        zeroPad = (signalLen+1)*fs - length(currentSignal);
        currentSignal = [currentSignal; zeros(zeroPad,1)];
        micSignal(:,mm) = micSignal(:,mm) + currentSignal;
        referenceSignal{ss}(:,mm) = currentSignal;
        referenceSTFT{ss}(:,:,mm) = stft(currentSignal, analysisWin, ...
            hop, nfft, fs);
    end
    sourceTemp = resample(tmp,fs,rirFs);
    sourceSignal(:,ss) = [sourceTemp; zeros(fs, 1)];
end

for mm = 1:nMic
    [micSTFT(:,:,mm), fAx, tAx] = stft(micSignal(:,mm), analysisWin, ...
        hop, nfft, fs);
end
tLen = length(tAx);         % Length of time axis
fLen = nfft/2+1;            % Length of frequency axis

%% Source separation through RS-MCNMF

% Ray-Space-Transformed reference signals
raySpaceRef = cell(sourceN,1);
% Source signal estimate based on basis and activation functions
sourceRecSTFT = cell(sourceN,1);
sourceRec = cell(1, sourceN);
% Ray Space source images estimated by the RS-MCNMF
raySpaceEstimateImage = cell(sourceN,1);

% Initialization of the RS-MCNMF algorithm
psdMix = 0.5 * (mean(abs(micSTFT(:,:,1)).^2 + abs(micSTFT(:,:,2)).^2, 2));
init.initA = 0.5 * (1.9 * abs(randn(nMic, sourceN)) + ...
    0.1 * ones(nMic, sourceN));
% W is intialized so that its enegy follows mixture PSD
init.initW = 0.5 * (abs(randn(fLen,nBasis)) + ones(fLen,nBasis)) .* ...
    (psdMix * ones(1,nBasis));
init.initH = 0.5 * (abs(randn(nBasis,tLen)) + ones(nBasis,tLen));
% init.initQ = abs(init.initA).^2;
% init.initQ = (0.5 * (1.9 * abs(randn(D*length(nu), sourceN)) + ...
%     0.1 * ones(D*length(nu), sourceN))).^2;
init.initQ = rand(size(greenFRaySpace,3), sourceN) + 1;
init.initM = abs(greenFRaySpace);


% Source separation using MCNMF in the ray space
[estimateImage, Q, basisF, activationF, xRaySpace, psi, ...
    invPsi, initQ, cost] = rayspacenmf(micSTFT, mubar, D, nubar, sigma, fAx, ...
    d, nMic, c, sourceN, nBasisSource,...
    nIter, tik, init,beta);

basisF = reshape(basisF,size(basisF,1),size(basisF,2)*size(basisF,3));
activationF = reshape(permute(activationF,[2 1 3]),size(activationF,2),size(activationF,1)*size(activationF,3))';

% Ray space reference and estimate
for ss = 1:sourceN
    raySpaceRef{ss} = spatialfilter(referenceSTFT{ss},psi,false);
    raySpaceRef{ss} = permute(raySpaceRef{ss}, [3,2,1]);
    sourceRecSTFT{ss} = sqrt(basisF(:,sourceNMFidx{ss}) * ...
        activationF(sourceNMFidx{ss} ,:)) .* angle(micSTFT(:,:,1));
    sourceRec{ss} = istft(sourceRecSTFT{ss}, ...
        analysisWin, synthesisWin, hop, nfft, fs);
end

% Metrics for the estimation of the sources only
sourceEstimation = cell2mat(sourceRec.');
[sourceRaySDR, sourceRaySIR, sourceRaySAR, perm] = ...
    bss_eval_sources(sourceEstimation, sourceSignal.');

estimateImage = estimateImage(:,:,perm,:);
sourceRecSTFT = sourceRecSTFT(perm);
estimateImage = estimateImage(:,:,perm,:);
for ss = 1:sourceN
    raySpaceEstimateImage{ss} = squeeze(estimateImage(:,:,ss,:));
    raySpaceEstimateImage{ss} = permute(raySpaceEstimateImage{ss}, ...
        [3,2,1]);
end

%% Plot the Ray space data
fprintf('Show Ray Space separation results...\n');
fIdx = 18;      % Frequency index to show
% Ray Space mixture
mixRay = squeeze(abs(xRaySpace(fIdx ,:,:))).'.^2;

% Find color range
% Temporary variables
tmpRef = cat(3, raySpaceRef{:});
tmpRef = abs(tmpRef(:,:,fIdx));
tmpEst = cat(3, raySpaceEstimateImage{:});
tmpEst = abs(tmpEst(:,:,fIdx));
maxC = max([mixRay(:); tmpRef(:); tmpEst(:)]);
minC = min([mixRay(:); tmpRef(:); tmpEst(:)]);

figure(2)
imagesc(tAx, t, mixRay, [minC, maxC])
xlabel('Time [s]'), ylabel('Ray space point t');
colorbar
title(['Ray Space Mixture at ' num2str(fAx(fIdx)), ' [Hz]']);

figure(3)
if sourceN == 2
    nCol = 2;
else
    nCol = 3;
end
nRow = 2;
for ss = 1:sourceN
    refRay = squeeze(abs(raySpaceRef{ss}(:,:,fIdx))).^2;
    estRay = squeeze(abs(raySpaceEstimateImage{ss}(:,:,fIdx))).^2;

    subplot(nRow,nCol,ss)
    imagesc(tAx, t, refRay, [minC, maxC])
    xlabel('Time [s]'), ylabel('Ray space point t');
    colorbar
    title(['Ray Space of source ' num2str(ss) ' image']);

    subplot(nRow,nCol,ss+sourceN)
    imagesc(tAx, t, estRay, [minC, maxC])
    xlabel('Time [s]'), ylabel('Ray space point t');
    colorbar
    title(['Estimated Ray Space source '  num2str(ss) ' image']);
end
% subplot(2,2,3)
% imagesc(tAx, t, estRay1, [minC, maxC])
% xlabel('Time [s]'), ylabel('Ray space point t');
% colorbar
% title('Estimated Ray Space source 1 image');
%
% subplot(2,2,4)
% imagesc(tAx, t, estRay2, [minC, maxC])
% xlabel('Time [s]'), ylabel('Ray space point t');
% colorbar
% title('Estimated Ray Space source 2 image');

%% Estimated source image at the microphone
fprintf('RS-MCNMF estimated source image at the microphones...\n');
istftParams.analysisWin = analysisWin;
istftParams.synthesisWin = synthesisWin;
istftParams.hop = hop;
istftParams.nfft = nfft;

% iRST of the estimated images
rsmcnmfEstimate = arraysignalreconstruction(estimateImage, ...
    referenceSignal, invPsi, sourceN, nMic, istftParams, fs);

% Dummy variables for bss_eval
% Matrix of the reference signals
ref = cat(3, referenceSignal{:});
% Matrix of the estimates given by the RS-MCNMF
est = cat(3, rsmcnmfEstimate{:});

raySDR = zeros(nMic,sourceN);
raySIR = zeros(nMic,sourceN);
raySAR = zeros(nMic,sourceN);

for mm = 1:nMic
    [raySDR(mm,:), raySIR(mm,:), raySAR(mm,:)] = bss_eval_sources( ...
        squeeze(est(:,mm,:)).', squeeze(ref(:,mm,:)).');
end

rsmcnmfSDR = mean(raySDR);
rsmcnmfSIR = mean(raySIR);
rsmcnmfSAR = mean(raySAR);

%% Plot results and listen to the signals
figure
dataPlot = [rsmcnmfSAR; rsmcnmfSDR; rsmcnmfSIR];
metrics = categorical({'SAR', 'SDR', 'SIR'});
bar(metrics, dataPlot);
grid on, ylabel('[dB]');
title('Separation metrics');
% disp(['Average array SDR: ' num2str(rsmcnmfSDR), '[dB]']);
% disp(['Average array SIR: ' num2str(rsmcnmfSIR), '[dB]']);
% disp(['Average array SAR: ' num2str(rsmcnmfSAR), '[dB]']);
fprintf('Listen to the separation results...\n');

micIdx = 1;

for ss = 1:sourceN
    fprintf(['RS-MCNMF Estimation source ' num2str(ss) ' in location '...
        num2str(sourceLabel(ss)) '\n']);
    fprintf('Press any key to start\n');
    pause();
    soundsc(rsmcnmfEstimate{ss}(:,micIdx), fs)
    pause(length(rsmcnmfEstimate{ss}(:,micIdx))/fs);
end

