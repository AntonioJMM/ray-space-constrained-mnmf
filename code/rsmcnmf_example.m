function [] = rsmcnmf_example(test)

%% Setup
beta = 1.4;
% Setting up the directory of the experiments
experiment = 'rev_';    % Choose the experiment type
referencePath = ['../../../mpezzoli/NMF_RS_separation/review/reference_' test filesep];

% Setting the directory for storing the results
resultPath= ['./results_NEW_sinREST' test];   % The folder with the array mixtures

fprintf('Setup...\n');
% addpath(genpath('..'))

% =========================================================================
% General parameters
fs = 8000;                                          % Sampling frequency [Hz]
c = 340;                                            % Speed of sound [m/s]
sourceN = 2;                                        % Number of sources

% =========================================================================
% STFT parameters
winLength = 512;                                    % Window length in samples
analysisWin = hamming(winLength,'periodic');
synthesisWin = hamming(winLength,'periodic');
hop  = winLength / 4;                                % Hop size (75% overlap)
nfft = 2^nextpow2(winLength);
fmin = 0;                                           % Minimum frequency for processing
fmax = 4000;                                        % Maximum frequency for processing
fCutoff = 10;                                       % High pass the signal

% =========================================================================
% Dictionary parameters
fbinsep = fs/nfft;
fbins   = (0:round(nfft/2))*fs/nfft;
minbin  = round(fmin/fbinsep + 1);
maxbin  = round(fmax/fbinsep + 1);
kbins   = 2*pi*fbins ./ c;
mybins  = minbin:maxbin;
nbins   = maxbin-minbin+1;

% =========================================================================
% ULA parameters
d           = 0.03;      % Distance between two adjacent microphones
nMic        = 32;               % Number of microphones
micPos      = [zeros(nMic,1), (0:d:d*(nMic-1))', zeros(nMic,1)].';

% =========================================================================
% Point grid parameters
gridPx = 0.1:0.1:1;
gridPy = 0.1:0.1:1;
[XX, YY] = meshgrid(gridPx, gridPy);
gridPts = [XX(:), YY(:), zeros(length(gridPx)*length(gridPy),1)].';

% =========================================================================
% Ray Space Transform parameters
mubar = 0.06;                                       % m axis sampling interval
D = 4;                                              % Length of m axis
nubar = 4*d;                                        % nu axis sampling interval
sigma = 8*d;                                        % Gaussian window standard deviation
mu = ((0:mubar:(D-1)*mubar)-((D-1)/2*mubar))';      % [D,1] mu axis
nu = (0:nubar:micPos(2,end))';                      % [L,1] nu axis
t = 1:nMic;
<<<<<<< Updated upstream
tik = [6 10 7];                                     % Regularization parameter for the psiTilde
=======
tik = [6 10 7];                                          % Regularization parameter for the psiTilde
% load('tik.mat');
>>>>>>> Stashed changes

% =========================================================================
% MNMF parameters
nBasisSource = 12;
nBasis = nBasisSource * sourceN;
sourceNMFidx = cell(1,sourceN);
for iSrc = 1:sourceN
    sourceNMFidx{iSrc} = (1:nBasisSource) + (iSrc-1)*nBasisSource;
end
nIter = 400;
% beta = 1.5;

% =========================================================================
%% RST computation
psi = rayspacetransformmatrix(fbins(mybins),c,d,length(micPos),mubar,D,nubar,sigma);
I = size(psi, 1);                       % Number of Ray space data points

dist = pdist2(micPos', gridPts');        % Source microphones distances
for ss = 1:length(gridPts)
    % Microphone signals (e.g. free-field Green's Function)
    greenF(ss, :,:) = (exp(-1i*2*pi*fbins(mybins)/c .*dist(:,ss)) ./ (4*pi*dist(:,ss))).';
    % RST computation --> [greenFRaySpace] = [mybins, IxD, gridPts]
    greenFRaySpace(:,:,ss) = squeeze(spatialfilter(permute(greenF(ss, :,:), [2,1,3]), psi,false));
end

%% Experiments
allTest = dir([referencePath, filesep, experiment, '*']);
cnt = 1;
for tt = 1:length(allTest)
    if allTest(tt).name(1) ~= '.'
        experimentList{cnt} = [allTest(tt).name];
        cnt = cnt+1;
    end
end

experimentList = experimentList(1,[1 5 9 35 42 46 61]);
save(['tmp_' test '.mat'])

for eIdx = 1:length(experimentList)
    parfor_rsmcnmf(eIdx,test)
end

%***************************************************************
%***************************************************************
%***************************************************************
%***************************************************************
% for eIdx = 1:length(experimentList)
% experimentName = experimentList{eIdx};
% str = split(experimentName, '_');
% fprintf(['(' num2str(eIdx) '/' num2str(length(experimentList)) ') ']);
% 
% checkComputed = [resultPath filesep experimentName, filesep, ...
%     'beta',num2str(beta *10), filesep 'beta',num2str(beta *10) '_2_.wav'];
% % if exist(checkComputed, 'file')
% %     fprintf(['SKIP: ', experimentName, ' already computed!\n'])
% % else
%     tic
%     referenceFolder = [referencePath, filesep, experimentName];
%     fprintf(['Loading array signals from: ', ...
%         referenceFolder, '\n']);
%     arrayFile = [referenceFolder, filesep, 'array_mix.wav'];
%     % Read the array signals
%     [micSignal, fs] = audioread(arrayFile);
%     micSignal = highpass(micSignal, fCutoff, fs);
% 
%     for ss = 1:sourceN
%         referenceFile = [referenceFolder, filesep, ...
%             'array_ref_', num2str(ss), '.wav'];
%         referenceSignal{ss} =  audioread(referenceFile);
%         referenceSignal{ss} = highpass(referenceSignal{ss}, ...
%             fCutoff, fs);
%     end
%     % Compute the STFT of the mic signals
%     for mm = 1:nMic
%         [micSTFT(:,:,mm), fAx, tAx] = stft(micSignal(:,mm), ...
%             analysisWin, hop, nfft, fs);
%     end
%     fAx = fAx(1:(nfft/2+1)).';  % Frequency axis
%     tLen = length(tAx);         % Length of time axis
%     fLen = nfft/2+1;            % Length of frequency axis
% 
%     % =========================================================================
%     %% Source separation through RS-MCNMF
%     % Ray-Space-Transformed reference signals
%     raySpaceRef = cell(sourceN,1);
%     % Source signal estimate based on basis and activation functions
%     sourceRecSTFT = cell(sourceN,1);
%     sourceRec = cell(1, sourceN);
%     % Ray Space source images estimated by the RS-MCNMF
%     raySpaceEstimateImage = cell(sourceN,1);
% 
%     % Initialization of the RS-MCNMF algorithm
%     psdMix = 0.5 * (mean(abs(micSTFT(:,:,1)).^2 + abs(micSTFT(:,:,2)).^2, 2));
%     init.initA = 0.5 * (1.9 * abs(randn(nMic, sourceN)) + ...
%         0.1 * ones(nMic, sourceN));
%     % W is intialized so that its enegy follows mixture PSD
%     init.initW = 0.5 * (abs(randn(fLen,nBasis)) + ones(fLen,nBasis));% .* ...
%     %     (psdMix * ones(1,nBasis));
%     init.initH = 0.5 * (abs(randn(nBasis,tLen)) + ones(nBasis,tLen));
%     % init.initQ = abs(init.initA).^2;
%     % init.initQ = (0.5 * (1.9 * abs(randn(D*length(nu), sourceN)) + ...
%     %     0.1 * ones(D*length(nu), sourceN))).^2;
%     init.initQ = rand(size(greenFRaySpace,3), sourceN) + 1;
%     init.initM = abs(greenFRaySpace);
% 
% 
%     % Source separation using MCNMF in the ray space
%     [estimateImage, Q, basisF, activationF, xRaySpace, invPsi, initQ, cost] =...
%         rayspacenmf(micSTFT, mubar, D, nubar, sigma, fAx, psi, d, nMic, c, ...
%         sourceN, nBasisSource, nIter, tik, init,beta);
% 
%     basisF = reshape(basisF,size(basisF,1),size(basisF,2)*size(basisF,3));
%     activationF = reshape(permute(activationF,[2 1 3]),size(activationF,2),size(activationF,1)*size(activationF,3))';
% 
%     % Ray space reference and estimate
%     for ss = 1:sourceN
%         sourceRecSTFT{ss} = sqrt(basisF(:,sourceNMFidx{ss}) * ...
%             activationF(sourceNMFidx{ss} ,:)) .* angle(micSTFT(:,:,1));
%         sourceRec{ss} = istft(sourceRecSTFT{ss}, ...
%             analysisWin, synthesisWin, hop, nfft, fs);
%     end
% 
%     %% Estimated source image at the microphone
%     istftParams.analysisWin = analysisWin;
%     istftParams.synthesisWin = synthesisWin;
%     istftParams.hop = hop;
%     istftParams.nfft = nfft;
% 
%     % iRST of the estimated images
%     rsmcnmfEstimate = arraysignalreconstruction(estimateImage, ...
%         referenceSignal, invPsi, sourceN, nMic, istftParams, fs);
% 
%     %% Signal reconstruction
%     saveDir = [resultPath filesep, experimentName, filesep, ...
%         'beta', num2str(beta*10)];
%     if ~exist(saveDir, 'dir')
%         mkdir(saveDir)
%     end
%     fprintf(['Saving rayspace estimates to: ', saveDir, '\n']);
% 
%     for ss = 1:sourceN
%         fileName = [saveDir, filesep, 'beta', num2str(beta*10),'_', ...
%             num2str(ss), '_'];
%         %             for mm = 1:nMic
%         filePath = [fileName, '.wav'];
%         audiowrite(filePath, 0.97*rsmcnmfEstimate{ss} ./ ...
%             max(max(rsmcnmfEstimate{ss})), fs);
%         %             end
%     end
%     toc
% end
end