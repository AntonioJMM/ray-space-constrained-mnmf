function [] = parfor_rsmcnmf(eIdx,test)
load(['tmp_' test '.mat'])

beta = beta;
sigma = sigma;
psi = psi;

experimentName = experimentList{eIdx};
str = split(experimentName, '_');
fprintf(['(' num2str(eIdx) '/' num2str(length(experimentList)) ') ']);

checkComputed = [resultPath filesep experimentName, filesep, ...
    'beta',num2str(beta *10), filesep 'beta',num2str(beta *10) '_2_.wav'];
% if exist(checkComputed, 'file')
%     fprintf(['SKIP: ', experimentName, ' already computed!\n'])
% else
    tic
    referenceFolder = [referencePath, filesep, experimentName];
    fprintf(['Loading array signals from: ', ...
        referenceFolder, '\n']);
    arrayFile = [referenceFolder, filesep, 'array_mix.wav'];
    % Read the array signals
    [micSignal, fs] = audioread(arrayFile);
    micSignal = highpass(micSignal, fCutoff, fs);
    
    for ss = 1:sourceN
        referenceFile = [referenceFolder, filesep, ...
            'array_ref_', num2str(ss), '.wav'];
        referenceSignal{ss} =  audioread(referenceFile);
        referenceSignal{ss} = highpass(referenceSignal{ss}, ...
            fCutoff, fs);
    end
    % Compute the STFT of the mic signals
    for mm = 1:nMic
        [micSTFT(:,:,mm), fAx, tAx] = stft(micSignal(:,mm), ...
            analysisWin, hop, nfft, fs);
    end
    fAx = fAx(1:(nfft/2+1)).';  % Frequency axis
    tLen = length(tAx);         % Length of time axis
    fLen = nfft/2+1;            % Length of frequency axis
    
    % =========================================================================
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
    init.initW = 0.5 * (abs(randn(fLen,nBasis)) + ones(fLen,nBasis));% .* ...
    %     (psdMix * ones(1,nBasis));
    init.initH = 0.5 * (abs(randn(nBasis,tLen)) + ones(nBasis,tLen));
    % init.initQ = abs(init.initA).^2;
    % init.initQ = (0.5 * (1.9 * abs(randn(D*length(nu), sourceN)) + ...
    %     0.1 * ones(D*length(nu), sourceN))).^2;
    init.initQ = rand(size(greenFRaySpace,3), sourceN) + 1;
    init.initM = abs(greenFRaySpace);
    
    
    % Source separation using MCNMF in the ray space
    [estimateImage, Q, basisF, activationF, xRaySpace, invPsi, initQ, cost] =...
        rayspacenmf(micSTFT, mubar, D, nubar, sigma, fAx, psi, d, nMic, c, ...
        sourceN, nBasisSource, nIter, tik, init,beta);
    
    basisF = reshape(basisF,size(basisF,1),size(basisF,2)*size(basisF,3));
    activationF = reshape(permute(activationF,[2 1 3]),size(activationF,2),size(activationF,1)*size(activationF,3))';
    
    % Ray space reference and estimate
    for ss = 1:sourceN
        sourceRecSTFT{ss} = sqrt(basisF(:,sourceNMFidx{ss}) * ...
            activationF(sourceNMFidx{ss} ,:)) .* angle(micSTFT(:,:,1));
        sourceRec{ss} = istft(sourceRecSTFT{ss}, ...
            analysisWin, synthesisWin, hop, nfft, fs);
    end
    
    %% Estimated source image at the microphone
    istftParams.analysisWin = analysisWin;
    istftParams.synthesisWin = synthesisWin;
    istftParams.hop = hop;
    istftParams.nfft = nfft;
    
    % iRST of the estimated images
    rsmcnmfEstimate = arraysignalreconstruction(estimateImage, ...
        referenceSignal, invPsi, sourceN, nMic, istftParams, fs);
    
    %% Signal reconstruction
    saveDir = [resultPath filesep, experimentName, filesep, ...
        'beta', num2str(beta*10)];
    if ~exist(saveDir, 'dir')
        mkdir(saveDir)
    end
    fprintf(['Saving rayspace estimates to: ', saveDir, '\n']);
    
    for ss = 1:sourceN
        fileName = [saveDir, filesep, 'beta', num2str(beta*10),'_', ...
            num2str(ss), '_'];
        %             for mm = 1:nMic
        filePath = [fileName, '.wav'];
        audiowrite(filePath, 0.97*rsmcnmfEstimate{ss} ./ ...
            max(max(rsmcnmfEstimate{ss})), fs);
        %             end
    end
    toc
end
% end