%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Evaluation of source separation algorithms                            %
%   Analyze the results obtain for different source setups.               %
%   v 0.2                                                                 %
%   Mirco Pezzoli                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc
addpath(genpath('lib'));

%% Setup
fprintf('Setup...\n');
homeFolder  = '../../../mpezzoli/NMF_RS_separation/review/';
homeFolder2 = './';
sourceSigType = 'music';
referencePath = [homeFolder 'reference_' sourceSigType];    % Path with the reference
resultsPath = [homeFolder2 'results_NEW_sinREST' sourceSigType];   % Path with estimates of the BSS
metricsPath = [homeFolder2 'metrics_NEW_sinREST' sourceSigType];%'metrics_music3';



% Algoritms
algorithms = {'beta0','beta1','beta2','beta3','beta4','beta5','beta6','beta7','beta8','beta9','beta10','beta11','beta12','beta13','beta14','beta15','beta16','beta17','beta18','beta19','beta20'};
algorithms = {'beta14'};
nAlgorithms = length(algorithms);
computeMetrics = true*ones(1,nAlgorithms);
% computeMetrics(4) = true;
% Source settings
nSource = 2;
idxExp = 1;
% Mic setting
nMic = 32;

if idxExp == 1
    experiment = 'rev_';
else
    experiment = 'test_';
end
if strcmp(experiment, 'test_')
    sourcePos{1} = [0.5 0.75 1];% 1.5 2];
    sourcePos{2} = [0.5 0.75 1];% 1.5 2];
else
    sourcePos{1} = [0.3 0.6 0.9]; %1.2 1.5 1.8 2.1];
    sourcePos{2} = [0.3 0.6 0.9]; %1.2 1.5 1.8 2.1];
end

[A,B] = meshgrid(sourcePos{2},sourcePos{1});
posCombination = cat(2, A', B');
posCombination = reshape(posCombination,[],2);
posCombination = fliplr(posCombination);

allExperiment = dir(resultsPath);
cnt = 1;
for tt = 1:length(allExperiment)
    if allExperiment(tt).name(1) ~= '.'
        experimentName = allExperiment(tt).name;
        str = split(experimentName, '_');
        if ~strcmp([str{1} '_'], experiment)
            continue
        else
            experimentList{cnt} = [allExperiment(tt).name];
            cnt = cnt+1;
        end
    end
end

%% Load data
dataSAR = cell(nAlgorithms, 1);%zeros(length(sourcePos{1})*3, length(sourcePos{2})*3);
dataSDR = cell(nAlgorithms, 1);%zeros(length(sourcePos{1})*3, length(sourcePos{2})*3);
dataSIR = cell(nAlgorithms, 1);%zeros(length(sourcePos{1})*3, length(sourcePos{2})*3);

for aa = 1:nAlgorithms
    dataSAR{aa} = zeros(length(sourcePos{1})*3, length(sourcePos{2})*3);
    dataSDR{aa} = zeros(length(sourcePos{1})*3, length(sourcePos{2})*3);
    dataSIR{aa} = zeros(length(sourcePos{1})*3, length(sourcePos{2})*3);
end

count = zeros(length(sourcePos{1})*3, length(sourcePos{2})*3);

for cc = 1:1%size(posCombination,1)
    src1Name = num2str(posCombination(cc,1) * 100);
    src2Name = num2str(posCombination(cc,2) * 100);
    for tt = 1:length(experimentList)
        str = split(experimentList{tt}, '_');
        %         if length(str) < 3
        %             continue
        %         end
        %         if ~((strcmp(str{2}, src1Name)) && (strcmp(str{3}, src2Name)))
        %             continue
        %         end
        %         if ~strcmp(str{1}, strrep(experiment, '_',  ''))
        %             continue
        %         end

        idx1 = round(str2num(str{4}) / 10);
        idx2 = mod(str2num(str{4}), 10);
        %         if idx1 > idx2
        %             continue
        %         end
        %% Load reference signals
        if computeMetrics == true
            % Directory of the experiment

        end
        %% Load estimed signals
        % Directory of the experiment
        experimentName = experimentList{tt};
        estimateExperimentPath = [resultsPath, filesep, experimentName];
        fprintf(['Loading estimates from: ' estimateExperimentPath '\n'])

        referenceExperimentPath = [referencePath, filesep, experimentName];


        if sum(computeMetrics) > 0
            fprintf(['Loading reference from: ' referenceExperimentPath '\n'])
            array_ref = cell(1,nSource);
            for ss = 1:nSource
                fileName = [referenceExperimentPath filesep, 'array_ref_' num2str(ss), '.wav'];
                [array_ref{ss}, fs] = audioread(fileName);
            end
            ref = permute(cat(3, array_ref{:}), [3,1,2]);
        end
        for aa = 1:nAlgorithms
            if computeMetrics(aa) == true
                dirPath = [estimateExperimentPath, filesep, algorithms{aa}];
                for ss = 1:nSource
                    %                     if strcmp(algorithms{aa}, 'beamSpace') || ...
                    %                             strcmp(algorithms{aa}, 'multiCh')
                    %                         filePath = [dirPath filesep algorithms{aa} ...
                    %                             num2str(ss), '_'];
                    %                     else
                    filePath = [dirPath filesep algorithms{aa} '_' ...
                        num2str(ss), '_'];
                    %                     end
                    %                     for mm = 1:nMic
                    try
                        %                             fileName = [filePath, num2str(mm) ,'.wav'];
                        fileName = [filePath, '.wav'];
                        %                             [signalEstimate{aa}(ss,:,mm), fs] = audioread(fileName);
                        [signalEstimate{aa}(ss,:,:), fs] = audioread(fileName);
                    catch
                        fprintf(['No file found for ' dirPath, '\n']);
                        signalEstimate{aa} = [];
                    end
                    %                     end
                end
            end
        end

        %% Evaluate separation performace
        %         fprintf(['BSS metrics computation...\n'])
        for aa = 1:nAlgorithms
            if computeMetrics(aa) == true

                % Save metrics on disk
                metricsFile = [algorithms{aa},'_metrics.mat'];
                thisPath = [metricsPath, filesep, experimentName];
                if ~exist(thisPath)
                    mkdir(thisPath)
                end
                metricsFile= [thisPath, filesep, metricsFile];
                fprintf(['Saving metrics on disk: ' metricsFile, '\n']);

                %if ~exist(metricsFile)
                 %   disp("exist")

                    est = signalEstimate{aa};
                    for mm = 1:nMic
                        if isempty(est)
                            SDR(mm,:) = NaN;
                            SIR(mm,:) = NaN;
                            SAR(mm,:) = NaN;
                            PERM(mm,:) = NaN;
                        else
                            refMic = ref(:,:,mm);
                            estMic = est(:,:,mm);
                            diffLen = length(refMic) - length(estMic);
                            if diffLen > 0
                                estMic = [estMic, zeros(nSource, diffLen)];
                            end

                            estMic = estMic/max(abs(estMic(:)))*max(abs(refMic(:)));

                            [sdr, sir, sar, perm] = bss_eval_sources(estMic, refMic);
                            SDR(mm,:) = sdr(perm);
                            SIR(mm,:) = sir(perm);
                            SAR(mm,:) = sar(perm);
                            PERM(mm,:) = perm;
                        end
                    end
                    avgSDR{aa} = mean(SDR);
                    avgSIR{aa} = mean(SIR);
                    avgSAR{aa} = mean(SAR);


                    save(metricsFile, 'SAR', 'SDR', 'SIR', 'PERM')
              %  end
            % else
            %     metricsFile = [algorithms{aa},'_metrics.mat'];
            %     thisPath = [metricsPath, filesep, experimentName];
            %     if ~exist(thisPath)
            %         mkdir(thisPath)
            %     end
            %     metricsFile= [thisPath, filesep, metricsFile];
            %     %                 fprintf(['Loading metrics from disk: ' metricsFile, '\n']);
            %     load(metricsFile);
            %     avgSDR{aa} = nanmean(SDR);
            %     avgSIR{aa} = nanmean(SIR);
            %     avgSAR{aa} = nanmean(SAR);
            % end
            % if isnan(avgSAR{aa}(1)) || isnan(avgSAR{aa}(2))
            %     disp(['NaN found for ' algorithms{aa}]);
            %     missed{tt} = experimentName;
            % end
            end
        end

        %         %% Store the data for the full evaluation
        %         s1 = sourcePos{1} == posCombination(cc,1);
        %         s1idx = 1:numel(sourcePos{1});
        %         s1idx = 3*(s1idx(s1)-1)+idx1;
        %         s2 = sourcePos{2} == posCombination(cc,2);
        %         s2idx = 1:numel(sourcePos{2});
        %         s2idx = 3*(s2idx(s2)-1)+idx2;
        %         for aa = 1:nAlgorithms
        %             % All the SAR results
        %             dataSAR{aa}(s1idx, s2idx) =  dataSAR{aa}(s1idx, s2idx) + ...
        %                 avgSAR{aa}(1);
        %             dataSAR{aa}(s2idx, s1idx) =  dataSAR{aa}(s2idx, s1idx) + ...
        %                 avgSAR{aa}(2);
        %             % All the SDR results
        %             dataSDR{aa}(s1idx, s2idx) =  dataSDR{aa}(s1idx, s2idx) + ...
        %                 avgSDR{aa}(1);
        %             dataSDR{aa}(s2idx, s1idx) =  dataSDR{aa}(s2idx, s1idx) + ...
        %                 avgSDR{aa}(2);
        %             % All the SIR results
        %             dataSIR{aa}(s1idx, s2idx) =  dataSIR{aa}(s1idx, s2idx) + ...
        %                 avgSIR{aa}(1);
        %             dataSIR{aa}(s2idx, s1idx) =  dataSIR{aa}(s2idx, s1idx) + ...
        %                 avgSIR{aa}(2);
        %         end
        %         count(s1idx, s2idx) = count(s1idx, s2idx) +1;
        %         count(s2idx, s1idx) = count(s2idx, s1idx) +1;
        fprintf('\n');
    end
end
%% Create the matrices
noZero = count ~= 0;
for aa = 1:nAlgorithms
    dataSAR{aa}(noZero) = dataSAR{aa}(noZero) ./ count(noZero);
    dataSDR{aa}(noZero) = dataSDR{aa}(noZero) ./ count(noZero);
    dataSIR{aa}(noZero) = dataSIR{aa}(noZero) ./ count(noZero);

    dataSAR{aa}(~noZero) = NaN;
    dataSDR{aa}(~noZero) = NaN;
    dataSIR{aa}(~noZero) = NaN;
end

%% Plot setup
titleStr = 'Setup';
allPosition = unique(posCombination);

micPos = [0:0.03:(31*0.03); zeros(1,32)].';
iFig = 1;
figure(iFig), iFig = iFig +1;
scatter(micPos(:,1), micPos(:,2), 20, 'filled')
hold on
srcPos = [];
xPos = [0.1 0.45 0.81];

for ii = 1:length(allPosition)
    yPos = allPosition(ii) * ones(3, 1);
    srcPos = [srcPos; [xPos.', yPos]];
end
labels = string(1:length(srcPos));
scatter(srcPos(:,1), srcPos(:,2),  40, 'filled', 'diamond');
text(srcPos(:,1),srcPos(:,2),labels,'VerticalAlignment','top', ...
    'HorizontalAlignment','left')
xlim([-0.1, 1]), ylim([-0.1, 2.2]), grid on
title(titleStr)
legend({'Array', 'Source'}, 'location', 'best')


saveDir = [metricsPath, filesep, 'images'];
if ~exist(saveDir, 'dir')
    mkdir([saveDir])
end

name = [experiment, 'setup'];
path = [saveDir, name];
saveas(gcf, path, 'epsc');

%% Show the metrics
fprintf('Plot the graphs with the metrics...\n');
allSAR = cat(3, dataSAR{:});
allSDR = cat(3, dataSDR{:});
allSIR = cat(3, dataSIR{:});

plotMatrix = false;
if plotMatrix == true
    for aa = 1:nAlgorithms
        % SAR PLOT
        minC = nanmin(allSAR(:));
        maxC = nanmax(allSAR(:));

        figure(iFig), iFig = iFig +1;
        b = imagesc(allSAR(:,:,aa), [minC, maxC]);
        axis equal, axis tight, colorbar
        set(gca,'xtick',[1:length(srcPos)],'xticklabel',labels);
        set(gca,'ytick',[1:length(srcPos)],'yticklabel',labels);
        cmap = colormap('parula');
        set(b, 'AlphaData', ~isnan(allSAR(:,:,aa)))
        title([algorithms{aa}, ' SAR']);

        name = [filesep experiment 'SAR_matrix_', algorithms{aa}];
        path = [saveDir, name];
        %     fprintf(['Saving plot on disk ' path '\n']);
        saveas(gcf, path, 'epsc');

        % SDR PLOT
        minC = min(allSDR(:));
        maxC = max(allSDR(:));

        figure(iFig), iFig = iFig +1;
        b = imagesc(allSDR(:,:,aa), [minC, maxC]);
        axis equal, axis tight, colorbar
        set(gca,'xtick',[1:length(srcPos)],'xticklabel',labels);
        set(gca,'ytick',[1:length(srcPos)],'yticklabel',labels);
        cmap = colormap('parula');
        set(b, 'AlphaData', ~isnan(allSDR(:,:,aa)))
        title([algorithms{aa}, ' SDR']);

        name = [filesep experiment 'SDR_matrix_', algorithms{aa}];
        path = [saveDir, name];
        %     fprintf(['Saving plot on disk ' path '\n']);
        saveas(gcf, path, 'epsc');

        % SIR PLOT
        minC = min(allSIR(:));
        maxC = max(allSIR(:));

        figure(iFig), iFig = iFig +1;
        b = imagesc(allSIR(:,:,aa), [minC, maxC]);
        axis equal, axis tight, colorbar
        set(gca,'xtick',[1:length(srcPos)],'xticklabel',labels);
        set(gca,'ytick',[1:length(srcPos)],'yticklabel',labels);
        cmap = colormap('parula');
        set(b, 'AlphaData', ~isnan(allSIR(:,:,aa)))
        title([algorithms{aa}, ' SIR']);

        name = [filesep experiment 'SIR_matrix_', algorithms{aa}];
        path = [saveDir, name];
        fprintf(['Saving plot on disk ' path '\n']);
        saveas(gcf, path, 'epsc');
    end
end

%% Box plot
%
% for aa = 1:nAlgorithms
%    tmp = allSAR(:,:,aa);
%    sarVec(:,aa) = tmp(:);
%
%    tmp = allSDR(:,:,aa);
%    sdrVec(:,aa) = tmp(:);
%
%    tmp = allSIR(:,:,aa);
%    sirVec(:,aa) = tmp(:);
% end
% betas = 1:nAlgorithms;
% figure
% subplot(3,1,1)
% errorbar(betas, nanmean(sarVec), nanstd(sarVec))
% set(gca,'xtick',betas,'xticklabel',algorithms);
%
% title('SAR')
% subplot(3,1,2)
% errorbar(betas, nanmean(sdrVec), nanstd(sdrVec))
% set(gca,'xtick',betas,'xticklabel',algorithms);
% title('SDR')
%
% subplot(3,1,3)
% errorbar(betas, nanmean(sirVec), nanstd(sirVec))
% set(gca,'xtick',betas,'xticklabel',algorithms);
% title('SIR')


%% Evaluate
winnerIdx = [2, 1;
    2, 3;
    2,4;
    2,5;
    2,6;
    3, 1;
    3,4;
    3,5;
    3,6;
    ];
allData = length(experimentList);
nWinner = size(winnerIdx,1);
for winIdx = 1:nWinner
    method1 = winnerIdx(winIdx,1);
    method2 = winnerIdx(winIdx,2);
    sarWin = sum(dataSAR{method1}(:) > dataSAR{method2}(:));
    sdrWin = sum(dataSDR{method1}(:) > dataSDR{method2}(:));
    sirWin = sum(dataSIR{method1}(:) > dataSIR{method2}(:));
    disp([algorithms{method1} ' vs ' algorithms{method2} ':'])
    disp(['SAR :' num2str(sarWin)...
        '/' num2str(allData) ' (' num2str(sarWin/allData*100,3) '%)']);
    disp(['SDR :' num2str(sdrWin)...
        '/' num2str(allData) ' (' num2str(sdrWin/allData*100,3) '%)']);
    disp(['SIR :' num2str(sirWin)...
        '/' num2str(allData) ' (' num2str(sirWin/allData*100,3) '%)']);
    fprintf('\n');
end

%%
figureDir = ['images' filesep metricsPath];
if ~exist(figureDir, 'dir')
    mkdir(figureDir);
end

idx = 1:9;
for aa = 1:nAlgorithms
    aux = allSAR(idx,idx,aa);
    aux = aux(~isnan(aux));
    alignSAR(:,aa) =  aux(:);

    aux = allSDR(idx,idx,aa);
    aux = aux(~isnan(aux));
    alignSDR(:,aa) =  aux(:);

    aux = allSIR(idx,idx,aa);
    aux = aux(~isnan(aux));
    alignSIR(:,aa) =  aux(:);
end

figure
set(gcf,'position',[10 150 2300 1000])
subplot(1,3,1)
boxplot(alignSAR, algorithms)
title('All sources SAR')
subplot(1,3,2)
boxplot(alignSDR, algorithms)
title('All sources SDR')
subplot(1,3,3)
boxplot(alignSIR, algorithms)
title('All sources SIR')
fileName = [figureDir filesep experiment 'all_source_metrics'];
set(gca,'LooseInset',get(gca,'TightInset'));

saveas(gcf, fileName, 'epsc');

%%
letterName = {'FastNMF', 'RS-MNMF_{\beta=0.9}', 'RS-MNMF_{\beta=0}', ...
    'DOA-MNMF', 'WN-MNMF', 'BS-MNMF' 'ILRMA'};
figureDir = ['images' filesep sourceSigType num2str(nSource) '/'];
if ~exist(figureDir, 'dir')
    mkdir(figureDir);
end
figSize = [10 150 700 150];
figure
set(gcf,'position',figSize)
errorbar(1:nAlgorithms, mean(alignSAR), std(alignSAR))
xlim([0,nAlgorithms+1]), ylim([-5,20]);
ylabel('[dB]'), grid on
set(gca,'xtick',1:nAlgorithms,'xticklabel',letterName);
title([num2str(nSource) ' ' sourceSigType ' sources'  ' SAR'])
fileName = [figureDir 'SAR_source_metrics'];
saveas(gcf, fileName, 'epsc');

figure
set(gcf,'position',figSize)
errorbar(1:nAlgorithms, mean(alignSDR), std(alignSDR))
xlim([0,nAlgorithms+1]);
set(gca,'xtick',1:nAlgorithms,'xticklabel',algorithms);
xlim([0,nAlgorithms+1]), ylim([-5,20]);
ylabel('[dB]'), grid on
set(gca,'xtick',1:nAlgorithms,'xticklabel',letterName);
title([num2str(nSource) ' ' sourceSigType ' sources'  ' SDR'])
fileName = [figureDir 'SDR_source_metrics'];
saveas(gcf, fileName, 'epsc');

figure
set(gcf,'position',figSize)
errorbar(1:nAlgorithms, mean(alignSIR), std(alignSIR))
xlim([0,nAlgorithms+1]), ylim([-5,20]);
ylabel('[dB]'), grid on
set(gca,'xtick',1:nAlgorithms,'xticklabel',letterName);
title([num2str(nSource) ' ' sourceSigType ' sources' ' SIR'])
fileName = [figureDir 'SIR_source_metrics'];
saveas(gcf, fileName, 'epsc');









