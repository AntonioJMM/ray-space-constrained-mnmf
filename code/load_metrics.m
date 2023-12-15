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
homeFolder = './';

sdr_total = [];
sir_total = [];
sar_total = [];

% sourceS = {'female','male','female3','male3'};
sourceS = {'female','male','music'};
sourceS = {'music'};
for zz = 1:length(sourceS)
    sourceSigType = sourceS{zz};
    resultsPath = [homeFolder 'results_NEW_sinREST' sourceSigType];   % Path with estimates of the BSS
    metricsPath = [homeFolder 'metrics_NEW_sinREST' sourceSigType];%'metrics_music3';

    idxExp = 1;
    if idxExp == 1
        experiment = 'rev_';
    else
        experiment = 'test_';
    end


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
    experimentList = {'rev_30_30_12','rev_30_30_31','rev_30_60_13','rev_60_60_13','rev_60_90_13','rev_60_90_31','rev_90_60_21'};
    for tt = 1:length(experimentList)
        % Load estimed signals
        % Directory of the experiment
        experimentName = experimentList{tt};
        metricsPathExperimentPath = [metricsPath, filesep, experimentName];

        load([metricsPathExperimentPath filesep 'beta14_metrics.mat'])

        sdr_total = [sdr_total;SDR(:)];
        sir_total = [sir_total;SIR(:)];
        sar_total = [sar_total;SAR(:)];

    end
end

median(sdr_total)
