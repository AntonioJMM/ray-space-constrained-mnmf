clear all; clc;

algorithm = {'beamSpace','beta0','beta9','CNMF','fast','ilrma','WNNTF'};
signals = {'female','male','music'};
sources = {'','3'};

path_mirco = '../../../mpezzoli/NMF_RS_separation/review/';


for ss = 1:length(signals)
    for aa = 1:length(algorithm)
        cont = 0;
        for cc = 1:length(sources)
            D = dir([path_mirco filesep 'metrics_' signals{ss} sources{cc} filesep 'rev*']);
            for mm = 1:length(D)
                if D(mm).name(1) ~= '.' && D(mm).name(1) ~= 'i'
                    load([path_mirco filesep 'metrics_' signals{ss} sources{cc} filesep D(mm).name filesep algorithm{aa} '_metrics.mat']);
                    cont =  cont + 1;

                    resSDR(cont,aa,ss) = mean(SDR(:),"omitnan");
                    resSIR(cont,aa,ss) = mean(SIR(:),"omitnan");
                    resSAR(cont,aa,ss) = mean(SAR(:),"omitnan");
                end
            end
        end
    end
end

for ss = 1:length(signals)
    cont = 0;
    for cc = 1:length(sources)
        D = dir(['metrics_NEW_' signals{ss} sources{cc}]);
        for mm = 1:length(D)
            if D(mm).name(1) ~= '.' && D(mm).name(1) ~= 'i'
                load(['metrics_NEW_' signals{ss} sources{cc} filesep D(mm).name filesep 'beta14' '_metrics.mat']);
                cont =  cont + 1;
                resSDR(cont,8,ss) = mean(SDR(:),"omitnan");
                resSIR(cont,8,ss) = mean(SIR(:),"omitnan");
                resSAR(cont,8,ss) = mean(SAR(:),"omitnan");
            end
        end
    end
end
