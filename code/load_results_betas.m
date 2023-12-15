clear all; clc;

algorithms = {'beta0','beta1','beta2','beta3','beta4','beta5','beta6','beta7','beta8','beta9','beta10','beta11','beta12','beta13','beta14','beta15','beta16','beta17','beta18','beta19','beta20'};
signals = {'female','male','music'};
signals = {'female'};
sources = {'','3'};
sources = {''};

for cc = 1:length(sources)
    for ss = 1:length(signals)
        cont = 0;
        D = dir(['metrics_NEW_' signals{ss} sources{cc}]);
        for mm = 1:length(D)
            if D(mm).name(1) ~= '.' && D(mm).name(1) ~= 'i'
                cont =  cont + 1;
                res{cont,1} = D(mm).name;
                for aa = 1:length(algorithms)
                    load(['metrics_NEW_' signals{ss} sources{cc} filesep D(mm).name filesep algorithms{aa} '_metrics.mat']);
                    res{cont,(aa-1)*3 + 2} = mean(SDR(:),"omitnan");
                    res{cont,(aa-1)*3 + 3} = mean(SIR(:),"omitnan");
                    res{cont,(aa-1)*3 + 4} = mean(SAR(:),"omitnan");
                end
            end
        end
    end
    writecell(res,['betas_NEW.xls'])
end
