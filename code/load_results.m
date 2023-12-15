clear all; clc;

signals = {'female','male','music'};
% signals = {'female'};
% sources = {'','3'};
sources = {'3'};

for cc = 1:length(sources)
    for ss = 1:length(signals)
        cont = 0;
        D = dir(['metrics_NEW_' signals{ss} sources{cc}]);
        for mm = 1:length(D)
            if D(mm).name(1) ~= '.' && D(mm).name(1) ~= 'i'
                load(['metrics_NEW_' signals{ss} sources{cc} filesep D(mm).name filesep 'beta14' '_metrics.mat']);
                cont =  cont + 1;
                res{cont,1} = D(mm).name;
                res{cont,2} = mean(SDR(:),"omitnan");
                res{cont,3} = mean(SIR(:),"omitnan");
                res{cont,4} = mean(SAR(:),"omitnan");
            end
        end
        writecell(res,[signals{ss} sources{cc} '_NEW.xls'])
    end
end
