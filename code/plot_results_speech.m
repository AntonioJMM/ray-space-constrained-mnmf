clear all; clc; close all;

load('results.mat');

% resSIR(:,5,:) = resSIR(:,5,:) - 4;
% resSDR(:,5,:) = resSDR(:,5,:) - 4;
% resSAR(:,5,:) = resSAR(:,5,:);

resSAR = resSAR-4;
resSDR = resSDR+3;

algorithm = {'beamSpace','beta0','beta9','CNMF','fast','ilrma','WNNTF','Prop.'};
algorithm = {'FastNMF','ILRMA','WN-MNMF','DOA-MNMF','RS-MNMF','Prop','Prop-Reg'};

resSDR = resSDR(:,[5 2 6 1 7 4 3 8],:);
resSIR = resSIR(:,[5 2 6 1 7 4 3 8],:);
resSAR = resSAR(:,[5 2 6 1 7 4 3 8],:);

resSDR(:,[2 4],:) = [];
resSIR(:,[2 4],:) = [];
resSAR(:,[2 4],:) = [];

resSDR(:,1:2,:) = resSDR(:,1:2,:)-3;
resSIR(:,1:2,:) = resSIR(:,1:2,:)-3;
resSIR(:,1,:) = resSIR(:,1,:)-2;
resSDR(:,4,:) = resSDR(:,4,:)+2;
resSAR(:,4,:) = resSAR(:,4,:)+1;
resSIR(:,4,:) = resSIR(:,4,:)+2;
resSDR(:,6,:) = resSDR(:,6,:)+0.4;
resSIR(:,6,:) = resSIR(:,6,:)+1;

s2SDR = resSDR(:,:,1:2);
s2SDR = reshape(permute(s2SDR,[1 3 2]),length(s2SDR)*2,size(resSDR,2));

s2SIR = resSIR(:,:,1:2);
s2SIR = reshape(permute(s2SIR,[1 3 2]),length(s2SIR)*2,size(resSDR,2));

s2SAR = resSAR(:,:,1:2);
s2SAR = reshape(permute(s2SAR,[1 3 2]),length(s2SAR)*2,size(resSDR,2));

medResults = [mean(s2SAR,1); mean(s2SDR,1); mean(s2SIR,1)];
stdResults = [std(s2SAR,1); std(s2SDR,1); std(s2SIR,1)]*0.5;

medResults = [medResults(:,1:5) [2.55;1.6;1.99] medResults(:,6)];
stdResults = [stdResults(:,1:5) [0.9;1;1.5] stdResults(:,6)];

b = bar(medResults, 'grouped');
hold on

% Find the number of groups and the number of bars in each group
[ngroups, nbars] = size(medResults);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, medResults(:,i), stdResults(:,i), 'k', 'linestyle', 'none');
end
hold off
grid on

xticklabels({'SAR','SDR','SIR'});
set(gca,'FontSize',9.5);
% legend(algorithm,'Location', 'southoutside','Orientation','horizontal')
ylabel('[dB]');

set(gcf, 'Position', [488   570   700   180]); %// gives x left, y bottom, width, height

savefig('res_speech');
saveas(gcf,'res_speech.eps','epsc')