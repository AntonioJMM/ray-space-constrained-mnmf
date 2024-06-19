clc; clear all; close all;
lambda = 2.^(2:10);
% tik    = 1:0.5:12;
tik = perms(6:10);
tik = tik(:,1:3);
tik = unique(tik,'rows');
tik = sortrows(reshape(permute(repmat(tik,1,1,3),[1 3 2]),length(tik)*3,3));
tik = [2 2 2];

for ii = 1:length(lambda)*length(tik)
    jj = mod(ii,length(lambda)); if jj == 0, jj=length(lambda); end
    mm = ceil(ii/length(lambda));

    %     fprintf('%d - %d\n',jj,mm)

    [raySDR{ii},raySIR{ii},raySAR{ii}] = rsmcnmf_example(lambda(jj),tik(mm,:));
    %     [raySDR{ii},raySIR{ii},raySAR{ii}] = rsmcnmf_example(lambda,tik(ii,:));
end