beta= 0:0.1:2;

aa = {'female','male','music'};

for bb = 1:length(aa)
    rsmcnmf_example(aa{bb},1.4)
end