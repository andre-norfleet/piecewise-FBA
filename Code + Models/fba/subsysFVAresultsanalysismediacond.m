load('Z:/Andre/From Andre/FBA-pipeline-master-Riya/Code + Models/data/recon/recon3d_qflux.mat');
results1 = tdfread('Z:/Andre/From Andre/FBA-pipeline-master-Riya/Code + Models/fba/results/WTC11_D0_Zhangetal_ATPallmTeSRplusAlbuMAX_relaxconstralltest_4/_TEMPLATE_1/fva/WTC_-dox_avg.tsv','\t');
results2 = tdfread('Z:/Andre/From Andre/FBA-pipeline-master-Riya/Code + Models/fba/results/WTC11_D0_Zhangetal_ATPallmTeSRplusAlbuMAX_allexchconstralltest_3/_TEMPLATE_1/fva/WTC_-dox_avg.tsv','\t');

minvals1 = results1.MIN_0x5Bmmol0x2FgDW0x2Fhr0x5D;
maxvals1 = results1.MAX_0x5Bmmol0x2FgDW0x2Fhr0x5D;
minvals2 = results2.MIN_0x5Bmmol0x2FgDW0x2Fhr0x5D;
maxvals2 = results2.MAX_0x5Bmmol0x2FgDW0x2Fhr0x5D;


orderrxnlist = results1.REACTION;
realrxnlist = model.rxns;
peeked = orderrxnlist(3,:);

%Ratios of treatment max vals: ctrl treatment max vals
ratiovals = maxvals2./maxvals1;

%Log2-fold Ratios of treatment max vals: ctrl treatment max vals
log2ratiovals = log2(ratiovals);

wholelist = [maxvals1 maxvals2];

dblorderrxnlist = convertCharsToStrings(orderrxnlist);

%Subsys-spec rxn indices and ID extraction
oxphosrxnIndex = find(contains(model.subSystems,'Fatty acid oxidation'));
subsysrxnIDvec = {};
for i = 1:length(oxphosrxnIndex)
    subsysrxnIDvec{end+1} = model.rxns(oxphosrxnIndex(i,1)); %returns 1 x i cell
end

oxphosmaxvals1 = [];
oxphosmaxvals2 = [];

for ii = 1:length(model.rxns)
    for jj = 1:length(subsysrxnIDvec)
        oxphosmaxvals1(jj,1) = maxvals1(oxphosrxnIndex(jj,1));
        oxphosmaxvals2(jj,1) = maxvals2(oxphosrxnIndex(jj,1));
    end
end

charIDvec = [];

for ii = 1:length(model.rxns)
    for jj = 1:length(subsysrxnIDvec)
        charIDvec(jj,:) = orderrxnlist(jj,:);
    end
end

ctrlsubsysavg = mean(oxphosmaxvals1);
treatsubsysavg = mean(oxphosmaxvals2);
ctrlsubsystddev = std(oxphosmaxvals1);
treatsubsystddev = std(oxphosmaxvals2);

%file save encoding scheme
f_fva = fopen('Z:/Andre/test2FVAsubsysautomateresults.tsv','w');
fprintf(f_fva,'REACTION\tMAX ctrl [mmol/gDW/hr]\tMAX treat[mmol/gDW/hr]\n');
%Print your transformed data to excel file
for j = 1:length(subsysrxnIDvec)
    fprintf(f_fva,'%s\t%0.9f\t%0.9f\n',model.rxns{oxphosrxnIndex(j)},oxphosmaxvals1(j),oxphosmaxvals2(j));
end
fclose(f_fva);
