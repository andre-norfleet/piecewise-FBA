%readin = readmatrix('Z:/Andre/From Andre/FBA-pipeline-master-Riya/Code + Models/fba/results/WTC11_D0_Zhangetal_ATPallmTeSRplusAlbuMAX_relaxconstralltest_4/_TEMPLATE_1/fva/WTC_-dox_avg.tsv');
%data = dlmread('Z:/Andre/From Andre/FBA-pipeline-master-Riya/Code + Models/fba/results/WTC11_D0_Zhangetal_ATPallmTeSRplusAlbuMAX_relaxconstralltest_4/_TEMPLATE_1/fva/WTC_-dox_avg.tsv','\t', 1, 1);
results1 = tdfread('Z:/Andre/From Andre/FBA-pipeline-master-Riya/Code + Models/fba/results/WTC11_D0_Zhangetal_ATPallmTeSRplusAlbuMAX_relaxconstralltest_4/_TEMPLATE_1/fva/WTC_-dox_avg.tsv','\t');
results2 = tdfread('Z:/Andre/From Andre/FBA-pipeline-master-Riya/Code + Models/fba/results/WTC11_D0_Zhangetal_ATPallmTeSRplusAlbuMAX_allexchconstralltest_3/_TEMPLATE_1/fva/WTC_-dox_avg.tsv','\t');

minvals1 = results1.MIN_0x5Bmmol0x2FgDW0x2Fhr0x5D;
maxvals1 = results1.MAX_0x5Bmmol0x2FgDW0x2Fhr0x5D;
minvals2 = results2.MIN_0x5Bmmol0x2FgDW0x2Fhr0x5D;
maxvals2 = results2.MAX_0x5Bmmol0x2FgDW0x2Fhr0x5D;


orderrxnlist = results1.REACTION;
peeked = orderrxnlist(3,:);

%Ratios of treatment max vals: ctrl treatment max vals
ratiovals = maxvals2./maxvals1;

%Log2-fold Ratios of treatment max vals: ctrl treatment max vals
log2ratiovals = log2(ratiovals);

wholelist = [maxvals1 maxvals2];

dblorderrxnlist = convertCharsToStrings(orderrxnlist);

%file save encoding scheme
f_fva = fopen('Z:/Andre/testFVAautomateresults.tsv','w');
fprintf(f_fva,'REACTION\tMAX ctrl [mmol/gDW/hr]\tMAX treat[mmol/gDW/hr]\tLog2 ratiometric vals\n');

%Print your transformed data to excel file
for j = 1:length(orderrxnlist)
    fprintf(f_fva,'%s\t%0.9f\t%0.9f\t%0.9f\n',orderrxnlist(j),maxvals1(j),maxvals2(j),log2ratiovals(j));
end
