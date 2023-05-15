load('Z:/Andre/From Andre/FBA-pipeline-master-Riya/Code + Models/data/recon/recon3d_qflux.mat');
%results1 = tdfread('Z:/Andre/From Andre/FBA-pipeline-master-Riya/Code + Models/fba/results/WTC11_D0_Zhangetal_ATPallmTeSRplusAlbuMAX_relaxconstralltest_4/_TEMPLATE_1/fva/WTC_-dox_avg.tsv','\t');
%results2 = tdfread('Z:/Andre/From Andre/FBA-pipeline-master-Riya/Code + Models/fba/results/WTC11_D0_Zhangetal_ATPallmTeSRplusAlbuMAX_allexchconstralltest_3/_TEMPLATE_1/fva/WTC_-dox_avg.tsv','\t');
%resultsdata = readmatrix('Z:/Andre/New FBA Results/WTC11_D0_mTeSRvsmTeSRplusAlbuMAX_ALLexchconstraintest2FVA1.csv');
%resultsdata = readmatrix('Z:/Andre/New FBA Results/WTC11D0VsD2FVA_ATPall_ctrl.csv');
resultsdata = readmatrix('Z:/Andre/New FBA Results/WTC11D0_Turneretal_NormVsmTeSRpluspalmitateConstrainFVA.csv');
resultsrxns = tdfread('Z:/Andre/From Andre/FBA-pipeline-master-Riya/Code + Models/fba/results/WTC11_D0_Zhangetal_ATPallmTeSRplusAlbuMAX_relaxconstralltest_4/_TEMPLATE_1/fva/WTC_-dox_avg.tsv','\t');
orderrxnlist = resultsrxns.REACTION;

maxvals1 = resultsdata(:,2);
maxvals2 = resultsdata(:,3);

%realrxnlist = model.rxns;
peeked = orderrxnlist(3,:);

%Ratios of treatment max vals: ctrl treatment max vals
ratiovals = maxvals2./maxvals1;

%Log2-fold Ratios of treatment max vals: ctrl treatment max vals
%Negative vals mean increase in flux for ctrl condition, dec for treat and vice versa 
log2ratiovals = log2(ratiovals);
log2ratiovalsreal = real(log2ratiovals);

wholelist = [maxvals1 maxvals2];

%dblorderrxnlist = convertCharsToStrings(orderrxnlist);

%Subsys-spec rxn indices and ID extraction
PPPrxnIndex = find(contains(model.subSystems,'Fatty acid synthesis'));
subsysrxnIDvec = {};
for i = 1:length(PPPrxnIndex)
    subsysrxnIDvec{end+1} = model.rxns(PPPrxnIndex(i,1)); %returns 1 x i cell
end

PPPmaxvals1 = [];
PPPmaxvals2 = [];
subsysLog2FCvals = [];
subsysLog2FCvalsreal = [];

for ii = 1:length(model.rxns)
    for jj = 1:length(subsysrxnIDvec)
        PPPmaxvals1(jj,1) = maxvals1(PPPrxnIndex(jj,1));
        PPPmaxvals2(jj,1) = maxvals2(PPPrxnIndex(jj,1));
        subsysLog2FCvals(jj,1) = log2ratiovals(PPPrxnIndex(jj,1));
        subsysLog2FCvalsreal(jj,1) = log2ratiovalsreal(PPPrxnIndex(jj,1));
    end
end

charIDvec = [];

for ii = 1:length(model.rxns)
    for jj = 1:length(subsysrxnIDvec)
        charIDvec(jj,:) = orderrxnlist(jj,:);
    end
end

ctrlsubsysavg = mean(PPPmaxvals1);
treatsubsysavg = mean(PPPmaxvals2);
ctrlsubsystddev = std(PPPmaxvals1);
treatsubsystddev = std(PPPmaxvals2);

fluxdata = [PPPmaxvals1 PPPmaxvals2];
[h,p] = ttest(PPPmaxvals1,PPPmaxvals2,'Alpha',0.05)

% %file save encoding scheme
% f_fva = fopen('Z:/Andre/FVA subsystem results analysis/mtesrvsmtesrnoglucglutOxPhos_CSVstats_FVAsubsysautomateresults.tsv','w');
% fprintf(f_fva,'REACTION\tMAX ctrl [mmol/gDW/hr]\tMAX treat[mmol/gDW/hr]\tctrl avg\tctrl stddev\ttreat avg\ttreat stddev\n');
% %Print your transformed data to excel file
% for j = 1:length(subsysrxnIDvec)
%     fprintf(f_fva,'%s\t%0.9f\t%0.9f\t%0.9f\t%0.9f\t%0.9f\t%0.9f\n',model.rxns{PPPrxnIndex(j)},PPPmaxvals1(j),PPPmaxvals2(j),ctrlsubsysavg,ctrlsubsystddev,treatsubsysavg,treatsubsystddev);
% end
%fprintf(f_fva,'t%0.9f\n',ctrlsubsysavg)
%fclose(f_fva);
yvalues = subsysrxnIDvec';
xvalues = {'Log2FC D0 vs D4'};
%h = heatmap(xvalues,yvalues,fluxdata);
%h = heatmap(xvalues,yvalues,subsysLog2FCvals);
% Use below if complex numbers in subsysLog2FCvals vector
h = heatmap(xvalues,yvalues,subsysLog2FCvalsreal,'Colormap', hot);