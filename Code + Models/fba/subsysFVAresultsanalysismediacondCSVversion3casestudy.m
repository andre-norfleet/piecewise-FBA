load('Z:/Andre/From Andre/FBA-pipeline-master-Riya/Code + Models/data/recon/recon3d_qflux.mat');
%results1 = tdfread('Z:/Andre/From Andre/FBA-pipeline-master-Riya/Code + Models/fba/results/WTC11_D0_Zhangetal_ATPallmTeSRplusAlbuMAX_relaxconstralltest_4/_TEMPLATE_1/fva/WTC_-dox_avg.tsv','\t');
%results2 = tdfread('Z:/Andre/From Andre/FBA-pipeline-master-Riya/Code + Models/fba/results/WTC11_D0_Zhangetal_ATPallmTeSRplusAlbuMAX_allexchconstralltest_3/_TEMPLATE_1/fva/WTC_-dox_avg.tsv','\t');
%resultsdata = readmatrix('Z:/Andre/New FBA Results/WTC11_D0_mTeSRvsmTeSRplusAlbuMAX_ALLexchconstraintest2FVA1.csv');
resultsdataA = readmatrix('Z:/Andre/New FBA Results/WTC11mTeSRvsmTeSRnogluc_allexchconstraintest2FVA1.csv');
%resultsdataB = readmatrix('Z:/Andre/New FBA Results/WTC11D0_Turneretal_NormVsmTeSRpluspalmitateConstrainFVA.csv');
resultsdataB = readmatrix('Z:/Andre/New FBA Results/WTC11mTeSRvsmTeSRminusglut_allexchconstraintest3.csv');
resultsdataC = readmatrix('Z:/Andre/New FBA Results/WTC11mTeSRvsmTeSRnoglucglut_allexchconstraintest2.csv');
resultsrxns = tdfread('Z:/Andre/From Andre/FBA-pipeline-master-Riya/Code + Models/fba/results/WTC11_D0_Zhangetal_ATPallmTeSRplusAlbuMAX_relaxconstralltest_4/_TEMPLATE_1/fva/WTC_-dox_avg.tsv','\t');
orderrxnlist = resultsrxns.REACTION;

%Study A max vals
maxvals1A = resultsdataA(:,2);
maxvals2A = resultsdataA(:,3);
%Study B max vals
maxvals1B = resultsdataB(:,2);
maxvals2B = resultsdataB(:,3);
%Study C max vals
maxvals1C = resultsdataB(:,2);
maxvals2C = resultsdataB(:,3);

%realrxnlist = model.rxns;
peeked = orderrxnlist(3,:);

%Ratios of treatment max vals: ctrl treatment max vals
ratiovalsA = maxvals2A./maxvals1A;
ratiovalsB = maxvals2B./maxvals1B;
ratiovalsC = maxvals2C./maxvals1C;

%Log2-fold Ratios of treatment max vals: ctrl treatment max vals
%Negative vals mean increase in flux for ctrl condition, dec for treat and vice versa 
log2ratiovalsA = log2(ratiovalsA);
log2ratiovalsrealA = real(log2ratiovalsA);
log2ratiovalsB = log2(ratiovalsB);
log2ratiovalsrealB = real(log2ratiovalsB);
log2ratiovalsC = log2(ratiovalsC);
log2ratiovalsrealC = real(log2ratiovalsC);

%wholelist = [maxvals1 maxvals2];

%dblorderrxnlist = convertCharsToStrings(orderrxnlist);

%Subsys-spec rxn indices and ID extraction
PPPrxnIndex = find(contains(model.subSystems,'Citric acid cycle'));
subsysrxnIDvec = {};
for i = 1:length(PPPrxnIndex)
    subsysrxnIDvec{end+1} = model.rxns(PPPrxnIndex(i,1)); %returns 1 x i cell
end

PPPmaxvals1A = [];
PPPmaxvals2A = [];
PPPmaxvals1B = [];
PPPmaxvals2B = [];
PPPmaxvals1C = [];
PPPmaxvals2C = [];
subsysLog2FCvalsA = [];
subsysLog2FCvalsrealA = [];
subsysLog2FCvalsB = [];
subsysLog2FCvalsrealB = [];
subsysLog2FCvalsC = [];
subsysLog2FCvalsrealC = [];

for ii = 1:length(model.rxns)
    for jj = 1:length(subsysrxnIDvec)
        PPPmaxvals1A(jj,1) = maxvals1A(PPPrxnIndex(jj,1));
        PPPmaxvals2A(jj,1) = maxvals2A(PPPrxnIndex(jj,1));
        PPPmaxvals1B(jj,1) = maxvals1B(PPPrxnIndex(jj,1));
        PPPmaxvals2B(jj,1) = maxvals2B(PPPrxnIndex(jj,1));
        PPPmaxvals1C(jj,1) = maxvals1C(PPPrxnIndex(jj,1));
        PPPmaxvals2C(jj,1) = maxvals2C(PPPrxnIndex(jj,1));
        subsysLog2FCvalsA(jj,1) = log2ratiovalsA(PPPrxnIndex(jj,1));
        subsysLog2FCvalsrealA(jj,1) = log2ratiovalsrealA(PPPrxnIndex(jj,1));
        subsysLog2FCvalsB(jj,1) = log2ratiovalsB(PPPrxnIndex(jj,1));
        subsysLog2FCvalsrealB(jj,1) = log2ratiovalsrealB(PPPrxnIndex(jj,1));
        subsysLog2FCvalsC(jj,1) = log2ratiovalsC(PPPrxnIndex(jj,1));
        subsysLog2FCvalsrealC(jj,1) = log2ratiovalsrealC(PPPrxnIndex(jj,1));
    end
end

subsysLog2FCvalsreal3 =[subsysLog2FCvalsrealA subsysLog2FCvalsrealB subsysLog2FCvalsrealC];

charIDvec = [ ];

for ii = 1:length(model.rxns)
    for jj = 1:length(subsysrxnIDvec)
        charIDvec(jj,:) = orderrxnlist(jj,:);
    end
end

% ctrlsubsysavg = mean(PPPmaxvals1);
% treatsubsysavg = mean(PPPmaxvals2);
% ctrlsubsystddev = std(PPPmaxvals1);
% treatsubsystddev = std(PPPmaxvals2);
% 
% fluxdata = [PPPmaxvals1 PPPmaxvals2];
% [h,p] = ttest(PPPmaxvals1,PPPmaxvals2,'Alpha',0.05)

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
xvalues = {'Log2FC Tohyama ctrl vs (-)glc','Log2FC Tohyama ctrl vs (-)gln', 'Log2FC Tohyama ctrl vs (-)glc/gln'};
%h = heatmap(xvalues,yvalues,fluxdata);
%h = heatmap(xvalues,yvalues,subsysLog2FCvals);
% Use below if complex numbers in subsysLog2FCvals vector
%h = heatmap(xvalues,yvalues,subsysLog2FCvalsrealA,subsysLog2FCvalsrealB);
h = heatmap(xvalues,yvalues,subsysLog2FCvalsreal3,'Colormap', hot);