function [GFAMexchrxns] = customexchangetransportGFAMDRE()
load('Y:\Riya\From Andre\FBA-pipeline-master-Riya\Code + Models\data\recon\recon3d_qflux.mat');
today = model.subSystems;
today2 = model.rxns;
today3 = model.mets;
today4 = model.ub;
today5 = model.S;
%allglucmetabs = find(contains(today3,'glc_D['));

GFAMchoice = fopen('../data/media/output/RPMI1640GFAM.csv','r');
GFAMvardata1 = textscan(GFAMchoice,'%s %s %f %f','Delimiter',',','headerLines',1);
fclose(GFAMchoice);

GFAMexchrxns = [];

for i = 1:length(model.rxns)
        for kk = 1:length(GFAMvardata1{1,2})
            if any(strcmp(GFAMvardata1{1,2}{kk},model.rxns{i}))
                GFAMexchrxns(end+1) = i; %index in model.rxns for all exchange rxns accompanying mTeSR media component
            end
        end
end

GFAMexchrxns