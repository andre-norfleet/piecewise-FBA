function [GFAMconcvaluevector] = customGFAMconcDRE()
load('Z:\Riya\From Andre\FBA-pipeline-master-Riya\Code + Models\data\recon\recon3d_qflux.mat');

GFAMchoice = fopen('../data/media/output/RPMI1640GFAM.csv','r');
GFAMvardata1 = textscan(GFAMchoice,'%s %s %f %f','Delimiter',',','headerLines',1);
fclose(GFAMchoice);

GFAMconcvaluevector = [];
GFAMexchrxns = [];

for i = 1:length(model.rxns)
        for kk = 1:length(GFAMvardata1{1,2})
            if any(strcmp(GFAMvardata1{1,2}{kk},model.rxns{i}))
                GFAMexchrxns(end+1) = i; %index in model.rxns for all exchange rxns accompanying mTeSR media component
            end
        end
end

for nn = 1:length(GFAMvardata1{1,3})
    for nnn = 1:length(GFAMexchrxns)
        GFAMconcvaluevector(1,nnn) =  GFAMvardata1{1,3}(nnn);
    end
end