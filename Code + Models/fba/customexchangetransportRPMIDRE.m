function [RPMIexchrxns] = customexchangetransportRPMIDRE()
load('Z:\Riya\From Andre\FBA-pipeline-master-Riya\Code + Models\data\recon\recon3d_qflux.mat');
today = model.subSystems;
today2 = model.rxns;
today3 = model.mets;
today4 = model.ub;
today5 = model.S;
%allglucmetabs = find(contains(today3,'glc_D['));

RPMIchoice = fopen('../data/media/output/RPMI1640.csv','r');
RPMIvardata1 = textscan(RPMIchoice,'%s %s %f %f','Delimiter',',','headerLines',1);
fclose(RPMIchoice);

RPMIexchrxns = [];

for i = 1:length(model.rxns)
        for kk = 1:length(RPMIvardata1{1,2})
            if any(strcmp(RPMIvardata1{1,2}{kk},model.rxns{i}))
                RPMIexchrxns(end+1) = i; %index in model.rxns for all exchange rxns accompanying mTeSR media component
            end
        end
end

RPMIexchrxns