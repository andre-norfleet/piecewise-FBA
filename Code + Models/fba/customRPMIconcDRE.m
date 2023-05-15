function [RPMIconcvaluevector] = customRPMIconcDRE()
load('Y:\Riya\From Andre\FBA-pipeline-master-Riya\Code + Models\data\recon\recon3d_qflux.mat');

RPMIchoice = fopen('../data/media/output/RPMI1640.csv','r');
RPMIvardata1 = textscan(RPMIchoice,'%s %s %f %f','Delimiter',',','headerLines',1);
fclose(RPMIchoice);

RPMIconcvaluevector = [];
RPMIexchrxns = [];

for i = 1:length(model.rxns)
        for kk = 1:length(RPMIvardata1{1,2})
            if any(strcmp(RPMIvardata1{1,2}{kk},model.rxns{i}))
                RPMIexchrxns(end+1) = i; %index in model.rxns for all exchange rxns accompanying mTeSR media component
            end
        end
end

for nn = 1:length(RPMIvardata1{1,3})
    for nnn = 1:length(RPMIexchrxns)
        RPMIconcvaluevector(1,nnn) =  RPMIvardata1{1,3}(nnn);
    end
end