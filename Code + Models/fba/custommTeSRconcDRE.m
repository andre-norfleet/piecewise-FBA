function [mTeSRconcvaluevector] = custommTeSRconcDRE()
load('Z:\Riya\From Andre\FBA-pipeline-master-Riya\Code + Models\data\recon\recon3d_qflux.mat');

mTeSRchoice = fopen('../data/media/output/mTeSR.csv','r');
mediavardata1 = textscan(mTeSRchoice,'%s %s %f %f','Delimiter',',','headerLines',1);
fclose(mTeSRchoice);

mTeSRconcvaluevector = [];
mTeSRexchrxns = [];

for i = 1:length(model.rxns)
        for kk = 1:length(mediavardata1{1,2})
            if any(strcmp(mediavardata1{1,2}{kk},model.rxns{i}))
                mTeSRexchrxns(end+1) = i; %index in model.rxns for all exchange rxns accompanying mTeSR media component
            end
        end
end

for nn = 1:length(mediavardata1{1,3})
    for nnn = 1:length(mTeSRexchrxns)
        mTeSRconcvaluevector(1,nnn) =  mediavardata1{1,3}(nnn);
    end
end