%function [mTeSRconcvaluevector,mTeSRexchrxns] = customextracellcytosoltransportDRE()
load('Y:\Riya\From Andre\FBA-pipeline-master-Riya\Code + Models\data\recon\recon3d_qflux.mat');
today = model.subSystems;
today2 = model.rxns;
today3 = model.mets;
today4 = model.ub;
today5 = model.S;
allglucmetabs = find(contains(today3,'glc_D['));

mTeSRchoice = fopen('../data/media/output/mTeSR.csv','r');
mediavardata1 = textscan(mTeSRchoice,'%s %s %f %f','Delimiter',',','headerLines',1);
fclose(mTeSRchoice);

mTeSR1stmetab = mediavardata1{1,1}{1,1}(1:end-3);
mTeSRmetabonlyvector = mediavardata1{1,1}(1:end-3);
glucose = mediavardata1{1,1}{12,1}(1:end-3);

mTeSRmetabvector = {};
allmTeSRmetabindices = [];
for n = 1:length(mediavardata1{1,1})
    mTeSRmetabvector{n,end+1} = find(contains(model.mets,mediavardata1{1,1}{n,1}(1:end-3) ));
end

syntcheck = mTeSRmetabvector{12,12};
syntcheck2 = mTeSRmetabvector{12,12}';

%find glc_D
glc_d = []; %saves indices of all glucose metabolites in m direction
for k = 1:length(model.mets) %parse through all metabolites
    if length(model.mets{k}) == 8 %metabolite is 8 characters in length (ex: glc_D[e])
        if strcmp(model.mets{k}(1:6),'glc_D[') %metabolite matches gluc 
            glc_d(end+1) = k; %add index 
        end
    end
end

mTeSRmetabos2 = [];
for k = 1:length(today3) %parse through all metabolites
    for kk = 1:length(mediavardata1{1,1})
            if strcmp(today3{k},mediavardata1{1,1}{kk,1}) %metabolite matches gluc (all compartment forms essentially)
                mTeSRmetabos2(end+1) = k; %index for extracell form of each mTeSR metab
            end
    end
end

m2babyy = [];

for k = 1:length(today2)
    m2babyy = find(full(today5(:,k)) < 0)';
end


m3babyy = [];

for k = 1:length(today2)
    if today4(k) > 0
        reactant_found2 = false;
        for m3babyy = find(full(today5(:,k)) < 0)'
            if any(m3babyy == mTeSRmetabos2)
                reactant_found2 = true;
                break;
            end
        end
    end
end

reactant_found2
