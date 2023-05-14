function [mTeSRconcvaluevector,mTeSRexchrxns] = customextracellcytosoltransportDRE()
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

%find glc_D
glc_d = []; %saves indices of all glucose metabolites in m direction
mTeSRexchrxns = [];
for k = 1:length(today3) %parse through all metabolites
    if length(today3{k}) == 8 %metabolite is 8 characters in length (ex: glc_D[e])
        if strcmp(today3{k}(1:6),'glc_D[') %metabolite matches gluc (all compartment forms essentially)
            glc_d(end+1) = k; %add index 
        end
    end
end

metabo1 = mediavardata1{1,2}{1,1};
metabo2 = mediavardata1{1,2}{2,1};

jean = model.rxns{13};

for i = 1:length(model.rxns)
    %for k = 1:length(model.mets)
        for kk = 1:length(mediavardata1{1,2})
            if any(strcmp(mediavardata1{1,2}{kk},model.rxns{i}))
                mTeSRexchrxns(end+1) = i; %index in model.rxns for all exchange rxns accompanying mTeSR media component
            end
        end
    %end
end

mTeSRexchrxns

mTeSRconcvaluevector = [];

for nn = 1:length(mediavardata1{1,3})
    for nnn = 1:length(mTeSRexchrxns)
        mTeSRconcvaluevector(1,nnn) =  mediavardata1{1,3}(nnn);
    end
end

mTeSRconcvaluevector
                   
 % find reactions with glc_D
 rxnlist = [];
 rightones2 = {};


 for k = 1:length(today2) %parse through rxns 
     % check forward direction
     if today4(k) > 0 %check if upper bound flux is positive -- not sure why we need to check this
         % find glc_d reactant
         reactant_found = false; %set glc_d as a reactant as false
         for m = find(full(today5(:,k)) < 0)' %m is vector with indices of rxns with glc_d as reactant being consumed
             if any(m == glc_d) %returns 1 for all rxns with glc_d as reactant 
                 reactant_found = true; %set glc_d as a reactant as true
                 break;
             end
         end
 
         % find product
         product_found = false;
         for m = find(full(today5(:,k)) > 0)'
             if any(m == glc_d)
                 product_found = true;
                 break;
             end
         end
         % add
         if reactant_found && product_found && ~any(k == rxnlist)
             rxnlist(end+1) = k;
         end
     end
     
     % check reverse direction
     if today4(k) < 0
         % find reactant
         reactant_found = false;
         for m = find(full(today5(:,k)) > 0)'
             if any(m == glc_d)
                 reactant_found = true;
                 break;
             end
         end
         % find product
         product_found = false;
         for m = find(full(today5(:,k)) < 0)'
             if any(m == glc_d)
                 product_found = true;
                 break;
             end
         end
         % add
         if reactant_found && product_found && ~any(k == rxnlist)
             rxnlist(end+1) = k;
         end
     end
 end
 rxnlist
 for ii = 1:length(rxnlist)
    rightones2{ii} = model.rxns{rxnlist(ii),1};
 end