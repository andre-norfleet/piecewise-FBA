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

mTeSRmetabvector = {};
allmTeSRmetabindices = [];
for n = 1:length(mediavardata1{1,1})
    mTeSRmetabvector{n,end+1} = find(contains(model.mets,mediavardata1{1,1}{n,1}(1:end-3) ));
end

mTeSRmetabvector2 = {};

for n = 1:length(mediavardata1{1,1})
    mTeSRmetabvector2{n,end+1} = mTeSRmetabvector{n,n}'; %transposing for downstream index comparison to reactant_found and product_found
end

%hoping = mTeSRmetabvector2{12,12}

%for nn = 1:length(mediavardata1{1,1})
%    for nnn = 1:length(mTeSRexchrxns)
%        mTeSRmetabvector(1,nnn) =  mediavardata1{1,1}(nnn);
%    end
%end

%find glc_D
glc_d = []; %saves indices of all glucose metabolites in m direction
mTeSRmetabos2 = [];
for k = 1:length(today3) %parse through all metabolites
    for kk = 1:length(mediavardata1{1,1})
        %if length(today3{k}) == 8 %metabolite is 8 characters in length (ex: glc_D[e])
            if any(strcmp(today3{k},mediavardata1{1,1}{kk})) %metabolite matches gluc (all compartment forms essentially)
                %glc_d(end+1) = k; %add index
                mTeSRmetabos2(end+1) = k;
            end
        %end
    end
end



                   
 % find reactions with glc_D
 rxnlist = [];
 testmTeSRrxnlist = [];
 mTeSRrightones = {};
 rightones2 = {};


 for k = 1:length(today2) %parse through rxns
     for jj = 1:length(mediavardata1{1,1})
     % check forward direction
        if today4(k) > 0 %check if upper bound flux is positive
         % find glc_d reactant
             reactant_found = false; %set glc_d as a reactant as false
             for m = find(full(today5(:,k)) < 0)' %m is vector with indices of rxns with glc_d as reactant being consumed
                 if any(m == mTeSRmetabvector2{jj,jj}) %returns 1 for all rxns with metab(j) as reactant 
                     reactant_found = true; %set glc_d as a reactant as true
                     break;
                 end
             end
 
         % find product
             product_found = false;
             for m = find(full(today5(:,k)) > 0)'
                 if any(m == mTeSRmetabvector2{jj,jj})
                     product_found = true;
                     break;
                 end
             end
             % add
             if reactant_found && product_found% && ~any(k == testmTeSRrxnlist)  non zeros must be present for all three scenarios model.rxn
                 testmTeSRrxnlist(jj,end+1) = k;
             end
         end
     
         % check reverse direction
         if today4(k) < 0
             % find reactant
             reactant_found = false;
             for m = find(full(today5(:,k)) > 0)'
                 if any(m == mTeSRmetabvector2{jj,jj})
                     reactant_found = true;
                     break;
                 end
             end
             % find product
             product_found = false;
             for m = find(full(today5(:,k)) < 0)'
                 if any(m == mTeSRmetabvector2{jj,jj})
                     product_found = true;
                     break;
                 end
             end
             % add
             if reactant_found && product_found% && ~any(k == testmTeSRrxnlist)
                 testmTeSRrxnlist(jj,end+1) = k;
             end
         end
     end
 end

glcfull = testmTeSRrxnlist(12,:);
glcfullsparase = sparse(glcfull);
glcindex = nonzeros(glcfullsparase)';

finalgroupindices = {};
metabsparse = [];

for i = 1:length(mediavardata1{1,1})
    metabsparse(i,:) = sparse(testmTeSRrxnlist(i,:));
    finalgroupindices{i} = nonzeros(metabsparse(i,:))';
end


