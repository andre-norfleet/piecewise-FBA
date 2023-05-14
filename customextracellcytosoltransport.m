function [rxnlist] = customextracellcytosoltransport()
load('Y:\Riya\From Andre\FBA-pipeline-master-Riya\Code + Models\data\recon\recon3d_qflux.mat');
rxnids = model.rxns;
metids = model.mets;
ub = model.ub;
stoichmatrix = model.S;

%find glc_D
glc_d = []; %saves indices of all glucose metabolites in m direction
for k = 1:length(metids) %parse through all metabolites
    if length(metids{k}) == 8 %metabolite is 8 characters in length (ex: glc_D[e])
        if strcmp(metids{k}(1:6),'glc_D[') %metabolite matches gluc 
            glc_d(end+1) = k; %add index 
        end
    end
end
                   
 % find reactions with glc_D
 rxnlist = []; %list of glc_d transport rxns
 for k = 1:length(rxnids) %parse through rxns 
     % check forward direction
     if ub(k) > 0 %check if upper bound flux is positive -- not sure why we need to check this
         % find glc_d reactant
         reactant_found = false; %set glc_d as a reactant as false
         for m = find(full(stoichmatrix(:,k)) < 0)' %m is vector with indices of rxns with glc_d as reactant being consumed
             if any(m == glc_d) %returns 1 for all rxns with glc_d as reactant 
                 reactant_found = true; %set glc_d as a reactant as true
                 break;
             end
         end
 
         % find product
         product_found = false;
         for m = find(full(stoichmatrix(:,k)) > 0)'
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
     if ub(k) < 0
         % find reactant
         reactant_found = false;
         for m = find(full(stoichmatrix(:,k)) > 0)'
             if any(m == glc_d)
                 reactant_found = true;
                 break;
             end
         end
         % find product
         product_found = false;
         for m = find(full(stoichmatrix(:,k)) < 0)'
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

end

 
 

  

