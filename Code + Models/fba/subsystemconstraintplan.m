load('Z:/Andre/From Andre/FBA-pipeline-master-Riya/Code + Models/data/recon/recon3d_qflux.mat');

%Read in subsystems
subsyst = model.subSystems;

rxnvector = model.rxns;
upperboundz = model.ub;

%Parse glyc rxns

%Glycolysis index vector
glycIndex = find(contains(subsyst,'Glycolysis'));

%Extract IDs
testIDvect = {};
for i = 1:length(glycIndex)
    testIDvect{end+1} = model.rxns(glycIndex(i,1)); %returns 1 x i cell
end

%Use to set upper bound
testupperboundz = model.ub; 

for j = 1:length(glycIndex)
     testupperboundz(glycIndex(j,1)) = 133;
end

