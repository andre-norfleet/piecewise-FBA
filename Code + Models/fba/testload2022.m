load('Z:/Andre/From Andre/FBA-pipeline-master-Riya/Code + Models/data/recon/recon3d_qflux.mat');
model2 = model;
B = model.b;
rxnvector = model.rxns;
metabos = model2.mets;
%B = model.mets;
%B = model.b;
Index = find(contains(rxnvector,'GLCt1r'));
indtrbshoot = find(contains(rxnvector,'glc'))
Index2 = find(contains(rxnvector,'EX_o2[e]'));
Index105 = find(contains(rxnvector,'EX_gln_L[e]'));
ascorbate = model.rxns{9017,1};
rightone = model.rxns{1854,1};
rightone2 = model.rxns{Index2,1};
%rightone3 = strcmp(model.rxns,'EX_ascb_L[e]')

%model22 = addReaction(model2, 'InsPractice', 'reactionFormula', 'M03170[e] -> M03170[c]');
%model22b = addReaction(model2, 'InsPractice', 'reactionFormula', 'Ins[e] -> Ins[c]');
%model22c = addReaction(model22b, 'EX_Ins[e]', 'reactionFormula', 'Ins[] -> Ins[e]');
%metabos2 = model22b.mets;
%Index3 = find(contains(metabos2,'Ins[e]'));

%save('Z:\Riya\From Andre\FBA-pipeline-master-Riya\Code + Models\data\recon\recon3d_qflux22.mat','model22')
%save('Z:\Riya\From Andre\FBA-pipeline-master-Riya\Code + Models\data\recon\recon3d_qflux22c.mat','model22c')
