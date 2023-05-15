load('Y:\Riya\From Andre\FBA-pipeline-master-Riya\Code + Models\data\recon\recon3d_qflux.mat');
growthRates = zeros(21);
for i = 0:20
	for j = 0:20
		model = changeRxnBounds(model,'EX_glc(e)',-i,'b');
		model = changeRxnBounds(model,'EX_o2(e)',-j,'b');
		FBAsolution = optimizeCbModel(model,'max');
		growthRates(i+1,j+1) = FBAsolution.f;
	end
end
