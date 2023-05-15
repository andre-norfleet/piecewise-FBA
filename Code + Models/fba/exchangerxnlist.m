function[wholeexchrxnlist1, exchlistallindex] = exchangerxnlist()
load('Z:\Andre\From Andre\FBA-pipeline-master-Riya\Code + Models\data\recon\recon3d_qflux.mat');

exchlistallindex = [];
% restrict input exchange reactions to those in media
    for i = 1:length(model.rxns)
        if length(model.rxns{i}) >= 3
            if strcmp(model.rxns{i}(1:3),'EX_') % if the first 3 characters of the rxn descriptive string ID contain exchange denotation
                exchlistallindex(end+1) = i;
            end
        end
    end

    wholeexchrxnlist1 = model.rxns(exchlistallindex(:));
