function[fincustomconcvaluevector, finalgroupindices, medianumbers, exchrxnindex] = mediaconstraintsDre(input_folder)
load('Z:\Andre\From Andre\FBA-pipeline-master-Riya\Code + Models\data\recon\recon3d_qflux.mat');
%input_folder = 'WTC11_D4CtrlFVA_RPMIGFAMMediaANDcustconstraintsTOTAL';


directory = dir(sprintf('input/%s/*.xlsx',input_folder));
input_files = {};
for i = 1:length(directory)
    fn = sprintf('input/%s/%s',input_folder,directory(i).name);
    test = strsplit(fn,'/');
    if ~strcmp(extractBetween(test{end},1,2),'~$')
        temp = strsplit(test{end},'.xlsx');
        %if ~any(strcmp(overlap,temp{1}))
        input_files{end+1} = fn;
        %end
    end
end

input_file = fn;
%medianumbers = length(data{2});

    media_dir = dir('../data/media/output/*.csv');
    available_media = {};
    for i = 1:length(media_dir)
        available_media{end+1} = media_dir(i).name;
    end

    trbsht = 1;
    [~,~,raw] = xlsread(input_file,'Media');

    media_choice = {};
    for row = 9:size(raw,1)
        value = raw(row,1);
        value = value{1};

        if ~isnan(value)
            if ischar(value)
                media_choice{end+1} = value;
            else
                error('ERROR - Media - Medium Makeup - Value must either be X or blank')
            end
        end
    end

    for i = 1:length(media_choice)
        if ~any(strcmp(available_media,sprintf('%s.csv',media_choice{i})))
            error('ERROR - Can''t find media file %s.csv',media_choice{i})
        end
    end

    for i = 1:length(media_choice)

        % load media file
        f = fopen(sprintf('../data/media/output/%s.csv',media_choice{i}),'r'); %actually opening selected media file 'i' and loaidng into Matlab
        data = textscan(f,'%s %s %f %f','Delimiter',',','headerLines',1);
        fclose(f);

        custommetabvector = {};
        for n = 1:length(data{2})
            custommetabvector{n,end+1} = find(contains(model.mets,data{1,1}{n,1}(1:end-3) ));
        end

        custommetabvectorfin = {};

        for n = 1:length(data{2})
            custommetabvectorfin{n,end+1} = custommetabvector{n,n}'; %transposing for downstream index comparison to reactant_found and product_found
        end
    end

    customrxnlist = [];

    for k = 1:length(model.rxns) %parse through rxns
        for jj = 1:length(data{2})
        % check forward direction
            if model.ub(k) > 0 %check if upper bound flux is positive
         % find glc_d reactant
                reactant_found = false; %set glc_d as a reactant as false
                for m = find(full(model.S(:,k)) < 0)' %m is vector with indices of rxns with glc_d as reactant being consumed
                    if any(m == custommetabvectorfin{jj,jj}) %returns 1 for all rxns with metab(j) as reactant 
                        reactant_found = true; %set glc_d as a reactant as true
                        break;
                    end
                end

                %find product
                product_found = false;
                for m = find(full(model.S(:,k)) > 0)'
                    if any(m == custommetabvectorfin{jj,jj})
                        product_found = true;
                        break;
                    end
                end

                %add
                if reactant_found && product_found% && ~any(k == testmTeSRrxnlist)
                    customrxnlist(jj,end+1) = k;
                end
            end

            %reverse direction
            if model.ub(k) < 0
                % find reactant
                reactant_found = false;
                for m = find(full(model.S(:,k)) > 0)'
                    if any(m == custommetabvectorfin{jj,jj})
                        reactant_found = true;
                        break;
                 end
             end
             % find product
             product_found = false;
             for m = find(full(model.S(:,k)) < 0)'
                 if any(m == custommetabvectorfin{jj,jj})
                     product_found = true;
                     break;
                 end
             end

             %add
             if reactant_found && product_found% && ~any(k == testmTeSRrxnlist)
                 customrxnlist(jj,end+1) = k;
             end
            end
        end
    end

finalgroupindices = {};
metabsparse = [];

for i = 1:length(data{2})
    metabsparse(i,:) = sparse(customrxnlist(i,:));
    finalgroupindices{i} = nonzeros(metabsparse(i,:))';
end

customconcvaluevector = [];

for nn = 1:length(data{1,3})
    for nnn = 1:length(data{1,3})
        customconcvaluevector(1,nnn) = data{1,3}(nnn);
    end
end

exchrxnindex = [];

for i = 1:length(model.rxns)
        for kk = 1:length(data{1,2})
            if any(strcmp(data{1,2}{kk},model.rxns{i}))
                exchrxnindex(end+1) = i; %index in model.rxns for all exchange rxns accompanying mTeSR media component
            end
        end
end

fincustomconcvaluevector = customconcvaluevector*1000; % usually mulitiplied by .001 for mM conv
medianumbers = length(data{2});

%custtransportrxnindices = {};

seein2 = model.rxns(finalgroupindices{1,12}); % (66x1 cell)
seein23 = model.rxns(finalgroupindices{1,12}(:)); % (66x1 cell)
seein233 = model.rxns(finalgroupindices{1,1});

%for j = 1:medianumbers
%    custtransportrxnindices{1,j} = model.rxns(finalgroupindices{1,j}(:));
%end

%custtransportrxnindices