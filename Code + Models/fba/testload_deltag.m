    load('Y:\Riya\From Andre\FBA-pipeline-master-Riya\Code + Models\data\recon\recon3d_qflux.mat');
    input_file = 'input/WTC11testconc/_TEMPLATE_1.xlsx';
    [~,~,raw] = xlsread(input_file,'Concentrations');

    % recalculate thermodynamics
    row = 3; %basically if yes was chosen for concentrations 
    value = raw(row,1);
    value = value{1};
    if ~isnan(value)
        if ischar(value)
            if strcmpi(value,'X')
                recalculate_thermodynamics = 1;
            else
                error('ERROR - Concentrations - Recalculate Thermodynamics - Value must either be X or blank')
            end
        else
            error('ERROR - Concentrations - Recalculate Thermodynamics - Value must either be X or blank')
        end
    else
        recalculate_thermodynamics = 0;
    end

    row = 4;
    value = raw(row,1);
    value = value{1};
    if ~isnan(value)
        if ischar(value)
            if strcmpi(value,'X')
                if recalculate_thermodynamics == 1;
                   error('ERROR - Concentrations - Recalculate Thermodynamics - Cannot select both Yes and No') 
                end
            else
                error('ERROR - Concentrations - Recalculate Thermodynamics - Value must either be X or blank')
            end
        else
            error('ERROR - Concentrations - Recalculate Thermodynamics - Value must either be X or blank')
        end
    elseif recalculate_thermodynamics == 0;
        error('ERROR - Concentrations - Recalculate Thermodynamics - Must select either Yes and No')
    end

        concentration_ranges_all = {};
    for row = 13:size(raw,1)
        value = raw(row,2);
        value = value{1};

        if ~isnan(value)
            if ischar(value)
                concentration_ranges_all{end+1} = value;
            else
                error('ERROR - Concentrations - Ranges Applied to All Samples - File name must be valid string')
            end
        end
    end
    
    % concentration values/ranges applied to individual samples
    concentration_ranges_sample_folder = {};
    concentration_ranges_sample = {};
    for row = 13:size(raw,1)

        value = raw(row,4);
        value = value{1};
        if ~isnan(value)
            if ischar(value)
            
                % folder name
                concentration_ranges_sample_folder{end+1} = value;
                
                % sample name
                value = raw(row,5);
                value = value{1};
                if ~isnan(value)
                    if ischar(value)
                        concentration_ranges_sample{end+1} = value;
                    else
                        error('ERROR - Concentrations - Ranges Applied to Individual Samples - Sample name must be valid string')
                    end
                else
                    error('ERROR - Concentrations - Ranges Applied to Individual Samples - Must provide sample name')
                end       
            else
                error('ERROR - Concentrations - Ranges Applied to Individual Samples - Folder name must be valid string')
            end
        end
    end
    
    % default concentration ranges - lower bound
    value = raw(8,1);
    value = value{1};
    if recalculate_thermodynamics == 1
        if ~isnan(value)
            if isfloat(value)
                if value >= 0
                    default_concentration_lb = value;
                else
                    error('ERROR - Concentrations - Default Lower Bound - Value must be >= 0')
                end
            else
                error('ERROR - Concentrations - Default Lower Bound - Value must be a number')
            end
        else
            error('ERROR - Concentrations - Default Lower Bound - Must provide value if recalculating thermodynamics')
        end
    else
        default_concentration_lb = NaN;
    end
    
    % default concentration ranges - upper bound
    value = raw(9,1);
    value = value{1};
    if recalculate_thermodynamics == 1
        if ~isnan(value)
            if isfloat(value)
                if value >= 0
                    if value >= default_concentration_lb
                        default_concentration_ub = value;
                    else
                        error('ERROR - Concentrations - Default Upper Bound - Value must be >= default lower bound')
                    end
                else
                    error('ERROR - Concentrations - Default Upper Bound - Value must be >= 0')
                end
            else
                error('ERROR - Concentrations - Default Upper Bound - Value must be a number')
            end
        else
            error('ERROR - Concentrations - Default Upper Bound - Must provide value if recalculating thermodynamics')
        end
    else
        default_concentration_ub = NaN;
    end        

conc_all_met = {};
        conc_all_lower = [];
        conc_all_upper = [];
        for i = 1:length(concentration_ranges_all)
        
            % load data
            f = fopen(sprintf('../data/concentration/_ALL_/%s.csv',concentration_ranges_all{i}),'r');
            data = textscan(f,'%s %f %f','Delimiter',',','headerLines',1);
            fclose(f);
            
            % add metabolite id, lower and upper bound
            for j = 1:length(data{1})
                
                % if metabolite already added
                if any(strcmp(conc_all_met,data{1}{j}))
                
                    % change lower bound if lower than original
                    if data{2}(j) < conc_all_lower(strcmp(conc_all_met,data{1}{j}))
                        conc_all_lower(strcmp(conc_all_met,data{1}{j})) = data{2}(j);
                    end
                    
                    % change upper bound if greater than original
                    if data{3}(j) > conc_all_upper(strcmp(conc_all_met,data{1}{j}))
                        conc_all_upper(strcmp(conc_all_met,data{1}{j})) = data{3}(j);
                    end  
                
                % if metabolite not already added
                else
                
                    % add info
                    conc_all_met{end+1} = data{1}{j};
                    conc_all_lower(end+1) = data{2}(j);
                    conc_all_upper(end+1) = data{3}(j);        
                end
            end
        end
        
        % load datasets for individual samples
        conc_individual_sample = {};
        conc_individual_dataset = {};
        conc_individual_met = {};
        conc_individual_lower = {};
        conc_individual_upper = {};
        for i = 1:length(concentration_ranges_sample)
        
            % add "all" concentrations to individual sample list
            conc_individual_sample{end+1} = concentration_ranges_sample{i};
            conc_individual_dataset{end+1} = concentration_ranges_sample_folder{i};
            conc_individual_met{end+1} = conc_all_met;
            conc_individual_lower{end+1} = conc_all_lower;
            conc_individual_upper{end+1} = conc_all_upper;
                
            % load data
            f = fopen(sprintf('Y:/Riya/From Andre/FBA-pipeline-master-Riya/Code + Models/data/concentration/%s/%s.csv',concentration_ranges_sample_folder{i},concentration_ranges_sample{i}),'r');
            data = textscan(f,'%s %f %f','Delimiter',',','headerLines',1);
            fclose(f);
                
            % add metabolite id, lower and upper bound
            for j = 1:length(data{1})

                % if metabolite already added
                if any(strcmp(conc_individual_met{end},data{1}{j}))

                    % change lower bound if lower than original
                    if data{2}(j) < conc_individual_lower{end}(strcmp(conc_individual_met{end},data{1}{j}))
                        conc_individual_lower{end}(strcmp(conc_individual_met{end},data{1}{j})) = data{2}(j);
                    end

                    % change upper bound if greater than original
                    if data{3}(j) > conc_individual_upper{end}(strcmp(conc_individual_met{end},data{1}{j}))
                        conc_individual_upper{end}(strcmp(conc_individual_met{end},data{1}{j})) = data{3}(j);
                    end    

                % if metabolite not already added
                else

                    % add info
                    conc_individual_met{end}{end+1} = data{1}{j};
                    conc_individual_lower{end}(end+1) = data{2}(j);
                    conc_individual_upper{end}(end+1) = data{3}(j);
                end
            end
        end
        
        % load deltaG file
        f = fopen('../data/deltag/deltag.csv','r');
        data = textscan(f,'%s %f %f %f %f %f','Delimiter',',','headerLines',1);
        fclose(f);
        
        % calculate deltag ranges based on concentrations
        deltag_lower = -Inf*ones(length(model.rxns),length(samples_all));
        deltag_upper = Inf*ones(length(model.rxns),length(samples_all));
        
        % implement concentrations for all samples
        for i = 1:length(model.rxns)
            
            % if deltag data available
            if any(strcmp(data{1},model.rxns{i}))
            
                % get concentrations of reactants and products
                reactant_power = [];
                reactant_lower = [];
                reactant_upper = [];
                product_power = [];
                product_lower = [];
                product_upper = [];

                % iterate over reactants
                for j = find(model.S(:,i) < 0)'
                    reactant_power = full(model.S(j,i));

                    % if metabolite found in concentration data
                    if any(strcmp(conc_all_met,model.mets{j}))
                        reactant_lower(end+1) = conc_all_lower(strcmp(conc_all_met,model.mets{j}));
                        reactant_upper(end+1) = conc_all_upper(strcmp(conc_all_met,model.mets{j}));

                    % if metabolite not found in concentration data
                    else
                        reactant_lower(end+1) = default_concentration_lb;
                        reactant_upper(end+1) = default_concentration_ub;
                    end
                end

                % iterate over products
                for j = find(model.S(:,i) > 0)'
                    product_power = full(model.S(j,i));

                    % if metabolite found in concentration data
                    if any(strcmp(conc_all_met,model.mets{j}))
                        reactant_lower(end+1) = conc_all_lower(strcmp(conc_all_met,model.mets{j}));
                        reactant_upper(end+1) = conc_all_upper(strcmp(conc_all_met,model.mets{j}));

                    % if metabolite not found in concentration data
                    else
                        reactant_lower(end+1) = default_concentration_lb;
                        reactant_upper(end+1) = default_concentration_ub;
                    end
                end

                % calculate deltag bounds
                Q_max = prod(product_upper.^product_power)/prod(reactant_lower.^reactant_power);
                Q_min = prod(product_lower.^product_power)/prod(reactant_upper.^reactant_power);
                deltag_lower(i,:) = data{4}(strcmp(data{1},model.rxns{i})) - 8.3144598/1000*310.15*log(Q_max);
                deltag_upper(i,:) = data{4}(strcmp(data{1},model.rxns{i})) - 8.3144598/1000*310.15*log(Q_min);
                %disp(Q_max)
            end
        end
       
        
        % implement concentrations for individual samples
        for a = 1:length(conc_individual_sample)
        
            % implement concentrations for all samples
            for i = 1:length(model.rxns)

                % if deltag data available
                if any(strcmp(data{1},model.rxns{i}))

                    % get concentrations of reactants and products
                    reactant_power = [];
                    reactant_lower = [];
                    reactant_upper = [];
                    product_power = [];
                    product_lower = [];
                    product_upper = [];

                    % iterate over reactants
                    for j = find(model.S(:,i) < 0)'
                        reactant_power = full(model.S(j,i));

                        % if metabolite found in concentration data
                        if any(strcmp(conc_individual_met{a},model.mets{j}))
                            reactant_lower(end+1) = conc_individual_lower{a}(strcmp(conc_individual_met{a},model.mets{j}));
                            reactant_upper(end+1) = conc_individual_upper{a}(strcmp(conc_individual_met{a},model.mets{j}));

                        % if metabolite not found in concentration data
                        else
                            reactant_lower(end+1) = default_concentration_lb;
                            reactant_upper(end+1) = default_concentration_ub;
                        end
                    end

                    % iterate over products
                    for j = find(model.S(:,i) > 0)'
                        product_power = full(model.S(j,i));

                        % if metabolite found in concentration data
                        if any(strcmp(conc_individual_met{a},model.mets{j}))
                            reactant_lower(end+1) = conc_individual_lower{a}(strcmp(conc_individual_met{a},model.mets{j}));
                            reactant_upper(end+1) = conc_individual_upper{a}(strcmp(conc_individual_met{a},model.mets{j}));

                        % if metabolite not found in concentration data
                        else
                            reactant_lower(end+1) = default_concentration_lb;
                            reactant_upper(end+1) = default_concentration_ub;
                        end
                    end

                    % calculate deltag bounds
                    Q_max = prod(product_upper.^product_power)/prod(reactant_lower.^reactant_power);
                    Q_min = prod(product_lower.^product_power)/prod(reactant_upper.^reactant_power);
                    deltag_lower(i,strcmp(samples_all,conc_individual_sample{a})) = data{4}(strcmp(data{1},model.rxns{i})) - 8.3144598/1000*310.15*log(Q_max)
                    deltag_upper(i,strcmp(samples_all,conc_individual_sample{a})) = data{4}(strcmp(data{1},model.rxns{i})) - 8.3144598/1000*310.15*log(Q_min)
                    
                end
            end
        end