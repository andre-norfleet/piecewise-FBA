function [] = run_fba333(input_folder,output_folder,input_file,model,type_of_analysis,samples_all,samples_all_source,samples_all_source_protein,vmax_forward,vmax_reverse,objective_function_direction,objective_function_weight,objective_function_reaction,exchange_rxn,deltag_lower,deltag_upper,constraints_fraction_count,constraints_fraction_type,constraints_fraction_reaction,constraints_fraction_value,constraints_fraction_id,constraints_maxmulti_fraction,constraints_maxmulti_count,pfba_fraction,fva_all_fraction,fva_select_fraction,fva_select_groups,fva_select_reactions,sampler_fraction,sampler_samples,sampler_steps,media_sensitivity_critical_fraction,media_sensitivity_critical_metabolite1,media_sensitivity_spectrum_fraction,media_sensitivity_spectrum_metabolite1,media_sensitivity_spectrum_metabolite2,media_sensitivity_spectrum_points,parallel_cores,knockdown_genes,knockdown_fractions,catalase_version,module_5fu,module_cis,module_cpa,module_dox);
%added module_glc 07/11/22 by Riya
%input_file = 'input/WTC11mTeSRmedsensATPm5/_TEMPLATE_1.xlsx';
%output_folder = 'WTC11mTeSRmedsensATPm5';
    % initialize results
    temp = strsplit(input_file,'.xlsx');
    temp = strsplit(temp{1},'/');
    input_filename = temp{end};
    mkdir(sprintf('results/%s/%s',output_folder,input_filename))
    if type_of_analysis == 2
        mkdir(sprintf('results/%s/%s/pfba',output_folder,input_filename))
    elseif (type_of_analysis == 3) || (type_of_analysis == 4)
        mkdir(sprintf('results/%s/%s/fva',output_folder,input_filename))
    elseif type_of_analysis == 5
        mkdir(sprintf('results/%s/%s/sampling',output_folder,input_filename))
    elseif type_of_analysis == 7
        mkdir(sprintf('results/%s/%s/media_sensitivity',output_folder,input_filename))
    end
    
    % create objective value file
    f_obj = fopen(sprintf('results/%s/%s/objval.tsv',output_folder,input_filename),'w');
    fprintf(f_obj,'SAMPLE\tSOURCE\tOBJVAL\n');
        
    % create critical media sensitivity file
    if type_of_analysis == 6
        f_media = fopen(sprintf('results/%s/%s/media_sensitivity_critical.tsv',output_folder,input_filename),'w');
        fprintf(f_media,'SAMPLE\tSOURCE\tCRITICAL UPTAKE [mmol/gDW/hr]\n');
    end

    % save initial model object
    modelOriginal = model;
    
    % iterate over each sample
    for i = 1:length(samples_all)
        
        % load original model
        model = modelOriginal;

        % implement vmax values
        for j = 1:length(model.rxns)

            % forward value
            if (model.ub(j) == max(model.ub)) && (~(isnan(vmax_forward(j,i))))
                model.ub(j) = vmax_forward(j,i);
            end

            % reverse value
            if (model.lb(j) == min(model.lb)) && (~(isnan(vmax_reverse(j,i))))
                model.lb(j) = -vmax_reverse(j,i);
            end   
        end

        % implement deltag constraints
        for j = 1:length(model.rxns)

            % if deltag range only positive, restrict forward direction
            if deltag_lower(j,i) > 0
                model.ub(j) = 0;
            end

            % if deltag range only negative, restrict reverse direction (if reaction is reversible)
            if deltag_upper(j,i) < 0
                model.lb(j) = 0;
            end
        end
        
        % implement catalase
        catalase_reactions = {'r0010','CATm','CATp','CATr'};
        if catalase_version == 0
            for j = 1:length(catalase_reactions)
                model.lb(strcmp(model.rxns,strcat(catalase_reactions{j},'_nadph'))) = 0;
                model.ub(strcmp(model.rxns,strcat(catalase_reactions{j},'_nadph'))) = 0;
            end
        else
            for j = 1:length(catalase_reactions)
                model.lb(strcmp(model.rxns,catalase_reactions{j})) = 0;
                model.ub(strcmp(model.rxns,catalase_reactions{j})) = 0;
            end
        end
        %ins = 3; %setting insulin concentration (uM)
        %insvalue = xlsread('Riya_2022.csv','Riya_2022','B136');
        %glucose_uptake_rate = (50-0)/(1+exp((3-insvalue)/(6e-3))); %sets gluc uptake rate based on ins conc
        %glcrxnindex = customextracellcytosoltransport();
        %glcrxnlist = model.rxns(glcrxnindex(:));
        %glut1rxn = {'GLCt1r'};
        %for j = 1:length(glcrxnlist)
        %    model.lb(strcmp(model.rxns, glcrxnlist(j))) = 0; %lower bound
        %    model.ub(strcmp(model.rxns, glcrxnlist(j))) = 0; %07/12/22 upper bound changed to 0 instead of glucose_uptake_rate
        %end

        %mTeSRexchindex = customexchangetransportDRE();
        %mTESEexchrxnlist = model.rxns(RPMIexchindex(:));
        %mTeSRconcvec = customRPMIconcDRE();

        %for j = 1:length(mTeSRexchrxnlist)
        %    model.lb(strcmp(model.rxns, mTeSRexchrxnlist(j))) = 0; %lower bound, no uptake
        %    model.ub(strcmp(model.rxns, mTeSRexchrxnlist(j))) = mTeSRconcvec(j); %07/29/22 upper bound changed to media concentration
        %end

        %RPMIexchindex = customexchangetransportRPMIDRE();
        %RPMIexchrxnlist = model.rxns(RPMIexchindex(:));
        %RPMIconcvec = customRPMIconcDRE();

        %for j = 1:length(RPMIexchrxnlist)
        %    model.lb(strcmp(model.rxns, RPMIexchrxnlist(j))) = 0; %lower bound, no uptake
        %    model.ub(strcmp(model.rxns, RPMIexchrxnlist(j))) = RPMIconcvec(j); %07/29/22 upper bound changed to media concentration
        %end
        
        %GFAMexchindex = customexchangetransportGFAMDRE(); %(1x61 double) contain indices
        %GFAMexchrxnlist = model.rxns(GFAMexchindex(:));
        %GFAMconcvec = customGFAMconcDRE();

        %for j = 1:length(GFAMexchrxnlist)
        %    model.lb(strcmp(model.rxns, GFAMexchrxnlist(j))) = 0; %lower bound, no uptake
        %    model.ub(strcmp(model.rxns, GFAMexchrxnlist(j))) = GFAMconcvec(j); %07/29/22 upper bound changed to media concentration
        %end

        % Dre Novel Media-based Constraint Alg
        [fincustomconcvaluevector, finalgroupindices, medianumbers, exchrxnindex] = mediaconstraintsDre(input_folder);
        custexchrxnlistnew = model.rxns(exchrxnindex(:))
        custtransportrxnindices = {};

        for j = 1:length(custexchrxnlistnew)
            custtransportrxnindices{j} = model.rxns(finalgroupindices{1,j});
        end

        for j = 1:length(custexchrxnlistnew)
            model.lb(strcmp(model.rxns, custexchrxnlistnew(j))) = 0; %lower bound, no uptake
            model.ub(strcmp(model.rxns, custexchrxnlistnew(j))) = fincustomconcvaluevector(j); %07/29/22 upper bound changed to media concentration
            custtransportrxnindices{j} = model.rxns(finalgroupindices{1,j});
            model.ub(strcmp(model.rxns, custtransportrxnindices(j))) = fincustomconcvaluevector(j); 
        end

        % implement 5-fluorouracil module
        if module_5fu == 0
            for j = 1:length(model.rxns)
                if strcmp(model.subSystems{j},'MODULE 5-FLUOROURACIL');
                    model.lb(j) = 0;
                    model.ub(j) = 0;
                end
            end
        end
        
        % implement cisplatin module
        if module_cis == 0
            for j = 1:length(model.rxns)
                if strcmp(model.subSystems{j},'MODULE CISPLATIN');
                    model.lb(j) = 0;
                    model.ub(j) = 0;
                end
            end
        end
        
        % implement cyclophosphamide module
        if module_cpa == 0
            for j = 1:length(model.rxns)
                if strcmp(model.subSystems{j},'MODULE CYCLOPHOSPHAMIDE');
                    model.lb(j) = 0;
                    model.ub(j) = 0;
                end
            end
        end
        
        % implement doxorubicin module
        if module_dox == 0
            for j = 1:length(model.rxns)
                if strcmp(model.subSystems{j},'MODULE DOXORUBICIN');
                    model.lb(j) = 0;
                    model.ub(j) = 0;
                end
            end
        end
        
        % implement glucose module
%        if module_glc == 0 %added glc 07/11/22
%            for j = 1:length(model.rxns)
%                if strcmp(model.subSystems{j},'MODULE GLUCOSE');
%                    model.lb(j) = 0;
%                    model.ub(j) = 0;
%                end
%            end
%        end
        
        % setup LP problem
        model.A = model.S;
        model.rhs = model.b;
        model.obj = model.c;
        model.sense = repmat('=',1,length(model.mets));
        model.vtype = repmat('C',1,length(model.rxns));
        model.varnames = model.rxns;
        model.modelsense = lower(objective_function_direction);

        % reset objective value
        model.c = zeros(size(model.c));
        model.obj = zeros(size(model.obj));

        % find fraction of maximum constraints
        max_constraint_values = [];
        for j = 1:constraints_fraction_count

            % implement objective value
            model.c(model.c_fraction == j) = 1;
            model.obj(model.c_fraction == j) = 1;
            model.modelsense = 'max';

            % get objective value
            params.outputflag = 0;
            result = gurobi(model,params);
            max_constraint_values(end+1) = result.objval;

            % reset objective value
            model.c = zeros(size(model.c));
            model.obj = zeros(size(model.obj));
        end

        % implement fraction of maximum constraints
        for j = 1:constraints_fraction_count

            % equal to
            if strcmp(constraints_fraction_type{j},'equal')
                model.lb(constraints_fraction_id(j)) = constraints_fraction_value{j}(1) * max_constraint_values(j);
                model.ub(constraints_fraction_id(j)) = constraints_fraction_value{j}(1) * max_constraint_values(j);

            % less than
            elseif strcmp(constraints_fraction_type{j},'less')
                model.ub(constraints_fraction_id(j)) = constraints_fraction_value{j}(1) * max_constraint_values(j);

            % greater than
            elseif strcmp(constraints_fraction_type{j},'greater')
                model.lb(constraints_fraction_id(j)) = constraints_fraction_value{j}(1) * max_constraint_values(j);

            % range
            elseif strcmp(constraints_fraction_type{j},'range')
                if constraints_fraction_value{j}(1) <= constraints_fraction_value{j}(2)
                    model.lb(constraints_fraction_id(j)) = constraints_fraction_value{j}(1) * max_constraint_values(j);
                    model.ub(constraints_fraction_id(j)) = constraints_fraction_value{j}(2) * max_constraint_values(j);
                else
                    error('ERROR - Constraints - Custom - For Range, min value must be <= max value');
                end
            end      
        end  

		% find max multi constraints
		if constraints_maxmulti_count > 0
		
			% implement objective value
			model.c = model.c_maxmulti;
			model.obj = model.c_maxmulti;
			model.modelsense = 'max';

            % get objective value
            params.outputflag = 0;
            result = gurobi(model,params);
            maxmulti_value = result.objval;

            % reset objective value
            model.c = zeros(size(model.c));
            model.obj = zeros(size(model.obj));
        end

        % implement objective function
        model.c = model.c_obj;
        model.obj = model.c_obj;
        model.modelsense = 'max';
		
		% implement max multi constraint
		if constraints_maxmulti_count > 0
			model.S(end+1,:) = model.c_maxmulti;
			model.b(end+1) = constraints_maxmulti_fraction*maxmulti_value;
			model.csense(end+1) = 'E';
			model.mets{end+1} = 'maxmulti';
			model.metCharges(end+1) = 0;
			model.metFormulas{end+1} = '';
			model.metSmiles{end+1} = '';
			model.metNames{end+1} = 'Max Multi Constraint';
			model.metHMDBID{end+1} = '';
			model.metInChIString{end+1} = '';
			model.metKEGGID{end+1} = '';
			model.metPubChemID{end+1} = '';
			model.metCHEBIID{end+1} = '';
			model.metPdMap{end+1} = '';
			model.metReconMap{end+1} = '';
		end
		
		% create parameters
		model.A = model.S;
		model.rhs = model.b;
		model.sense = repmat('=',1,length(model.csense));
		for j = 1:length(model.csense)
			if model.csense(i) == 'G'
				model.sense(i) = '>';
			elseif model.csense(i) == 'L'
				model.sense(i) = '<';
			end
		end
		
        % objective function value
        params.outputflag = 0;
        result = gurobi(model,params);
        fprintf(f_obj,'%s\t%s\t%0.9f\n',samples_all{i},samples_all_source_protein{i},result.objval);

        %%% TROUBLESHOOTING

        % media sensitivity - spectrum
        if type_of_analysis == 7

            % ensure sum(model.c) = 1
            %if sum(model.c) ~= 1
                %error('ERROR - General - Media Sensitivity - Can only perform media sensitivity with one objective function reaction')

            % ensure only one maximization objective
                %if length(objective_function_reaction) ~= 1
                %    error('ERROR - General - Type of Analysis - Can only perform media sensitivity with maximizing one objective function reaction');
                %elseif ~strcmp(objective_function_direction{1},'MAX')
                %    error('ERROR - General - Type of Analysis - Can only perform media sensitivity with maximizing one objective function reaction');
                %end

            % check metabolite 1 name
            if ~any(strcmp(model.mets,sprintf('%s[e]',media_sensitivity_spectrum_metabolite1)))
                error('ERROR - General - Media Sensitivity Options - Metabolite - Metabolite 1 not found in extracellular compartment')

            % check if exchange reaction exists
            elseif ~any(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1)))
                error('ERROR - General - Media Sensitivity Options - Metabolite - Metabolite 1 not involved in exchange reaction')

            % ensure that exchange reaction is in exchange_rxn
            elseif ~any(strcmp(exchange_rxn,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1)))
                error('ERROR - General - Media Sensitivity Options - Metabolite - Exchange reaction for metabolite 1 not found in media composition')
            end  

            if ~isnan(media_sensitivity_spectrum_metabolite2)

                % check metabolite 2 name
                if ~any(strcmp(model.mets,sprintf('%s[e]',media_sensitivity_spectrum_metabolite2)))
                    error('ERROR - General - Media Sensitivity Options - Metabolite - Metabolite 2 not found in extracellular compartment')

                % check if exchange reaction exists
                elseif ~any(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2)))
                    error('ERROR - General - Media Sensitivity Options - Metabolite - Metabolite 2 not involved in exchange reaction')

                % ensure that exchange reaction is in exchange_rxn
                elseif ~any(strcmp(exchange_rxn,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2)))
                    error('ERROR - General - Media Sensitivity Options - Metabolite - Exchange reaction for metabolite 2 not found in media composition')
                end
            end

            % if 1 metabolite (end 569)
            if isnan(media_sensitivity_spectrum_metabolite2)

                % find minimum exchange value while maintaining objective value
                original_objective = find(model.c);
                original_objective_value = result.objval;
                original_objective_lb = model.lb(original_objective);
                model.lb(original_objective) = media_sensitivity_spectrum_fraction*original_objective_value;
                model.c = zeros(size(model.c));
                model.c(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))) = 1;

                model.A = model.S;
                model.rhs = model.b;
                model.obj = model.c;
                model.sense = repmat('=',1,length(model.mets));
                model.vtype = repmat('C',1,length(model.rxns));
                model.varnames = model.rxns;
                model.modelsense = 'max';

                params.outputflag = 0;
                result = gurobi(model,params);
                start_val = result.objval

                % if critical value is negative
                if start_val < 0

                    % return original objective function
                    model.c = zeros(size(model.c));
                    model.c(original_objective) = 1;
                    model.lb(original_objective) = original_objective_lb;

                    % initiailize output file
                    f_media = fopen(sprintf('results/%s/%s/media_sensitivity/%s.tsv',output_folder,input_filename,samples_all{i}),'w');
                    fprintf(f_media,'%s UPTAKE LB [mmol/gDW/hr]\tOBJVAL\n',sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1));

                    % progressively increase flux lb
                    uptake = [];
                    objval = [];
                    for j = 1:media_sensitivity_spectrum_points

                        % set lower bound
                        model.lb(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))) = start_val*(media_sensitivity_spectrum_points-j)/(media_sensitivity_spectrum_points-1); %Dre change media)

                        % find objective value
                        model.A = model.S;
                        model.rhs = model.b;
                        model.obj = model.c;
                        model.sense = repmat('=',1,length(model.mets));
                        model.vtype = repmat('C',1,length(model.rxns));
                        model.varnames = model.rxns;
                        model.modelsense = 'max';

                        params.outputflag = 0;
                        result = gurobi(model,params);

                        % output objective value
                        fprintf(f_media,'%f\t%f\n',model.lb(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))),result.objval);
                        uptake(end+1) = -model.lb(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1)))
                        objval(end+1) = result.objval;
                    end
                    %uptake

                    % close output file
                    fclose(f_media);

                    % create image
                    figure; hold on;
                    plot(uptake,objval,'k-');
                    plot([0,max(uptake)],[max(objval),max(objval)],'k--');
                    hold off;
                    xlabel(sprintf('%s Uptake Rate [mmol/gDW/hr]',media_sensitivity_spectrum_metabolite1),'Interpreter','none');
                    ylabel('Objective Value');
                    title(sprintf('Critical %s Uptake Rate = %g mmol/gDW/hr',media_sensitivity_spectrum_metabolite1,max(uptake)),'Interpreter','none');
                    saveas(gcf,sprintf('results/%s/%s/media_sensitivity/%s.png',output_folder,input_filename,samples_all{i}));        

                % if critical value is positive or zero
                else

                    % create output file
                    f_media = fopen(sprintf('results/%s/%s/media_sensitivity/%s.tsv',output_folder,input_filename,samples_all{i}),'w');
                    fprintf(f_media,'%s UPTAKE LB [mmol/gDW/hr]\tOBJVAL\n',sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1));    
                    fprintf(f_media,'nan\t%f\n',original_objective_value);
                    fclose(f_media);          
                end

            % if 2 metabolites
            elseif (any(media_sensitivity_spectrum_metabolite1)) && (any(media_sensitivity_spectrum_metabolite2)) && (~any(media_sensitivity_spectrum_metabolite3))

                % find minimum exchanges value while maintaining objective value
                original_objective = find(model.c);
                %peepOForig = find(model.c)
                original_objective_value = result.objval;
                original_objective_lb = model.lb(original_objective);
                model.lb(original_objective) = original_objective_value;
                model.c = zeros(size(model.c));

                % % 1-1 ratio (case 1/3)
                model.c(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))) = 1; % match model.rxns input with exchange rxn for metab of interest
                %bund = model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1);
                %bund2 = strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1));
                model.c(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2))) = 1;

                model.A = model.S;
                model.rhs = model.b;
                model.obj = model.c;
                model.sense = repmat('=',1,length(model.mets));
                model.vtype = repmat('C',1,length(model.rxns));
                model.varnames = model.rxns;
                model.modelsense = 'max';

                params.outputflag = 0;
                result = gurobi(model,params);
                start_val1 = result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1)));
                peepino2 = result.x(strcmp(model.rxns,'EX_o2[e]'))
                start_val2 = result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2)));
                peepinglc2 = result.x(strcmp(model.rxns,'EX_glc_D[e]'))


                % % 10-1 ratio (Case 2/3)
                model.c(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))) = 10;
                model.c(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2))) = 1;

                model.A = model.S;
                model.rhs = model.b;
                model.obj = model.c;
                model.sense = repmat('=',1,length(model.mets));
                model.vtype = repmat('C',1,length(model.rxns));
                model.varnames = model.rxns;
                model.modelsense = 'max';

                params.outputflag = 0;
                result = gurobi(model,params);
                if result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))) < start_val1
                    start_val1 = result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))); % was there improvement in max uptake comp to 1 to 1 case
                end
                %curious = result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1)))
                if result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2))) < start_val2
                    start_val2 = result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2)));
                end

                % % 1-10 ratio (Case 3/3)
                model.c(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))) = 1;
                model.c(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2))) = 10;

                model.A = model.S;
                model.rhs = model.b;
                model.obj = model.c;
                model.sense = repmat('=',1,length(model.mets));
                model.vtype = repmat('C',1,length(model.rxns));
                model.varnames = model.rxns;
                model.modelsense = 'max';

                params.outputflag = 0;
                result = gurobi(model,params);
                if result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))) < start_val1
                    start_val1 = result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1)));
                end
                if result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2))) < start_val2
                    start_val2 = result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2)));
                end
                %start_val1

                % if both critical values are negative
                if (start_val1 < 0) && (start_val2 < 0)

                    % result original objective function
                    model.c = zeros(size(model.c));
                    model.c(original_objective) = 1;
                    model.lb(original_objective) = original_objective_lb;

                    % initiailize output file
                    f_media = fopen(sprintf('results/%s/%s/media_sensitivity/%s.tsv',output_folder,input_filename,samples_all{i}),'w');
                    fprintf(f_media,'%s UPTAKE LB [mmol/gDW/hr]\t%s UPTAKE LB [mmol/gDW/hr]\tOBJVAL\n',sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1),sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2));

                    % progressively increase flux lb
                    uptake1 = [];
                    uptake2 = [];
                    objval = [];
                    %binoculars = [];
                    for j = 1:media_sensitivity_spectrum_points
                        for k = 1:media_sensitivity_spectrum_points

                            % set lower bound done prior to sim so obj func
                            % val can be assessed after
                            model.lb(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))) = start_val1*(media_sensitivity_spectrum_points-j)/(media_sensitivity_spectrum_points-1);
                            %binoculars(end+1) = start_val1*(media_sensitivity_spectrum_points-j)/(media_sensitivity_spectrum_points-1);
                            model.lb(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2))) = start_val2*(media_sensitivity_spectrum_points-k)/(media_sensitivity_spectrum_points-1);

                            % find objective value
                            model.A = model.S;
                            model.rhs = model.b;
                            model.obj = model.c;
                            model.sense = repmat('=',1,length(model.mets));
                            model.vtype = repmat('C',1,length(model.rxns));
                            model.varnames = model.rxns;
                            model.modelsense = 'max';

                            params.outputflag = 0;
                            result = gurobi(model,params);

                            % output objective value
                            fprintf(f_media,'%f\t%f\t%f\n',model.lb(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))),model.lb(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2))),result.objval);
                            uptake1(end+1) = -model.lb(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1)));
                            uptake2(end+1) = -model.lb(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2)));
                            objval(end+1) = result.objval;
                        end
                    end            
                    %print(uptake1)
                    % close output file
                    fclose(f_media);

                    % create image
                    Z = reshape(objval,media_sensitivity_spectrum_points,media_sensitivity_spectrum_points);
                    uptake2 = uptake2;
                    figure;
                    imagesc(uptake1,uptake2,Z);
                    xlabel(sprintf('%s Uptake Rate [mmol/gDW/hr]',media_sensitivity_spectrum_metabolite1),'Interpreter','none');
                    ylabel(sprintf('%s Uptake Rate [mmol/gDW/hr]',media_sensitivity_spectrum_metabolite2),'Interpreter','none');
                    c = colorbar;
                    c.Label.String = 'Objective Value';
                    c.Label.FontSize = 12;
                    set(gca,'YDir','normal');
                    saveas(gcf,sprintf('results/%s/%s/media_sensitivity/%s.png',output_folder,input_filename,samples_all{i}));

                % if one critical value is not negative
                else

                    % create output file
                    f_media = fopen(sprintf('results/%s/%s/media_sensitivity/%s.tsv',output_folder,input_filename,samples_all{i}),'w');
                    fprintf(f_media,'%s UPTAKE LB [mmol/gDW/hr]\t%s UPTAKE LB [mmol/gDW/hr]\tOBJVAL\n',sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1),sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2));
                    fprintf(f_media,'nan\tnan\t%f\n',original_objective_value);
                    fclose(f_media);
                end
            
            else % 3 Metabolites Dre
                original_objective = find(model.c);
                original_objective_value = result.objval;
                original_objective_lb = model.lb(original_objective);
                model.lb(original_objective) = original_objective_value;
                model.c = zeros(size(model.c));

                % Case 1 (1-1-1 ratio)
                model.c(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))) = 1;
                model.c(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2))) = 0;
                model.c(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite3))) = 0;

                model.A = model.S;
                model.rhs = model.b;
                model.obj = model.c;
                model.sense = repmat('=',1,length(model.mets));
                model.vtype = repmat('C',1,length(model.rxns));
                model.varnames = model.rxns;
                model.modelsense = 'max';

                params.outputflag = 0;
                result = gurobi(model,params);
                start_val1 = result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1)))
                start_val2 = result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2)))
                start_val3 = result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite3)))

                % Case 2 (10-1-1 ratio)
                model.c(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))) = 1;
                %model.c(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2))) = 1;
                model.c(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite3))) = 1;

                model.A = model.S;
                model.rhs = model.b;
                model.obj = model.c;
                model.sense = repmat('=',1,length(model.mets));
                model.vtype = repmat('C',1,length(model.rxns));
                model.varnames = model.rxns;
                model.modelsense = 'max';

                params.outputflag = 0;
                result = gurobi(model,params);
                if result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))) < start_val1
                    start_val1 = result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))); % was there improvement in max uptake comp to 1 to 1 case
                end
                if result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2))) < start_val2
                    start_val2 = result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2)));
                end
                if result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite3))) < start_val3
                    start_val3 = result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite3)));
                end

                % Case 3 (1-10-1 ratio)
                %model.c(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))) = 1;
                model.c(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2))) = 1;
                model.c(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite3))) = 1;

                model.A = model.S;
                model.rhs = model.b;
                model.obj = model.c;
                model.sense = repmat('=',1,length(model.mets));
                model.vtype = repmat('C',1,length(model.rxns));
                model.varnames = model.rxns;
                model.modelsense = 'max';

                params.outputflag = 0;
                result = gurobi(model,params);
                if result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))) < start_val1
                    start_val1 = result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))); % was there improvement in max uptake comp to 1 to 1 case
                end
                if result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2))) < start_val2
                    start_val2 = result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2)));
                end
                if result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite3))) < start_val3
                    start_val3 = result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite3)));
                end

                % Case 4 (1-1-10 ratio)
                model.c(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))) = 1;
                model.c(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2))) = 1;
                model.c(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite3))) = 10;

                model.A = model.S;
                model.rhs = model.b;
                model.obj = model.c;
                model.sense = repmat('=',1,length(model.mets));
                model.vtype = repmat('C',1,length(model.rxns));
                model.varnames = model.rxns;
                model.modelsense = 'max';

                params.outputflag = 0;
                result = gurobi(model,params);
                if result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))) < start_val1
                    start_val1 = result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))); % was there improvement in max uptake comp to 1 to 1 case
                end
                if result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2))) < start_val2
                    start_val2 = result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2)));
                end
                if result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite3))) < start_val3
                    start_val3 = result.x(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite3)));
                end
                
                % if ALL 3 critical values are negative
                if (start_val1 < 0) && (start_val2 < 0) && (start_val3 < 0)

                    % result original objective function
                    model.c = zeros(size(model.c));
                    model.c(original_objective) = 1;
                    model.lb(original_objective) = original_objective_lb;

                    % initiailize output file
                    f_media = fopen(sprintf('results/%s/%s/media_sensitivity/%s.tsv',output_folder,input_filename,samples_all{i}),'w');
                    fprintf(f_media,'%s UPTAKE LB [mmol/gDW/hr]\t%s UPTAKE LB [mmol/gDW/hr]\tOBJVAL\n',sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1),sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2));

                    % progressively increase flux lb
                    uptake1 = [];
                    uptake2 = [];
                    uptake3 = [];
                    objval = [];
                    %binoculars = [];
                    for j = 1:media_sensitivity_spectrum_points
                        for k = 1:media_sensitivity_spectrum_points
                            for kk = 1:media_sensitivity_spectrum_points

                                % set lower bound done prior to sim so obj func
                                % val can be assessed after
                                model.lb(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))) = start_val1*(media_sensitivity_spectrum_points-j)/(media_sensitivity_spectrum_points-1);
                                model.lb(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2))) = start_val2*(media_sensitivity_spectrum_points-k)/(media_sensitivity_spectrum_points-1);
                                model.lb(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite3))) = start_val3*(media_sensitivity_spectrum_points-k)/(media_sensitivity_spectrum_points-1);

                                % find objective value
                                model.A = model.S;
                                model.rhs = model.b;
                                model.obj = model.c;
                                model.sense = repmat('=',1,length(model.mets));
                                model.vtype = repmat('C',1,length(model.rxns));
                                model.varnames = model.rxns;
                                model.modelsense = 'max';

                                params.outputflag = 0;
                                result = gurobi(model,params);

                                % output objective value
                                fprintf(f_media,'%f\t%f\t%f\n',model.lb(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1))),model.lb(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2))),result.objval);
                                uptake1(end+1) = -model.lb(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1)));
                                uptake2(end+1) = -model.lb(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2)));
                                uptake3(end+1) = -model.lb(strcmp(model.rxns,sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite3)))
                                objval(end+1) = result.objval;
                            end
                        end
                    end
                        %print(uptake1)
                        % close output file
                        fclose(f_media);

                        % if one critical value is not negative
                    else

                    % create output file
                        f_media = fopen(sprintf('results/%s/%s/media_sensitivity/%s.tsv',output_folder,input_filename,samples_all{i}),'w');
                        fprintf(f_media,'%s UPTAKE LB [mmol/gDW/hr]\t%s UPTAKE LB [mmol/gDW/hr]\tOBJVAL\n',sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite1),sprintf('EX_%s[e]',media_sensitivity_spectrum_metabolite2));
                        fprintf(f_media,'nan\tnan\t%f\n',original_objective_value);
                        fclose(f_media);
                    end
                end


            end % if 1 metab end
            %end
        end % if type of analysis is 7 end
    end