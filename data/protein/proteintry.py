import numpy as np
import pandas as pd
import os
import pickle
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import seaborn as sns

input_folder = 'WTC11_ATPall33' # ****

if not os.path.isdir('output/%s' % input_folder):
    output_folder = input_folder
    os.mkdir('output/%s' % output_folder)
else:
    val = 1
    made_folder = False
    while not made_folder:
        if not os.path.isdir('output/%s_%d' % (input_folder,val)):
            output_folder = '%s_%d' % (input_folder,val)
            os.mkdir('output/%s' % output_folder)
            made_folder = True
        else:
            val += 1

os.mkdir('output/%s/_files_' % output_folder)
fn_stdout = 'output/%s/_files_/stdout.txt' % output_folder
with open(fn_stdout,'w') as f:
    pass

fn_stderr = 'output/%s/_files_/stderr.txt' % output_folder
with open(fn_stderr,'w') as f:
    pass

found_error = False

# load sheet
df_options = pd.read_excel('input/%s/_OPTIONS_.xlsx' % input_folder, sheet_name='Input', header=None)

# data type
datatype_option1 = str(df_options.loc[1][0])
datatype_option2 = str(df_options.loc[2][0])
datatype_option3 = str(df_options.loc[3][0])
if datatype_option1 == 'X' and datatype_option2 != 'X' and datatype_option3 != 'X':
    datatype = 'rnaseq'
elif datatype_option1 != 'X' and datatype_option2 == 'X' and datatype_option3 != 'X':
    datatype = 'gpl570'
elif datatype_option1 != 'X' and datatype_option2 != 'X' and datatype_option3 == 'X':
    datatype = 'gpl10558'
elif datatype_option1 != 'X' and datatype_option2 != 'X' and datatype_option3 != 'X':
    with open(fn_stderr,'a') as f:
        f.write('ERROR - Options File - Data Type - Must select data type with "X"\n')
    found_error = True
else:
    with open(fn_stderr,'a') as f:
        f.write('ERROR - Options File - Data Type - Can only select one data type\n')
    found_error = True

if not found_error:
    if datatype == 'rnaseq':
        
        # load options
        df_options = pd.read_excel('input/%s/_OPTIONS_.xlsx' % input_folder, sheet_name='RNA-Seq', header=None)
        
        # file name
        options_fn = str(df_options.loc[0][1])
        
        # excel sheet name
        options_excel_sheet = str(df_options.loc[3][2])
        
        # file delimiter
        options_delimiter = str(df_options.loc[6][2])
        
        # rna-seq normalization
        options_rnaseq = str(df_options.loc[9][3])
        
        # header row
        options_header = int(df_options.loc[17][2])
        
        # gene/transcript ID column
        options_geneid_column = int(df_options.loc[24][3])
        
        # gene ID type
        options_geneid_type = str(df_options.loc[27][3])
        
        # samples format
        samplesformat_option1 = str(df_options.loc[35][0])
        samplesformat_option2 = str(df_options.loc[36][0])
        if samplesformat_option1 == 'X' and samplesformat_option2 != 'X':
            samplesformat = 'manual'
        elif samplesformat_option1 != 'X' and samplesformat_option2 == 'X':
            samplesformat = 'automatic'
        elif samplesformat_option1 != 'X' and samplesformat_option2 != 'X':
            with open(fn_stderr,'a') as f:
                f.write('ERROR - Options File - RNA-Seq - Samples Format - Must select samples format with "X"\n')
            found_error = True
        else:
            with open(fn_stderr,'a') as f:
                f.write('ERROR - Options File - RNA-Seq - Samples Format - Can only select one samples format\n')
            found_error = True

if not found_error:
    if datatype == 'gpl570':
        
        # load options
        df_options = pd.read_excel('input/%s/_OPTIONS_.xlsx' % input_folder, sheet_name='GPL570', header=None)
        
        # file extension
        options_extension = str(df_options.loc[0][2])
        if options_extension not in ['.cel','.CEL']:
            with open(fn_stderr,'a') as f:
                f.write('ERROR - Options File - GPL570 - Files Extension - Must be either ".cel" or ".CEL"\n')
            found_error = True
        
        # samples format
        samplesformat_option1 = str(df_options.loc[4][0])
        samplesformat_option2 = str(df_options.loc[5][0])
        if samplesformat_option1 == 'X' and samplesformat_option2 != 'X':
            samplesformat = 'manual'
        elif samplesformat_option1 != 'X' and samplesformat_option2 == 'X':
            samplesformat = 'automatic'
        elif samplesformat_option1 != 'X' and samplesformat_option2 != 'X':
            with open(fn_stderr,'a') as f:
                f.write('ERROR - Options File - GPL570 - Samples Format - Must select samples format with "X"\n')
            found_error = True
        else:
            with open(fn_stderr,'a') as f:
                f.write('ERROR - Options File - GPL570 - Samples Format - Can only select one samples format\n')
            found_error = True

if not found_error:
    if datatype == 'gpl10558':
        
        # load options
        df_options = pd.read_excel('input/%s/_OPTIONS_.xlsx' % input_folder, sheet_name='GPL10558', header=None)
        
        # file name
        options_fn = str(df_options.loc[0][1])
        
        # excel sheet name
        options_excel_sheet = str(df_options.loc[3][2])
        
        # file delimiter
        options_delimiter = str(df_options.loc[6][2])

        # header row
        options_header = int(df_options.loc[9][2])
        
        # probe ID column
        options_probeid_column = int(df_options.loc[16][2])
        
        # samples format
        samplesformat_option1 = str(df_options.loc[20][0])
        samplesformat_option2 = str(df_options.loc[21][0])
        if samplesformat_option1 == 'X' and samplesformat_option2 != 'X':
            samplesformat = 'manual'
        elif samplesformat_option1 != 'X' and samplesformat_option2 == 'X':
            samplesformat = 'automatic'
        elif samplesformat_option1 != 'X' and samplesformat_option2 != 'X':
            with open(fn_stderr,'a') as f:
                f.write('ERROR - Options File - GPL10558 - Samples Format - Must select samples format with "X"\n')
            found_error = True
        else:
            with open(fn_stderr,'a') as f:
                f.write('ERROR - Options File - GPL10558 - Samples Format - Can only select one samples format\n')
            found_error = True

if not found_error:
    if samplesformat == 'manual':
        
        # load samples
        df_sample = pd.read_excel('input/%s/_OPTIONS_.xlsx' % input_folder, sheet_name='Manual Samples', header=None)
        
        # initialize list
        if datatype == 'gpl570':
            sample_files_ = []
        if datatype in ['rnaseq','gpl10558']:
            sample_columns_ = []
        if datatype == 'gpl10558':
            pvalue_columns_ = []
        sample_names_ = []

        # extract data
        for i in range(1,df_sample.shape[0]):
            if (isinstance(df_sample.loc[i][3],unicode)):
                
                # sample name
                sample_names_.append(str(df_sample.loc[i][3]))
                
                # sample file
                if datatype == 'gpl570':
                    if (isinstance(df_sample.loc[i][0],(unicode,float))):
                        sample_files_.append(str(df_sample.loc[i][0]))
                    else:
                        with open(fn_stderr,'a') as f:
                            f.write('ERROR - Options File - Manual Samples - Missing sample file for sample "%s"\n' % str(df_sample.loc[i][3]))
                        found_error = True
                        
                # sample column
                if datatype in ['rnaseq','gpl10558']:
                    if (isinstance(df_sample.loc[i][1],(unicode,float,int))):
                        sample_columns_.append(int(df_sample.loc[i][1]))
                    else:
                        with open(fn_stderr,'a') as f:
                            f.write('ERROR - Options File - Manual Samples - Missing sample column for sample "%s"\n' % str(df_sample.loc[i][3]))
                        found_error = True
                        
                # p-value column
                if datatype == 'gpl10558':
                    if (isinstance(df_sample.loc[i][2],(unicode,float,int))):
                        pvalue_columns_.append(int(df_sample.loc[i][2]))
                    else:
                        with open(fn_stderr,'a') as f:
                            f.write('ERROR - Options File - Manual Samples - Missing p-value column for sample "%s"\n' % str(df_sample.loc[i][3]))
                        found_error = True
            
        # error if no samples
        if len(sample_names_) == 0:
            with open(fn_stderr,'a') as f:
                f.write('ERROR - Options File - Manual Samples - Must have at least one sample\n')
            found_error = True

if not found_error:
    if samplesformat == 'manual':

        # get unique sample names
        sample_names = list(set(sample_names_))

        # sample file
        if datatype == 'gpl570':
            sample_files = []
            for sample in sample_names:
                sample_files.append([])
                for i in range(len(sample_names_)):
                    if sample_names_[i] == sample:
                        sample_files[-1].append(sample_files_[i])

        # sample column
        if datatype in ['rnaseq','gpl10558']:
            sample_columns = []
            for sample in sample_names:
                sample_columns.append([])
                for i in range(len(sample_names_)):
                    if sample_names_[i] == sample:
                        sample_columns[-1].append(sample_columns_[i])

        # p-value column
        if datatype == 'gpl10558':
            pvalue_columns = []
            for sample in sample_names:
                pvalue_columns.append([])
                for i in range(len(sample_names_)):
                    if sample_names_[i] == sample:
                        pvalue_columns[-1].append(pvalue_columns_[i])

if not found_error:

    # load data
    with open('_data_/geneids/geneids.pkl','r') as f:
        ensembltranscript_length, genesymbol_length, convert_ensembltranscript_genesymbol, convert_geneid_genesymbol, convert_ensemblgeneid_genesymbol, convert_refseqtranscript_ensembltranscript = pickle.load(f)

if not found_error:
    
    # spline function
    def piecewise_linear(x, x0, y0, k1, k2):
        return np.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])
    
    # gpl570
    if datatype == 'gpl570':
        with open('_data_/microarray/models/affy-gpl570.pkl','r') as f:
            p = pickle.load(f)
    
    # gpl570
    if datatype == 'gpl10558':
        with open('_data_/microarray/models/illumina-gpl10558.pkl','r') as f:
            p = pickle.load(f)

recon_genes = sorted(pd.read_table('_data_/recon/genes.tsv')['SYMBOL'].tolist())

if not found_error:

    # k_sp and k_dp parameters
    df_schwanhausser = pd.read_csv('_data_/schwanhausser/parameters.csv')
    schwanhausser_genes = df_schwanhausser['GENE'].tolist()
    schwanhausser_ksp = df_schwanhausser['KSP [1/hr]'].tolist()
    schwanhausser_kdp = df_schwanhausser['KDP [1/hr]'].tolist()

    # total cellular protein number
    with open('_data_/schwanhausser/protein_number.txt','r') as f:
        schwanhausser_protein_number = float(f.readlines()[0])

    # total cellular mRNA number
    with open('_data_/schwanhausser/mrna_number.txt','r') as f:
        schwanhausser_mrna_number = float(f.readlines()[0])
        
    # conversion factor for each recon gene
    conversion_factor = []
    for gene in recon_genes:
        if gene in schwanhausser_genes:
            conversion_factor.append(1./1000000*schwanhausser_mrna_number*schwanhausser_ksp[schwanhausser_genes.index(gene)]/schwanhausser_kdp[schwanhausser_genes.index(gene)]/schwanhausser_protein_number*1000000)
        else:
            conversion_factor.append(np.nan)

if not found_error:

    # load data
    with open('_data_/paxdb/paxdb.pkl','r') as f:
        recongenes_averageabundance = pickle.load(f)

if not found_error:
    if datatype == 'rnaseq':
    
        # file extension
        file_extension = options_fn.split('.')[-1]

        # load data
        if file_extension in ['xls','xlsx']:
            df_data = pd.read_excel('input/%s/%s' % (input_folder,options_fn), sheet_name=options_excel_sheet, skiprows=range(options_header-1), index_col=None)
        else:
            df_data = pd.read_table('input/%s/%s' % (input_folder,options_fn), delimiter=options_delimiter, skiprows=range(options_header-1), index_col=None)

        # extract gene ID's
        data_ids = [str(x) for x in df_data[df_data.columns[options_geneid_column-1]].tolist()] 
        
        # automatic sample names
        if samplesformat == 'automatic':
            sample_names = df_data.columns.tolist()[options_header:]
            sample_columns = [[x] for x in range(options_header+1,df_data.shape[1]+1)]
        
        # extract data
        data_gene = []
        for i in range(len(sample_names)):
            data_gene.append(df_data[df_data.columns[[x-1 for x in sample_columns[i]]]])

            # set index to gene ID's
            data_gene[-1].index = data_ids

            # sum values with same gene ID
            data_gene[-1] = data_gene[-1].groupby(data_gene[-1].index).sum()

        # new list of gene ID's
        data_ids = data_gene[-1].index.tolist()

if not found_error:
    if datatype == 'rnaseq':
        if options_geneid_type in ['Gene ID','Ensembl Gene ID']:
            
            # initialize gene symbols
            ids_query = []
            ids_symbols = []
            
            # gene ID
            if options_geneid_type == 'Gene ID':
            
                # iterate over data ids
                for i in range(len(data_ids)):
                    if int(data_ids[i]) in convert_geneid_genesymbol:
                        ids_query.append(data_ids[i])
                        ids_symbols.append(convert_geneid_genesymbol[int(data_ids[i])])
                    else:
                        with open(fn_stdout,'a') as f:
                            f.write('Gene ID Conversion - Gene ID "%s" not available and removed from dataset\n' % str(data_ids[i]))  
                        
            # ensembl gene ID
            elif options_geneid_type == 'Ensembl Gene ID':
            
                # iterate over data ids
                for i in range(len(data_ids)):
                    if str(data_ids[i]) in convert_ensemblgeneid_genesymbol:
                        ids_query.extend([data_ids[i] for x in convert_ensemblgeneid_genesymbol[str(data_ids[i])]])
                        ids_symbols.extend(convert_ensemblgeneid_genesymbol[str(data_ids[i])])
                    else:
                        with open(fn_stdout,'a') as f:
                            f.write('Gene ID Conversion - Ensembl Gene ID "%s" not available and removed from dataset\n' % str(data_ids[i]))
            
            # iterate over samples
            for i in range(len(sample_names)):
                
                # replace gene ID's with symbols
                data_gene[i] = data_gene[i].loc[ids_query]
                data_gene[i].index = ids_symbols
                
                # sum values with same gene symbol
                data_gene[i] = data_gene[i].groupby(data_gene[i].index).sum()
                
            # new list of gene ID's
            data_ids = data_gene[-1].index.tolist()

if not found_error:
    if datatype == 'rnaseq':
        if options_geneid_type == 'Gene Symbol':
            
            # symbols to keep
            data_ids = [x for x in data_ids if x in genesymbol_length]
            
            # iterate over samples
            for i in range(len(sample_names)):
                
                # replace gene ID's with symbols
                data_gene[i] = data_gene[i].loc[data_ids]

if not found_error:
    if datatype == 'rnaseq':
        if options_geneid_type in ['Gene Symbol','Gene ID','Ensembl Gene ID']:
            
            # counts
            if options_rnaseq == 'Count':
                
                # get length of largest transcript for that gene in Kb
                data_lengths = [genesymbol_length[gene]/1000. for gene in data_ids]
                
                # iterate over samples
                for i in range(len(sample_names)):

                    # divide gene values by largest length
                    data_gene[i] = data_gene[i].div(data_lengths, axis=0)
                    
            # RPKM/FPKM
            if options_rnaseq in ['Count','RPKM','FPKM','TPM']:
                
                # iterate over samples
                for i in range(len(sample_names)):

                    # convert to TPM: divide by sum of sample values, multiply by 1000000
                    data_gene[i] = data_gene[i].div(data_gene[i].sum(axis=0), axis=1) * 1000000

if not found_error:
    if datatype == 'rnaseq':
        if options_geneid_type == 'Ensembl Transcript ID':
            
            # symbols to keep
            data_ids = [x for x in data_ids if x in ensembltranscript_length]
            
            # iterate over samples
            for i in range(len(sample_names)):
                
                # replace gene ID's with symbols
                data_gene[i] = data_gene[i].loc[data_ids]

if not found_error:
    if datatype == 'rnaseq':
        if options_geneid_type == 'Ensembl Transcript ID':
            if options_rnaseq == 'Count':
                
                # get length of transcript in Kb
                data_lengths = [ensembltranscript_length[transcript]/1000. for transcript in data_ids]
                
                # iterate over samples
                for i in range(len(sample_names)):

                    # divide gene values by largest length
                    data_gene[i] = data_gene[i].div(data_lengths, axis=0)

if not found_error:
    if datatype == 'rnaseq':
        if options_geneid_type == 'Ensembl Transcript ID':
            
            # initialize gene symbols
            ids_query = []
            ids_symbols = []
           
            # ensembl transcript ID
            if options_geneid_type == 'Ensembl Transcript ID':
            
                # iterate over data ids
                for i in range(len(data_ids)):
                    if str(data_ids[i]) in convert_ensembltranscript_genesymbol:
                        ids_query.extend([data_ids[i] for x in convert_ensembltranscript_genesymbol[str(data_ids[i])]])
                        ids_symbols.extend(convert_ensembltranscript_genesymbol[str(data_ids[i])])
                    else:
                        with open(fn_stdout,'a') as f:
                            f.write('Gene ID Conversion - Ensembl Transcript ID "%s" not available and removed from dataset\n' % str(data_ids[i]))
            
            # iterate over samples
            for i in range(len(sample_names)):
                
                # replace transcript ID's with gene symbols
                data_gene[i] = data_gene[i].loc[ids_query]
                data_gene[i].index = ids_symbols
                
                # sum values with same gene symbol
                data_gene[i] = data_gene[i].groupby(data_gene[i].index).sum()
                
            # new list of gene ID's
            data_ids = data_gene[-1].index.tolist()

if not found_error:
    if datatype == 'rnaseq':
        if options_geneid_type == 'Ensembl Transcript ID':
                    
            # RPKM/FPKM
            if options_rnaseq in ['Count','RPKM','FPKM','TPM']:
                
                # iterate over samples
                for i in range(len(sample_names)):

                    # convert to TPM: divide by sum of sample values, multiply by 1000000
                    data_gene[i] = data_gene[i].div(data_gene[i].sum(axis=0), axis=1) * 1000000

if not found_error:
    if datatype == 'gpl570':

        # write file list
        os.mkdir('output/%s/_processing_' % output_folder)
        fn_input = 'output/%s/_processing_/input.txt' % output_folder
        fn_output = 'output/%s/_processing_/output.txt' % output_folder
        with open(fn_input,'w') as f_in:
            with open(fn_output,'w') as f_out:
                for i in range(len(sample_names)):
                    for fn in sample_files[i]:

                        # check .CEL file extension
                        if options_extension in fn:
                            f_in.write('input/%s/%s\n' % (input_folder, fn))
                            f_out.write('output/%s/_processing_/scan/%s\n' % (output_folder, fn))
                        else:
                            f_in.write('input/%s/%s%s\n' % (input_folder, fn, options_extension))
                            f_out.write('output/%s/_processing_/scan/%s%s\n' % (output_folder, fn, options_extension))
        
        # initialize scan folder
        os.mkdir('output/%s/_processing_/scan' % output_folder)
        
        # run scan on file list
        os.system('Rscript --vanilla _data_/microarray/scan/affy-gpl570.R %s %s' % (fn_input, fn_output))

if not found_error:
    if datatype == 'gpl570':

        # iterate over samples
        data_gene = []
        for i in range(len(sample_names)):
            data = []
            
            # load replicates
            for j in range(len(sample_files[i])):
                if options_extension in sample_files[i][j]:
                    data.append(pd.read_table('output/%s/_processing_/scan/%s' % (output_folder, sample_files[i][j])))
                else:
                    data.append(pd.read_table('output/%s/_processing_/scan/%s%s\n' % (output_folder, sample_files[i][j], options_extension)))

            # merge replicates
            data_gene.append(pd.concat(data, axis=1, sort=False))
            data_gene[-1].index = [x.split('_at')[0].split('.')[0] for x in data_gene[-1].index.tolist()]
            
        # list of gene ID's
        data_ids = data_gene[-1].index.tolist()

if not found_error:
    if datatype == 'gpl570':
        
        # iterate over samples
        for i in range(len(sample_names)):

            # subset
            data_gene[i] = data_gene[i].loc[[x for x in data_ids if x in ensembltranscript_length]]
            
        # list of gene ID's
        data_ids = data_gene[-1].index.tolist()

if not found_error:
    if datatype == 'gpl570':
        
        # iterate over samples
        for i in range(len(sample_names)):
            
            # iterate over replicates
            for column in data_gene[i].columns.tolist():
                
                # convert data
                data_gene[i].loc[:,column] = [10**x for x in piecewise_linear(data_gene[i][column].values.tolist(), *p)]

if not found_error:
    if datatype == 'gpl570':
            
            # initialize gene symbols
            ids_query = []
            ids_symbols = []
            
            # iterate over data ids
            for i in range(len(data_ids)):
                if str(data_ids[i]) in convert_ensembltranscript_genesymbol:
                    ids_query.extend([data_ids[i] for x in convert_ensembltranscript_genesymbol[str(data_ids[i])]])
                    ids_symbols.extend(convert_ensembltranscript_genesymbol[str(data_ids[i])])
                else:
                    with open(fn_stdout,'a') as f:
                        f.write('Gene ID Conversion - Ensembl Transcript ID "%s" not available and removed from dataset\n' % str(data_ids[i]))
            
            # iterate over samples
            for i in range(len(sample_names)):
                
                # replace transcript ID's with gene symbols
                data_gene[i] = data_gene[i].loc[ids_query]
                data_gene[i].index = ids_symbols
                
                # sum values with same gene symbol
                data_gene[i] = data_gene[i].groupby(data_gene[i].index).sum()
                
            # new list of gene ID's
            data_ids = data_gene[-1].index.tolist()

if not found_error:
    if datatype == 'gpl570':
        
        # iterate over samples
        for i in range(len(sample_names)):

            # convert to TPM: divide by sum of sample values, multiply by 1000000
            data_gene[i] = data_gene[i].div(data_gene[i].sum(axis=0), axis=1) * 1000000

if not found_error:
    
    # iterate over samples
    for i in range(len(sample_names)):
        
        # restrict to recon genes
        data_gene[i] = data_gene[i].reindex(recon_genes)

if not found_error:

    # iterate over samples
    for i in range(len(sample_names)):
        
        # replace 0's with nan's
        data_gene[i] = data_gene[i].replace(0, np.nan)
        
        # take average of replicates
        data_gene[i] = data_gene[i].mean(axis=1, skipna=True)

if not found_error:

    data_gene = pd.concat(data_gene, axis=1, sort=False)
    data_gene.columns = sample_names

if not found_error:

    # multiply TPM values by conversion factor to get PPM
    data_protein = data_gene.multiply(conversion_factor, axis=0)

if not found_error:
    
    # iterate over samples
    for i in range(len(sample_names)):
        
        # extract gene and protein expression
        protein_all = data_protein[sample_names[i]].tolist()
        gene_all = data_gene[sample_names[i]].tolist()
        keep_index = [a for a in range(len(protein_all)) if (not np.isnan(protein_all[a])) and (not np.isnan(gene_all[a]))]
        protein = [protein_all[a] for a in keep_index]
        gene = [gene_all[a] for a in keep_index]
        
        # linear regression of gene vs. protein expression
        reg = LinearRegression().fit([[x] for x in np.log10(gene)], np.log10(protein))
        
        # predict missing protein values with given gene values
        calculate_index = [a for a in range(len(protein_all)) if (np.isnan(protein_all[a])) and (not np.isnan(gene_all[a]))]
        gene_predict = data_gene.loc[[recon_genes[a] for a in calculate_index]][sample_names[i]].tolist()
        data_protein.loc[[recon_genes[a] for a in calculate_index],sample_names[i]] = [10**x for x in reg.predict([[x] for x in np.log10(gene_predict)])]

if not found_error:
    
    # iterate over samples
    for i in range(len(sample_names)):
        
        # iterate over missing values
        for gene in [x for x in data_protein.index if np.isnan(data_protein.loc[x][sample_names[i]])]:
            
            # if gene in PaxDB data
            if gene in recongenes_averageabundance:
            
                # fill in average PaxDB expression for that gene
                data_protein.loc[gene,sample_names[i]] = recongenes_averageabundance[gene]
            
            # otherwise. fill with average PaxDB expression for all genes
            else:
                data_protein.loc[gene,sample_names[i]] = recongenes_averageabundance['_ALL_']

if not found_error:
    
    # iterate over samples
    for i in range(len(sample_names)):
        
        # export data
        data_protein[sample_names[i]].to_csv('output/%s/%s.csv' % (output_folder, sample_names[i].replace('/','-')))