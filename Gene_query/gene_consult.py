import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import argparse
import json
import logging
import pandas as pd
import mysql.connector
import numpy as np

def barplot_FC(diff_df, barplot_FC_path, genename):
    # Start lists for managing
    names = []
    FC_lst = []
    colors = []

    # Get the fold change and color list filled up
    for i in diff_df:
        names.append(i)
        FC_lst.append(diff_df[i]['foldchange'])
        if diff_df[i]['foldchange'] > 0:
            colors.append('red')
        else:
            colors.append('blue')

    # Establish the x positions for the barplot according to the length of the data
    N = len(FC_lst)
    ind = np.flip(np.arange(N))

    # Create the plot
    fig, ax = plt.subplots(figsize=(10,10))

    # Plot the data in an horizontal barplot
    ax.barh(ind, FC_lst, 0.4, color = colors)

    # Give the categorical ticks their names
    plt.yticks(ind, names)

    # Set x limits according to the amount of data using minimum values
    if max(FC_lst) > 1:
        ax.set_xlim((max(FC_lst)+1)*-1,(max(FC_lst)+1))
    elif min(FC_lst) < -1:
        ax.set_xlim((min(FC_lst)-1),(min(FC_lst)-1)*-1)
    else:
        ax.set_xlim(-1,1)

    # Set title and labels
    ax.set_title(f'{genename} Fold change across Diff. Expr.', fontsize=25)

    ax.set_xlabel('Fold change', fontsize=23)

    # Tighten layout and save
    plt.tight_layout()
    plt.savefig(barplot_FC_path)

def counts_plot(counts_df, counts_path, genename):
    # Start lists for managing
    names = []
    mean_lst = []
    std_lst = []

    # Get the fold change, mean and std list filled up
    for i in counts_df:
        names.append(i)
        mean_lst.append(counts_df[i]['mean'])
        std_lst.append(counts_df[i]['std'])

    # Establish the x positions for the barplot according to the length of the data
    N = len(mean_lst)
    ind = np.flip(np.arange(N))

    # Create the plot
    fig, ax = plt.subplots(figsize=(10,10))

    # Create the plot
    ax.barh(ind, mean_lst, 0.35, xerr=std_lst, color = 'red')

    # Give the categorical ticks their names
    plt.yticks(ind, names)
    plt.xscale('log')

    # Set title, labels and save
    ax.set_title(f'{genename} counts across samples', fontsize=25)

    ax.set_xlabel('Normalized counts', fontsize=23)

    plt.tight_layout()
    plt.savefig(counts_path)

def course_plot(counts_df, course_path, genename):
    # Start lists for managing
    names = []
    mean_lst = []
    std_lst = []

    # Get the fold change, mean and std list filled up
    for i in counts_df:
        names.append(i)
        mean_lst.append(counts_df[i]['mean'])
        std_lst.append(counts_df[i]['std'])

    # Create the plot
    fig, ax = plt.subplots(figsize=(10,10))

    # Include the line with the error bars
    ax.errorbar(names, mean_lst, yerr=std_lst, marker="o", color='red', ecolor='black')

    # Set title, labels and save
    ax.set_title(f'{genename} in a time course', fontsize=25)

    ax.set_xlabel('Stage of the time course', fontsize=23)
    ax.set_ylabel('Normalized counts', fontsize=23)

    plt.tight_layout()
    plt.savefig(course_path)

def get_gene_ids(genename):

    # Load the gene name file for mouse
    mouse_ref_df = pd.read_csv("/DATA/mouse_genes.tsv", sep='\t', index_col=None, header=None)
    mouse_ref_df.columns = ['ensembl','entrez','genename']

    # Turn genename column to uppercase
    mouse_ref_df['genename'] = mouse_ref_df['genename'].str.upper()

    # Transform name into ensembl
    mouse_gene = mouse_ref_df[mouse_ref_df['genename'] == genename.upper()]['ensembl'].tolist()[0]

    # Load the gene name file for human
    human_ref_df = pd.read_csv("/DATA/human_genes.tsv", sep='\t', index_col=None, header=None)
    human_ref_df.columns = ['ensembl','entrez','genename']

    # Turn genename column to uppercase
    human_ref_df['genename'] = human_ref_df['genename'].str.upper()

    # Transform name into ensembl
    human_gene = human_ref_df[human_ref_df['genename'] == genename.upper()]['ensembl'].tolist()[0]

    return mouse_gene, human_gene

def report_plots(config, tool_name, mycursor, genename, gene_display, comparisons, organism):
    # Define common columns
    DefaultColumns = ['id','EnsGenes','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj','Genes']

    for project in comparisons[organism]:
        # Start important objects for each iteration
        diff_df = {}
        counts_df = {}
        tab_df = []
        model_list = []

        # import the design the project
        design  = pd.read_csv(comparisons[organism][project]['design'], sep='\t')
        design.columns = ['sample_ID', 'Tr']

        # Get the models
        for control in comparisons[organism][project]['comparisons']:
            samples = comparisons[organism][project]['comparisons'][control].split(',')
            for sample in samples:
                # Obtain model name
                model = f'{project}_{sample}_{control}'

                # Obtain sample names for each treatment
                designSampleFilt = design[design['Tr'] == sample]['sample_ID'].tolist()
                designControlFilt = design[design['Tr'] == control]['sample_ID'].tolist()

                # Get header from database
                header_command = f'DESCRIBE {model};'

                mycursor.execute(header_command)

                header = []
                for row in mycursor:
                    header.append(row[0])

                # Get gene data from database in each model
                command = f"""select * from {model} where EnsGenes='{genename}';"""

                mycursor.execute(command)

                df_set = []
                for row in mycursor:
                    df_set.append(row)

                # Transform to dataframe w/ header
                df = pd.DataFrame(df_set)

                # Obtain the pertinent tables if df is not empty
                if not df.empty:
                    df.columns = header

                    # Add values of the fold change and p value to dict, for the barplot
                    diff_df[f'{sample}_{control}'] = {}
                    diff_df[f'{sample}_{control}']['foldchange'] = df['log2FoldChange'][0]
                    diff_df[f'{sample}_{control}']['pval'] = df['pvalue'][0]

                    # Select columns with the normalized counts
                    filter_col_control = df[designControlFilt]
                    filter_col_sample = df[designSampleFilt]

                    # Get mean and standard deviation of counts and store in the dict
                    if control not in counts_df:
                        counts_df[control] = {}
                        counts_df[control]['mean'] = np.nanmean(filter_col_control.iloc[0].tolist())
                        counts_df[control]['std'] = np.nanstd(filter_col_control.iloc[0].tolist())
                    if sample not in counts_df:
                        counts_df[sample] = {}
                        counts_df[sample]['mean'] = np.nanmean(filter_col_sample.iloc[0].tolist())
                        counts_df[sample]['std'] = np.nanstd(filter_col_sample.iloc[0].tolist())

                    # Add values to table for report table
                    tab_df.append(df[DefaultColumns].values.tolist()[0])

                # Make a batch of null data if the gene is not found
                else:
                    # Get 0 values into diff_df and counts_df
                    diff_df[f'{sample}_{control}'] = {}
                    diff_df[f'{sample}_{control}']['foldchange'] = 0
                    diff_df[f'{sample}_{control}']['pval'] = 0

                    if control not in counts_df:
                        counts_df[control] = {}
                        counts_df[control]['mean'] = 0
                        counts_df[control]['std'] = 0
                    if sample not in counts_df:
                        counts_df[sample] = {}
                        counts_df[sample]['mean'] = 0
                        counts_df[sample]['std'] = 0

                    # Add null content to table
                    tab_df.append(["-",gene,0,0,0,0,1,1,"-"])

                # Add model to list
                model_list.append(f'{sample}_{control}')

        # Set model names for Diff Expr table
        tab_df = pd.DataFrame(tab_df)
        tab_df.columns = ['id','EnsGenes','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj','Genes']

        # Add the model to the table
        tab_df['model'] = model_list

        # Save the table to a file
        FC_table_path = config['tools_conf'][tool_name]['output']['FC_table'].replace('Mouse_models', project)

        tab_df.to_csv(FC_table_path, sep='\t', index=False)

        # Retrieve the name of the gene again
        #genename =

        # Check if the project is a time course or a normal
        if comparisons[organism][project]['type'] == 'normal':
            # 1.- Counts plot
            counts_path = config['tools_conf'][tool_name]['output']['counts_plot'].replace('Mouse_models', project)
            counts_plot(counts_df, counts_path, gene_display)
        elif comparisons[organism][project]['type'] == 'timecourse':
            # 1.- Time course plot
            counts_path = config['tools_conf'][tool_name]['output']['counts_plot'].replace('Mouse_models', project)
            course_plot(counts_df, counts_path, gene_display)

        # 2.- Diff expression plot
        barplot_FC_path = config['tools_conf'][tool_name]['output']['FC_barplot'].replace('Mouse_models', project)
        barplot_FC(diff_df, barplot_FC_path, gene_display)

def gene_consult(config, tool_name):
    """"""

    # Extract gene name and route of the outfolder
    genename = config['genename']

    outfolder = "/".join(config['tools_conf'][tool_name]['output']['FC_table'].split('/')[0:-1])

    if not os.path.exists(outfolder):
        os.system(f"mkdir {outfolder}")

    # Transform gene name into ensembl id
    mouse_gene, human_gene = get_gene_ids(genename)

    # Establish size of the ticks for the plots
    plt.rc('xtick', labelsize=18)
    plt.rc('ytick', labelsize=21)

    # Create a connection and cursor to the database
    mydb = mysql.connector.connect(
          host="localhost",
          user="root",
          passwd="Plater1a",
          database="RNAseq"
        )
    mycursor = mydb.cursor()

    # Create a dict with the comparisons that have to be in the report (hand selected)
    comparisons = {
        "mouse": {
            "Mouse_models": {
                "design": "/VAULT/Thesis_proj/design.txt",
                "comparisons": {
                    'Cerldc': 'cDSSdc,DSSdc,OxCdc,RKOdc',
                    'RKOdc': 'TCdc'
                },
                "type": "normal"
            },
            'DSS_TimeCourse': {
                "design": "/VAULT/DSS_rec_evolution/design.txt",
                "comparisons": {
                    "Healthy": "Inf_mid,Inf_hi,Rec_mod,Rec_ful"
                },
                "type": "timecourse"
            },
            'WoundHealing': {
                "design": "/VAULT/20200629_Wound_Healing_TC/design.txt",
                "comparisons": {
                    "h0": "h6,h24,h48"
                },
            "type": "timecourse"
            }
        },
        "human": {
            "WashUCohort_EMTAB5783": {
                "design": "/VAULT/Human_data/E_MTAB_5783_WashU_Cohort/design.txt",
                "comparisons": {
                    "normal": "CD"
                },
                "type": "normal"
            },
            "RISK_GSE57945": {
                "design": "/VAULT/Human_data/GSE57945_IBD_RISK_Cohort_Ileum/design.txt",
                "comparisons": {
                    "Not_IBD_Male": "UC_Male,CD_M_MaInf_DUlcer,CD_M_MaInf_NDUlcer,CD_M_MiInf_DUlcer,CD_M_MiInf_NDUlcer,CD_M_NMiMaInf_NDUlcer",
                    "Not_IBD_Female": "UC_Female,CD_F_MaInf_DUlcer,CD_F_MaInf_NDUlcer,CD_F_MiInf_NDUlcer,CD_F_NMiMaInf_NDUlcer"
                },
                "type": "normal"
            },
            "PSC_EMTAB7915": {
                "design": "/VAULT/Human_data/E_MTAB_7915_PSC_cohort/design.txt",
                "comparisons": {
                    "normal": "sclerosing_cholangitis,ulcerative_colitis"
                },
                "type": "normal"
            },
            "RISK_GSE117993": {
                "design": "/VAULT/Human_data/GSE117993_IBD_RISK_cohort_Rectum/design.txt",
                "comparisons": {
                    "NotIBD": "cCD,iCD,UC"
                },
                "type": "normal"
            },
            "PROTECT_GSE109142": {
                "design": "/VAULT/Human_data/GSE109142_Ulcerative_colitis_PROTECT_Cohort/design.txt",
                "comparisons": {
                    "Control_Male": "UC_Male_5ASA,UC_Male_CSIV,UC_Male_CSOral",
                    "Control_Female": "UC_Female_5ASA,UC_Female_CSIV,UC_Female_CSOral"
                },
                "type": "normal"
            }
        }
    }

    # Get the mouse models data and plots
    report_plots(config, tool_name, mycursor, mouse_gene, genename, comparisons, "mouse")

    # Get the human models data and plots
    report_plots(config, tool_name, mycursor, human_gene, genename, comparisons, "human")

    # Close connection to database
    mycursor.close()
    mydb.close()

def get_arguments():
    """
    Function that parse arguments given by the user, returning a dictionary
    that contains all the values.
    """

    # Create the top-level parser
    parser = argparse.ArgumentParser()

    # Mandatory variables
    parser.add_argument('--config', required=True, help='Configuration file in json format')

    # Test and debug variables
    parser.add_argument('--dry_run', action='store_true', default=False, help='debug')
    parser.add_argument('--debug', '-d', action='store_true', default=False, help='dry_run')
    parser.add_argument('--test', '-t', action='store_true', default=False, help='test')

    # parse some argument lists
    args = parser.parse_args()

    return args

def main():
    """
    Main function of the script. Launches the rest of the process
    """

    # Get arguments from user input
    args = get_arguments()

    with open(args.config, 'r') as f:
        config_dict = json.load(f)

    logfile = config_dict["output"]["FC_models_table"].replace('.tsv','') + '_query.log'
    logging.basicConfig(filename=logfile, level=logging.DEBUG, format='#[%(levelname)s]: - %(asctime)s - %(message)s')
    logging.info(f'Starting gene_consult')

    config = {'tools_conf': {'gene_consult': config_dict}}
    config['options'] = config['tools_conf']['gene_consult']['options']

    gene_consult(config, 'gene_consult')

    logging.info(f'Finished gene_consult')

if __name__ == "__main__":
    main()
