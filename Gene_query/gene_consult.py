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

def barplot_counts(counts_df, barplot_counts_path, genename):
    names = []
    mean_lst = []
    std_lst = []

    for i in counts_df:
        names.append(i)
        mean_lst.append(counts_df[i]['mean'])
        std_lst.append(counts_df[i]['std'])

    N = len(mean_lst)
    ind = np.arange(N)

    fig, ax = plt.subplots(figsize=(10,10))

    ax.barh(ind, mean_lst, 0.35, xerr=std_lst, color = 'red')

    plt.yticks(ind, names)
    plt.xscale('log')

    ax.set_title(f'{genename} counts across models', fontsize=25)

    ax.set_xlabel('Normalized counts', fontsize=23)

    plt.tight_layout()
    plt.savefig(barplot_counts_path)

def barplot_FC(diff_df, barplot_FC_path, genename):
    names = []
    FC_lst = []
    colors = []

    for i in diff_df:
        names.append(i)
        FC_lst.append(diff_df[i]['foldchange'])
        if diff_df[i]['foldchange'] > 0:
            colors.append('red')
        else:
            colors.append('blue')

    N = len(FC_lst)
    ind = np.arange(N)

    fig, ax = plt.subplots(figsize=(10,10))

    ax.barh(ind, FC_lst, 0.35, color = colors)

    plt.yticks(ind, names)

    if max(FC_lst) > 1:
        ax.set_xlim((max(FC_lst)+1)*-1,(max(FC_lst)+1))
    elif min(FC_lst) < -1:
        ax.set_xlim((min(FC_lst)-1),(min(FC_lst)-1)*-1)
    else:
        ax.set_xlim(-1,1)

    ax.set_title(f'{genename} Fold change across models', fontsize=25)

    ax.set_xlabel('Fold change', fontsize=23)

    plt.tight_layout()
    plt.savefig(barplot_FC_path)

def mouse_models_plots(config, tool_name, mycursor, gene):
    mouse_models = ["Mouse_models_DSSdc_Cerldc", "Mouse_models_cDSSdc_Cerldc",
                "Mouse_models_OxCdc_Cerldc", "Mouse_models_RKOdc_Cerldc",
                "Mouse_models_TCdc_RKOdc", "Mouse_models_Janvdc_Cerldc",
                "Mouse_models_O12dc_KFdc", "Mouse_models_SPFdc_KFdc",
                "Mouse_models_O12D4_KFD4", "Mouse_models_SPFD4_KFD4"]

    diff_df = {}
    counts_df = {}
    tab_df = pd.DataFrame()
    model_list = []

    for model in mouse_models:
        header_command = f'DESCRIBE {model};'

        mycursor.execute(header_command)

        header=[]
        for row in mycursor:
            header.append(row[0])

        command = f"""select * from {model} where EnsGenes='{gene}';"""

        mycursor.execute(command)

        df_set = []
        for row in mycursor:
            df_set.append(row)

        df = pd.DataFrame(df_set)
        df.columns = header

        control = model.replace('Mouse_models_','').split('_')[1]
        if control == 'Cerldc':
            control = 'CErldc'
        sample = model.replace('Mouse_models_','').split('_')[0]

        diff_df[f'{sample}_{control}'] = {}
        diff_df[f'{sample}_{control}']['foldchange'] = df['log2FoldChange'][0]
        diff_df[f'{sample}_{control}']['pval'] = df['pvalue'][0]

        filter_col_control = [col for col in df if col.startswith(control)]
        filter_col_sample = [col for col in df if col.startswith(sample)]

        if control not in counts_df:
            counts_df[control] = {}
            counts_df[control]['mean'] = np.mean(df[filter_col_control].iloc[0].tolist())
            counts_df[control]['std'] = np.std(df[filter_col_control].iloc[0].tolist())
        if sample not in counts_df:
            counts_df[sample] = {}
            counts_df[sample]['mean'] = np.mean(df[filter_col_sample].iloc[0].tolist())
            counts_df[sample]['std'] = np.std(df[filter_col_sample].iloc[0].tolist())
        if tab_df.empty:
            tab_df = df[[col for col in df if not col.startswith(control) and not col.startswith(sample)]]
        else:
            tab_df = tab_df.append(df[[col for col in df if not col.startswith(control) and not col.startswith(sample)]])

        model_list.append(f'{sample}_{control}')

    tab_df['model'] = model_list

    FC_models_table_path = config['tools_conf'][tool_name]['output']['FC_models_table']
    tab_df.to_csv(FC_models_table_path, sep='\t', index=False)

    genename = config['genename']

    # 1.- Counts plot
    barplot_counts_path = config['tools_conf'][tool_name]['output']['barplot_counts']
    barplot_counts(counts_df, barplot_counts_path, genename)

    # 2.- Diff expression plot
    barplot_FC_path = config['tools_conf'][tool_name]['output']['barplot_FC']
    barplot_FC(diff_df, barplot_FC_path, genename)

def course_plot(counts_df, course_plot_path, genename):
    names = []
    mean_lst = []
    std_lst = []

    for i in counts_df:
        names.append(i)
        mean_lst.append(counts_df[i]['mean'])
        std_lst.append(counts_df[i]['std'])

    fig, ax = plt.subplots(figsize=(10,10))

    ax.errorbar(names, mean_lst, yerr=std_lst, marker="o", color='red', ecolor='black')

    ax.set_title(f'{genename} in a DSS time course', fontsize=25)

    ax.set_xlabel('Stage of inflammation', fontsize=23)
    ax.set_ylabel('Normalized counts', fontsize=23)

    plt.tight_layout()
    plt.savefig(course_plot_path)

def mouse_course_plots(config, tool_name, mycursor, gene):

    mouse_models_DSS = ["DSS_Time_course_Inf_mid_Healthy", "DSS_Time_course_Inf_hi_Healthy",
                        "DSS_Time_course_Rec_mod_Healthy", "DSS_Time_course_Rec_ful_Healthy"]

    counts_df = {}

    for model in mouse_models_DSS:
        header_command = f'DESCRIBE {model};'

        mycursor.execute(header_command)

        header=[]
        for row in mycursor:
            header.append(row[0])

        command = f"""select * from {model} where EnsGenes='{gene}';"""

        mycursor.execute(command)

        df_set = []
        for row in mycursor:
            df_set.append(row)

        df = pd.DataFrame(df_set)
        df.columns = header

        control = model.replace('DSS_Time_course_','').split('_')[-1]
        sample = '_'.join(model.replace('DSS_Time_course_','').split('_')[0:2])

        filter_col_control = [col for col in df if col.startswith(control)]
        filter_col_sample = [col for col in df if col.startswith(sample)]

        if control not in counts_df:
            counts_df[control] = {}
            counts_df[control]['mean'] = np.mean(df[filter_col_control].iloc[0].tolist())
            counts_df[control]['std'] = np.std(df[filter_col_control].iloc[0].tolist())
        if sample not in counts_df:
            counts_df[sample] = {}
            counts_df[sample]['mean'] = np.mean(df[filter_col_sample].iloc[0].tolist())
            counts_df[sample]['std'] = np.std(df[filter_col_sample].iloc[0].tolist())

    course_plot_path = config['tools_conf'][tool_name]['output']['course_plot']
    genename = config['genename']
    course_plot(counts_df, course_plot_path, genename)

def get_gene_id(genename):
    mouse_ref_df = pd.read_csv("/DATA/mouse_genes.tsv", sep='\t', index_col=None, header=None)
    mouse_ref_df.columns = ['ensembl','entrez','genename']

    gene = mouse_ref_df[mouse_ref_df['genename'] == genename]['ensembl'].tolist()[0]

    return gene

def gene_consult(config, tool_name):
    """"""

    genename = config['genename']

    outfolder = "/".join(config['tools_conf'][tool_name]['output']['FC_models_table'].split('/')[0:-1])

    if not os.path.exists(outfolder):
        os.system(f"mkdir {outfolder}")

    gene = get_gene_id(genename)

    plt.rc('xtick', labelsize=18)
    plt.rc('ytick', labelsize=21)

    mydb = mysql.connector.connect(
          host="localhost",
          user="root",
          passwd="Plater1a",
          database="RNAseq"
        )

    mycursor = mydb.cursor()

    mouse_models_plots(config, tool_name, mycursor, gene)

    mouse_course_plots(config, tool_name, mycursor, gene)

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
