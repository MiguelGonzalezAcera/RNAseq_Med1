import mysql.connector
import argparse
import logging
import os
import glob
import json
import subprocess
import pandas as pd

def query_database(genelist, tab_name, outfile):
    """
    """

    mydb = mysql.connector.connect(
      host="localhost",
      user="root",
      passwd="Plater1a",
      database="RNAseq"
    )

    mycursor = mydb.cursor()

    with open(genelist) as f:
        genes = f.read().splitlines()
    genes = "\',\'".join(genes)

    header_command = f'DESCRIBE RNAseq.{tab_name};'
    logging.info(header_command)
    mycursor.execute(header_command)

    header=[]
    for row in mycursor:
        header.append(row[0])

    command = f"""select * from RNAseq.{tab_name} where EnsGenes in ('{genes}');"""
    logging.info(command)
    mycursor.execute(command)

    df_set = []
    for row in mycursor:
        df_set.append(row)

    df = pd.DataFrame(df_set)
    df.columns = header

    df.to_csv(outfile, sep='\t', index=False)

def fix_genelists(genelists_path,  organism):
    """"""
    genelist_files = genelists_path.split(',')

    genelist_path = "/".join(genelist_files[0].split('/')[0:-1])

    all_genes = pd.read_csv(f"/DATA/{organism}_genes.tsv", sep='\t', index_col=None, header=None)
    all_genes.columns = ['ensembl','entrez','genename']

    resdf = pd.DataFrame()

    for file in genelist_files:
        with open(file, 'r') as filehandle:
            glist = [i.rstrip() for i in filehandle.readlines()]

        entrez_genes = all_genes[all_genes['ensembl'].isin(glist)]['entrez'].tolist()

        label = [file.split('/')[-1].replace('_ensembl.txt','')]*len(entrez_genes)

        tmp_df = pd.DataFrame([label,entrez_genes]).transpose()

        if resdf.empty:
            resdf = tmp_df
        else:
            resdf = pd.concat([resdf,tmp_df]).dropna()

    resdf.columns = ['group','entrez']
    resdf['entrez'] = resdf['entrez'].astype('int')

    genegroups_path = f"{genelist_path}/genegroups_test.txt"

    resdf.to_csv(genegroups_path, sep='\t', index=False, header=False)

    return genegroups_path


def GSEA(config, tool_name):
    """Get the counts of a number of bam files in a directory
    """
    logging.info(f'Starting {tool_name} process')

    in_obj = config['tools_conf'][tool_name]['input']['in_obj']
    in_obj_tab = in_obj.replace('.Rda','.tsv')
    in_obj_path = "/".join(in_obj.split('/')[0:-1])

    organism = config['options']['organism']

    gseaplot = config['tools_conf'][tool_name]['output']['out_plot']
    id_tab_DExpr = gseaplot.replace(".png","_DExpr.tsv")
    id_tab_Ncounts = gseaplot.replace(".png","_norm_counts.tsv")
    out_dir = "/".join(gseaplot.split('/')[0:-1])

    genegroup_path = config['tools_conf'][tool_name]['input']['genegroup']
    genegroup_path_spac = genegroup_path.split(",")
    genegroup = fix_genelists(genegroup_path, organism)

    # Create the command to run the pca R script
    command = ""

    if not os.path.exists(out_dir):
        command += f"mkdir {out_dir};"

    command += f'Rscript Rscripts/GSEA.r --genegroup {genegroup} --in_obj {in_obj} --gseaplot {gseaplot} --organism {organism}; '

    for glist in genegroup_path_spac:
        glist_name = glist.split('/')[-1].replace(".txt",'')
        id_tab_DExpr_gl = id_tab_DExpr.replace('.tsv',f'_{glist_name}.tsv')
        in_obj_name = in_obj.split('/')[-1].replace('.Rda','')
        query_database(glist,in_obj_name,id_tab_DExpr_gl)

    logging.info(f'Running command: {command}')
    for cmd in command.split('; '):
        output = subprocess.run(cmd, shell=True, executable='/bin/bash', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stout = output.stdout.decode('utf-8')
        error = output.stderr.decode('utf-8')
        if output.returncode == 1:
            logging.error(f'{cmd}\n\n{error}')
        elif output.returncode == 0:
            logging.info(stout)
            logging.info(error)

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

    logfile = config_dict["output"]["out_plot"].replace('.png','') + '_GSEA.log'
    logging.basicConfig(filename=logfile, level=logging.DEBUG, format='#[%(levelname)s]: - %(asctime)s - %(message)s')
    logging.info(f'Starting GSEA')

    config = {'tools_conf': {'GSEA': config_dict}}
    config['options'] = config['tools_conf']['GSEA']['options']

    GSEA(config, 'GSEA')

    logging.info(f'Finished GSEA')

if __name__ == "__main__":
    main()
