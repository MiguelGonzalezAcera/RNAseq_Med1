import argparse
import logging
import os
import glob
import json
import subprocess
import pandas as pd
import mysql.connector

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


def clustering_FC_heatmap(config, tool_name):
    """Get the
    """
    logging.info(f'Starting {tool_name} process')

    organism = config['options']['organism']

    genes_groups_reference = {
        "Goblet": {
            'human': "/DATA/references/Gene_markers/Goblet_human_markers.txt",
            'mouse': "/DATA/references/Gene_markers/Goblet_markers.txt"
            },
        "Paneth": {
            'human': "/DATA/references/Gene_markers/Paneth_human_markers.txt",
            'mouse': "/DATA/references/Gene_markers/Paneth_markers.txt"
            },
        "Enteroendocrine": {
            'human': "/DATA/references/Gene_markers/Enteroendocrine_human_markers.txt",
            'mouse': "/DATA/references/Gene_markers/Enteroendocrine_markers.txt"
            },
        "Enterocyte_mat_distal": {
            'human': "/DATA/references/Gene_markers/Enterocyte_mat_distal_human_markers.txt",
            'mouse': "/DATA/references/Gene_markers/Enterocyte_mat_distal_markers.txt"
            },
        "Enterocyte_mat_proximal": {
            'human': "/DATA/references/Gene_markers/Enterocyte_mat_proximal_human_markers.txt",
            'mouse': "/DATA/references/Gene_markers/Enterocyte_mat_proximal_markers.txt"
            },
        "Tuft": {
            'human': "/DATA/references/Gene_markers/Tuft_human_markers.txt",
            'mouse': "/DATA/references/Gene_markers/Tuft_markers.txt"
            }
    }

    if config['pipeline'] == 'RNAseq' or config['pipeline'] == 'RNAseq_update':
        genes_dict = genes_groups_reference

        config['tools_conf'][tool_name]['input']['project'] = config['project']
        heatmaptouched = config['tools_conf'][tool_name]['output']['heatmaptouched']
        heatmap = "/".join(heatmaptouched.split('/')[0:-1]) + '/CellType_heatmap.png'
    else:
        genes_dict = {"genelist": {organism: config['tools_conf'][tool_name]['input']['genelist']}}
        heatmaptouched = '/DATA/tmp/ht.txt'
        heatmap = config['tools_conf'][tool_name]['output']['heatmap']

    out_dir = "/".join(heatmap.split('/')[0:-1])

    # Create the command to run the pca R script
    command = ""
    if not os.path.exists(out_dir):
        command += f"mkdir {out_dir};"


    mydb = mysql.connector.connect(
      host="localhost",
      user="root",
      passwd="Plater1a",
      database="Projects"
    )

    mycursor = mydb.cursor()

    for pway in genes_dict:
        heatmap_pway = heatmap.replace('.png',f'_{pway}.png')

        if 'project' in config['tools_conf'][tool_name]['input']:
            project = config['tools_conf'][tool_name]['input']['project']

            command += f'Rscript /DATA/RNAseq_test/Scripts/Rscripts/clustering_FC.r --heatmap {heatmap_pway} --project {project} --genelist {genes_dict[pway][organism]} --organism {organism}; '

            mycursor.execute(f"select Robj_path from {project}")

            RData_lst = []
            for x in mycursor:
                RData_lst.append(x[0])
            RData_fix = ",".join(RData_lst)

        else:
            RData = config['tools_conf'][tool_name]['input']['RData']
            colnames = config['tools_conf'][tool_name]['input']['colnames']
            command += f'Rscript /DATA/RNAseq_test/Scripts/Rscripts/clustering_FC.r --heatmap {heatmap_pway} --Rdata {RData} --colnames {colnames} --genelist {genes_dict[pway][organism]} --organism {organism}; '
            RData_fix = RData.replace('.Rda','.tsv')

        for RData_file in RData_fix.split(','):
            RData_in = RData_file.split('/')[-1].replace('.tsv','').replace('.Rda','')
            RData_out = heatmap.replace('heatmap.png',f'{RData_in}_{pway}_DE.tsv')
            query_database(genes_dict[pway][organism], RData_in, RData_out)

    mycursor.close()
    mydb.close()

    command += f'touch {heatmaptouched}; '

    logging.info(f'Running command: {command}')
    for cmd in command.split('; '):
        output = subprocess.run(cmd, shell=True, executable='/bin/bash', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stout = output.stdout.decode('utf-8')
        error = output.stderr.decode('utf-8')
        if output.returncode == 1:
            logging.error(f'{cmd}\n\n{error}')
            raise OSError(f'Error in command: {cmd}\n\n{error}')
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

    logfile = config_dict["output"]["heatmap"].replace('.png','') + '_clustering_FC.log'
    logging.basicConfig(filename=logfile, level=logging.DEBUG, format='#[%(levelname)s]: - %(asctime)s - %(message)s')
    logging.info(f'Starting clustering_FC')

    config = {'tools_conf': {'clustering_FC_heatmap': config_dict}}
    config['options'] = config['tools_conf']['clustering_FC_heatmap']['options']

    config['pipeline'] = 'Clustering_fc'

    clustering_FC_heatmap(config, 'clustering_FC_heatmap')

    logging.info(f'Finished clustering_FC')

if __name__ == "__main__":
    main()
