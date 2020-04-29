import argparse
import logging
import os
import glob
import json
import pandas as pd
import mysql.connector


def clustering_FC_heatmap(config, tool_name):
    """Get the
    """

    out_dir = "/".join(config['tools_conf'][tool_name]['output']['heatmap'].split('/')[0:-1])

    heatmap = config['tools_conf'][tool_name]['output']['heatmap']
    organism = config['options']['organism']

    # Create the command to run the pca R script
    command = ""
    if not os.path.exists(out_dir):
        command += f"mkdir {out_dir};"

    gene_file = config['tools_conf'][tool_name]['input']['genelist']

    if 'project' in config['tools_conf'][tool_name]['input']:
        project = config['tools_conf'][tool_name]['input']['project']
        command += f'Rscript /DATA/RNAseq_test/Scripts/Rscripts/clustering_FC.r --heatmap {heatmap} --project {project} --genelist {gene_file} --organism {organism}; '

        mydb = mysql.connector.connect(
          host="localhost",
          user="root",
          passwd="Plater1a",
          database="Projects"
        )

        mycursor = mydb.cursor()

        mycursor.execute(f"select Table_path from {project}")

        RData_lst = []
        for x in mycursor:
            RData_lst.append(x[0])
        RData_fix = ",".join(RData_lst)

    else:
        RData = config['tools_conf'][tool_name]['input']['RData']
        colnames = config['tools_conf'][tool_name]['input']['colnames']
        command += f'Rscript /DATA/RNAseq_test/Scripts/Rscripts/clustering_FC.r --heatmap {heatmap} --Rdata {RData} --colnames {colnames} --genelist {gene_file} --organism {organism}; '
        RData_fix = RData.replace('.Rda','.tsv')

    RData_wsp = RData_fix.replace(","," ")
    RData_fst = RData_fix.split(',')[0]
    RData_path = "/".join(RData_fst.split('/')[0:-1])
    RData_path_fix = RData_path.replace("/","\/") + "\/"
    outtab_DE = heatmap.replace('_heatmap.png','_DiffExpr.tsv')
    outtab_NC = heatmap.replace('_heatmap.png','_NormCounts.tsv')

    command += f'head -n +1 {RData_fst} | awk \'{{print \"Model\\t\" $0}}\' > {outtab_DE}; grep -f {gene_file} {RData_wsp} | sed "s/.tsv:/\t/g" | sed "s/{RData_path_fix}//g" >> {outtab_DE}; '
    command += f'head -n +1 {RData_path}/*_norm_counts.tsv | awk \'{{print \"EnsemblID\\t\" $0}}\' > {outtab_NC}; grep -f {gene_file} {RData_path}/*_norm_counts.tsv >> {outtab_NC}; '


    print(command)

    os.system(command)


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

    config = {'tools_conf': {'clustering_FC_heatmap': config_dict}}
    config['options'] = config['tools_conf']['clustering_FC_heatmap']['options']

    clustering_FC_heatmap(config, 'clustering_FC_heatmap')


if __name__ == "__main__":
    main()