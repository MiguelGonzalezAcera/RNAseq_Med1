import argparse
import logging
import os
import glob
import pandas as pd
import python_scripts.python_functions as pf


def clustering_heatmap(config, tool_name):
    """Get the
    """

    out_dir = "/".join(config['tools_conf'][tool_name]['output']['heatmap'].split('/')[0:-1])

    norm_counts = config['tools_conf'][tool_name]['input']['norm_counts']
    heatmap = config['tools_conf'][tool_name]['output']['heatmap']
    organism = config['options']['organism']

    # Create the command to run the pca R script
    command = ""
    if not os.path.exists(out_dir):
        command += f"mkdir {out_dir};"

    genes_list = []

    if 'project' in config['tools_conf'][tool_name]['input']:
        project = config['tools_conf'][tool_name]['input']['project']
        command += f'Rscript /DATA/RNAseq_test/Scripts/Rscripts/clustering_FC.r --heatmap {heatmap} --project --genelist {gene_file} --organism {organism}; '
    else:
        RData = config['tools_conf'][tool_name]['input']['RData']
        colnames = config['tools_conf'][tool_name]['input']['colnames']
        command += f'Rscript /DATA/RNAseq_test/Scripts/Rscripts/clustering_FC.r --heatmap {heatmap} --RData {RData} --colnames {colnames} --genelist {gene_file} --organism {organism}; '

    print(command)

    pf.run_command(command)


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

    pca(config, 'clustering_FC_heatmap')


if __name__ == "__main__":
    main()
