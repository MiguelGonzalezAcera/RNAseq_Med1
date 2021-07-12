import argparse
import logging
import os
import json
import subprocess
import pandas as pd


def GSVA(config, tool_name):
    """Get the
    """
    # Start log
    logging.info(f'Starting {tool_name} process')

    # Get input parameters
    norm_counts = config['tools_conf'][tool_name]['input']['norm_counts']
    design = config['tools_conf'][tool_name]['input']['design']

    # define the out directory
    out_dir = "/".join(config['tools_conf'][tool_name]['output']['heatmap'].split('/')[0:-1])

    # Define the rest of the output parameters
    heatmap = config['tools_conf'][tool_name]['output']['heatmap']
    heatmap_FC = config['tools_conf'][tool_name]['output']['heatmap'].replace(".png", "_FC.png")

    # Get options
    organism = config['options']['organism']
    project = config['project']
    samples = config['comparisons']

    # Create the command to run the pca R script
    command = ""

    # Make the outfolder if it doesnt wxist
    if not os.path.exists(out_dir):
        command += f"mkdir {out_dir};"

    # Loop through controls
    for control in samples:
        # Add the main command
        command += f'Rscript Rscripts/GSVA.r --heatmap {heatmap} --counts {norm_counts} --design {design} --organism {organism} --out_obj {out_dir}/{project}_GSVA.Rda --control {control} --comparisons {samples[control]}; '

    # Create the heatmap with the Fold Changes by limma
    res_list = []
    names_list = []
    for control in samples:
        for sample in samples[control].split(","):
            # Add result to list
            res_list.append(f"{out_dir}/{project}_GSVA_{sample}_{control}.Rda")
            names_list.append(f"{sample}_{control}")

    # Merge results in a single string
    res_str = ",".join(res_list)
    names_str = ",".join(names_list)

    # Add the command
    if len(names_list) > 1:
        command += f"Rscript Rscripts/GSVA_plot_FC.r --Rdata {res_str} --colnames {names_str} --heatmap {heatmap_FC} --organism {organism}; "

    # Run the command and log it
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


    logfile = config_dict["output"]["heatmap"].replace('.png','') + '_GSVA.log'
    logging.basicConfig(filename=logfile, level=logging.DEBUG, format='#[%(levelname)s]: - %(asctime)s - %(message)s')
    logging.info(f'Starting GSVA')

    config = {'tools_conf': {'GSVA': config_dict}}
    config['options'] = config['tools_conf']['GSVA']['options']

    GSVA(config, 'GSVA')

    logging.info(f'Finished GSVA')

if __name__ == "__main__":
    main()
