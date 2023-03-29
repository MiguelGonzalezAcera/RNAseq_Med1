import argparse
import logging
import os
import glob
import subprocess
import json
import pandas as pd
import python_scripts.python_functions as pf


def KEGG_enrichment(config, tool_name):
    """Get the counts of a number of bam files in a directory
    """

    logging.info(f'Starting {tool_name} process')

    # Get the organism
    organism = config['options']['organism']

    # Create the command to run the pca R script
    command = ""

    # Choose to run inside the pipeline or outside
    if 'keggtouched' in config['tools_conf'][tool_name]['output']:
        # Get directory of the DE files
        out_dir_DE = "/".join(config['tools_conf'][tool_name]['input']['DEtouched'].split('/')[0:-1])

        # Define the directory where the results will be stored
        out_dir = "/".join(config['tools_conf'][tool_name]['output']['keggtouched'].split('/')[0:-1])

        # Get the marker for end
        keggtouched = config['tools_conf'][tool_name]['output']['keggtouched']

        # Get the comparisons to use
        samples = config['comparisons']

        # Make the out dir if it doesnt exist
        if not os.path.exists(out_dir):
            command += f"mkdir {out_dir};"

        # Iter through the samples to make the analysis
        for control in samples:
            sample_ids = samples[control].split(",")
            for sample in sample_ids:
                # Get an identifyer for the object
                id_obj = f"{sample}_{control}"

                # define a separate folder for each analysis
                id_dir = out_dir + "/" + id_obj

                # Define the name of the out table
                id_tab = id_dir + "/" + f"{id_obj}_KEGG.tsv"

                # Obtain the path for the R object containing the DE result
                id_sample = out_dir_DE + "/" + config['project'] + "_" + f"{id_obj}.Rda"

                # Make the commands
                command += f"mkdir {id_dir}; "
                command += f'Rscript Rscripts/KEGG_enrichment.r --out_tab {id_tab} --in_obj {id_sample} --id {id_obj} --organism {organism}; '
        
        # make the command for the ending control file
        command += f'touch {keggtouched}'

    else:
        # Get the name of the out file
        out_tab = config['tools_conf'][tool_name]['output']['out_tab']

        # Get the out directory in case the folder does not exist
        out_dir = "/".join(config['tools_conf'][tool_name]['output']['out_tab'].split('/')[0:-1])

        # The RData object from the DE analysis
        in_obj = config['tools_conf'][tool_name]['input']['in_obj']

        # The name you want to give the analysis
        id = config['tools_conf'][tool_name]['input']['id']

        # The list of genes file, if any
        genelist = config['tools_conf'][tool_name]['input']['genelist']

        if not os.path.exists(out_dir):
            command += f"mkdir {out_dir};"

        command += f'Rscript Rscripts/KEGG_enrichment.r --out_tab {out_tab} --in_obj {in_obj} --id {id} --organism {organism} --genelist {genelist}; '

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

    logfile = config_dict["output"]["out_tab"].replace('.tsv','') + '_KEEG_enrichment.log'
    logging.basicConfig(filename=logfile, level=logging.DEBUG, format='#[%(levelname)s]: - %(asctime)s - %(message)s')
    logging.info(f'Starting KEEG_enrichment')

    config = {'tools_conf': {'KEEG_enrichment': config_dict}}
    config['options'] = config['tools_conf']['KEEG_enrichment']['options']

    KEGG_enrichment(config, 'KEEG_enrichment')

    logging.info(f'Finished KEEG_enrichment')

if __name__ == "__main__":
    main()
