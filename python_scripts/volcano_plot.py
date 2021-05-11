import argparse
import logging
import os
import glob
import json
import subprocess
import pandas as pd


def volcano_plot(config, tool_name):
    """Get the
    """

    logging.info(f'Starting {tool_name} process')

    organism = config['options']['organism']

    # Create the command to run the pca R script
    command = ""

    if 'volcanotouched' in config['tools_conf'][tool_name]['output']:
        out_dir_DE = "/".join(config['tools_conf'][tool_name]['input']['DEtouched'].split('/')[0:-1])
        out_dir = "/".join(config['tools_conf'][tool_name]['output']['volcanotouched'].split('/')[0:-1])
        volcanotouched = config['tools_conf'][tool_name]['output']['volcanotouched']
        samples = config['comparisons']


        path_list = []

        for control in samples:
            sample_ids = samples[control].split(",")
            for sample in sample_ids:
                id_tab = out_dir + "/" + f"{sample}_{control}_volcano.png"
                id_sample = out_dir_DE + "/" + config['project'] + "_" + f"{sample}_{control}.Rda"
                id_obj = f"{sample}_{control}"

                command += f'Rscript Rscripts/volcano_plot.r --out_plot {id_tab} --res {id_sample} --organism {organism}; '

                path_list.append([id_sample, id_tab])

        command += f'touch {volcanotouched}'

        # Add volcano plots paths to the design table
        path_df = pd.DataFrame(path_list, columns = ['Robj_path', 'Volcano_path'])
        design_df = pd.read_csv(config['tools_conf'][tool_name]['input']['design_tab'], sep='\t', index_col=None)
        design_df = pd.merge(design_df, path_df, on=['Robj_path'])
        design_df.to_csv(config['tools_conf'][tool_name]['output']['design_tab'], sep='\t', index=False)
    else:
        id_tab = config['tools_conf'][tool_name]['output']['volcano']
        id_sample = config['tools_conf'][tool_name]['input']['RData']
        out_dir = "/".join(config['tools_conf'][tool_name]['output']['volcano'].split('/')[0:-1])

        command += f'Rscript Rscripts/volcano_plot.r --out_plot {id_tab} --res {id_sample} --organism {organism}'

        if 'genelist' in config['tools_conf'][tool_name]['input']:
            genelist = config['tools_conf'][tool_name]['input']['genelist']

            command += f' --genelist {genelist}'

        command += '; '

    if not os.path.exists(out_dir):
        command = f"mkdir {out_dir};" + command

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

    logfile = config_dict["output"]["volcano"].replace('.png','') + '_volcano_plot.log'
    logging.basicConfig(filename=logfile, level=logging.DEBUG, format='#[%(levelname)s]: - %(asctime)s - %(message)s')
    logging.info(f'Starting volcano_plot')

    config = {'tools_conf': {'volcano_plot': config_dict}}
    config['options'] = config['tools_conf']['volcano_plot']['options']

    config['pipeline'] = 'volcano_plot'

    volcano_plot(config, 'volcano_plot')

    logging.info(f'Finished volcano_plot')

if __name__ == "__main__":
    main()
