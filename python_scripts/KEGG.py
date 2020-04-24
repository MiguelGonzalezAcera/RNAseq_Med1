import argparse
import logging
import os
import glob
import json


def KEGG_enrichment(config, tool_name):
    """Get the counts of a number of bam files in a directory
    """

    organism = config['options']['organism']

    # Create the command to run the pca R script
    command = ""

    if 'keggtouched' in config['tools_conf'][tool_name]['output']:
        out_dir_DE = "/".join(config['tools_conf'][tool_name]['input']['DEtouched'].split('/')[0:-1])

        out_dir = "/".join(config['tools_conf'][tool_name]['output']['keggtouched'].split('/')[0:-1])
        keggtouched = config['tools_conf'][tool_name]['output']['keggtouched']
        samples = config['comparisons']

        if not os.path.exists(out_dir):
            command += f"mkdir {out_dir};"

        for control in samples:
            sample_ids = samples[control].split(",")
            for sample in sample_ids:
                id_dir = out_dir + "/" + f"{sample}_{control}"
                id_tab = out_dir + "/" + f"{sample}_{control}" + "/" + f"{sample}_{control}_KEGG.tsv"
                id_sample = out_dir_DE + "/" + config['project'] + "_" + f"{sample}_{control}.Rda"
                id_obj = f"{sample}_{control}"

                command += f"mkdir {id_dir}; "
                command += f'Rscript /DATA/RNAseq_test/Scripts/Rscripts/KEGG_enrichment.r --out_tab {id_tab} --in_obj {id_sample} --id {id_obj} --organism {organism}; '
        command += f'touch {keggtouched}'
    else:
        out_tab = config['tools_conf'][tool_name]['output']['out_tab']
        out_dir = config['tools_conf'][tool_name]['output']['out_tab'].split('/')[0:-1]

        in_obj = config['tools_conf'][tool_name]['input']['in_obj']
        id = config['tools_conf'][tool_name]['input']['id']
        genelist = config['tools_conf'][tool_name]['input']['genelist']

        if not os.path.exists(out_dir):
            command += f"mkdir {out_dir};"

        command += f'Rscript /DATA/RNAseq_test/Scripts/Rscripts/KEGG_enrichment.r --out_tab {out_tab} --in_obj {in_obj} --id {id} --organism {organism} --genelist {genelist}'

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

    config = {'tools_conf': {'KEEG_enrichment': config_dict}}
    config['options'] = config['tools_conf']['KEEG_enrichment']['options']

    KEEG_enrichment(config, 'KEEG_enrichment')


if __name__ == "__main__":
    main()
