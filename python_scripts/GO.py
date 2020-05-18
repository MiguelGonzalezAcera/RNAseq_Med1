import argparse
import logging
import os
import subprocess
import glob
import json

def GO_enrichment(config, tool_name):
    """Get the counts of a number of bam files in a directory
    """

    logging.info(f'Starting {tool_name} process')

    organism = config['options']['organism']

    # Create the command to run the pca R script
    command = ""
    if 'gotouched' in config['tools_conf'][tool_name]['output']:
        out_dir_DE = "/".join(config['tools_conf'][tool_name]['input']['DEtouched'].split('/')[0:-1])

        out_dir = "/".join(config['tools_conf'][tool_name]['output']['gotouched'].split('/')[0:-1])
        GOtouched = config['tools_conf'][tool_name]['output']['gotouched']
        samples = config['comparisons']

        if not os.path.exists(out_dir):
            command += f"mkdir {out_dir};"

        for control in samples:
            sample_ids = samples[control].split(",")
            for sample in sample_ids:
                id_dir = out_dir + "/" + f"{sample}_{control}"
                id_tab = out_dir + "/" + f"{sample}_{control}" + "/" + f"{sample}_{control}_GO.tsv"
                id_tab_BP = out_dir + "/" + f"{sample}_{control}" + "/" + f"{sample}_{control}_GO_BP.Rda"
                id_tab_MF = out_dir + "/" + f"{sample}_{control}" + "/" + f"{sample}_{control}_GO_MF.Rda"
                id_tab_CC = out_dir + "/" + f"{sample}_{control}" + "/" + f"{sample}_{control}_GO_CC.Rda"
                id_geneids = out_dir + "/" + f"{sample}_{control}" + "/" + f"{sample}_{control}_GO_entrezgeneids.Rda"
                id_sample = out_dir_DE + "/" + config['project'] + "_" + f"{sample}_{control}.Rda"
                id_universe = out_dir_DE + "/" + config['project'] + "_universe.Rda"

                command += f"mkdir {id_dir}; "
                command += f'Rscript /DATA/RNAseq_test/Scripts/Rscripts/GO_enrichment.r --out_tab {id_tab} --obj {id_sample} --universe {id_universe} --organism {organism}; '
                command += f'Rscript /DATA/RNAseq_test/Scripts/Rscripts/GO_enrichment_plots.r --out_tab {id_tab_BP} --organism {organism} --geneids {id_geneids}; '
                command += f'Rscript /DATA/RNAseq_test/Scripts/Rscripts/GO_enrichment_plots.r --out_tab {id_tab_MF} --organism {organism} --geneids {id_geneids}; '
                command += f'Rscript /DATA/RNAseq_test/Scripts/Rscripts/GO_enrichment_plots.r --out_tab {id_tab_CC} --organism {organism} --geneids {id_geneids}; '
        command += f'touch {GOtouched}'
    else:
        out_tab = config['tools_conf'][tool_name]['output']['out_tab']
        out_dir = "/".join(config['tools_conf'][tool_name]['output']['out_tab'].split('/')[0:-1])
        id_tab_BP = out_tab.replace(".tsv","_BP.Rda")
        id_tab_MF = out_tab.replace(".tsv","_MF.Rda")
        id_tab_CC = out_tab.replace(".tsv","_CC.Rda")
        id_geneids = out_tab.replace(".tsv","_entrezgeneids.Rda")
        id_tab_DExpr = out_tab.replace(".tsv","_DExpr.tsv")
        id_tab_Ncounts = out_tab.replace(".tsv","_norm_counts.tsv")

        in_obj = config['tools_conf'][tool_name]['input']['in_obj']
        in_obj_tab = in_obj.replace('.Rda','.tsv')
        in_obj_path = "/".join(config['tools_conf'][tool_name]['input']['in_obj'].split('/')[0:-1])
        universe = config['tools_conf'][tool_name]['input']['universe']
        genelist = config['tools_conf'][tool_name]['input']['genelist']

        if not os.path.exists(out_dir):
            command += f"mkdir {out_dir};"

        command += f'Rscript /DATA/RNAseq_test/Scripts/Rscripts/GO_enrichment.r --out_tab {id_tab} --obj {in_obj} --universe {universe} --organism {organism} --genelist {genelist}; '
        command += f'Rscript /DATA/RNAseq_test/Scripts/Rscripts/GO_enrichment_plots.r --out_tab {id_tab_BP} --organism {organism} --geneids {id_geneids}; '
        command += f'Rscript /DATA/RNAseq_test/Scripts/Rscripts/GO_enrichment_plots.r --out_tab {id_tab_MF} --organism {organism} --geneids {id_geneids}; '
        command += f'Rscript /DATA/RNAseq_test/Scripts/Rscripts/GO_enrichment_plots.r --out_tab {id_tab_CC} --organism {organism} --geneids {id_geneids}; '
        command += f"head -n +1 {in_obj_tab} > {id_tab_DExpr}; grep -f {genelist} {in_obj_tab} >> {id_tab_DExpr}; "
        command += f"head -n +1 {in_obj_path}/*_norm_counts.tsv | awk \'{{print \"Ensembl\\t\" $0}}\' > {id_tab_Ncounts}; grep -f {genelist} {in_obj_path}/*_norm_counts.tsv >> {id_tab_Ncounts}; "

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
    parser.add_argument('--counts', required=True, help='Table with the counts of the assay, straight from featurecounts (so far)')
    parser.add_argument('--design', required=True, help='Table with the design of the experiment')
    parser.add_argument('--out_dir', required=True, help='Directory for all of the plots)')

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

    logfile = config_dict["output"]["out_tab"].replace('.tsv','') + '_GO_enrichment.log'
    logging.basicConfig(filename=logfile, level=logging.DEBUG, format='#[%(levelname)s]: - %(asctime)s - %(message)s')
    logging.info(f'Starting GO_enrichment')

    config = {'tools_conf': {'GO_enrichment': config_dict}}
    config['options'] = config['tools_conf']['GO_enrichment']['options']

    GO_enrichment(config, 'GO_enrichment')

    logging.info(f'Finished GO_enrichment')

if __name__ == "__main__":
    main()
