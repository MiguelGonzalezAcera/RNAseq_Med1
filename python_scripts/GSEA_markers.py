import mysql.connector
import argparse
import logging
import os
import glob
import json
import subprocess
import pandas as pd

def fix_genelists(genelist, outpath, organism):
    """"""
    # Get gene reference table
    all_genes = pd.read_csv(f"/DATA/{organism}_genes.tsv", sep='\t', index_col=None, header=None)
    all_genes.columns = ['ensembl','entrez','genename']

    # Open the list of genes
    with open(genelist, 'r') as filehandle:
        glist = [i.rstrip() for i in filehandle.readlines()]

    # Select lines with the genes from the list
    entrez_genes = all_genes[all_genes['ensembl'].isin(glist)]['entrez'].tolist()

    # Use the label from the genelist
    label_id = genelist.split('/')[-1].replace('_ensembl.txt','')
    label = [label_id]*len(entrez_genes)

    # Transfrom into dataframe
    resdf = pd.DataFrame([label,entrez_genes]).transpose()

    # Change column names and data types
    resdf.columns = ['group','entrez']

    # Remove rows with Nan
    resdf = resdf.dropna()
    resdf['entrez'] = resdf['entrez'].astype('int')

    # Dump onto file
    genegroups_path = f"{outpath}/{label_id}_genegroups_test.txt"
    resdf.to_csv(genegroups_path, sep='\t', index=False, header=False)

    return genegroups_path

def GSEA_markers_plots(in_obj, outpath, organism, command):
    """
    Do plots for markers
    """
    # Define markers
    gene_markers = {
        "mouse": {
            "EnterocyteDist": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/EnterocyteDist_ensembl.txt",
            "EnterocyteProx": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/EnterocyteProx_ensembl.txt",
            "Enteroendocrine": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/Enteroendocrine_ensembl.txt",
            "Goblet": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/Goblet_ensembl.txt",
            "M_cells": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/M_cells_ensembl.txt",
            "Paneth": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/Paneth_ensembl.txt",
            "StemProg": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/StemProg_ensembl.txt",
            "TAProg": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/TAProg_ensembl.txt",
            "Tuft": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/Tuft_ensembl.txt",
            "Fibroblasts": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/Fibroblasts_ensembl.txt",
            "MO_DC": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/MO_DC_ensembl.txt",
            "Plasma_cells": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/Plasma_cells_ensembl.txt",
            "T_cells": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/T_cells_ensembl.txt",
            "B_cells": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/B_cells_ensembl.txt",
            "Mast_cells": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/Mast_cells_ensembl.txt",
            "NK_ILC1": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/NK_ILC1_ensembl.txt",
            "Endothelial": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/Endothelial_ensembl.txt",
            "Neutrophils": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/Neutrophils_ensembl.txt",
            "Smooth_muscle": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/Smooth_muscle_ensembl.txt",
            "EntericGlial": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/EntericGlial_ensembl.txt",
            "EntericNeuron": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/EntericNeuron_ensembl.txt"
        }
    }

    # Loop through the markers
    for marker in gene_markers['mouse']:
        marker_path = gene_markers['mouse'][marker]

        # Make the marker outfile
        gseaplot_mark = outpath + "/" + in_obj.split('/')[-1].replace(".Rda",f"_{marker}_GSEA.png")

        # Fix the markers into a gene group
        genegroup = fix_genelists(marker_path, outpath, organism)

        # Add command
        command += f'Rscript Rscripts/GSEA.r --genegroup {genegroup} --in_obj {in_obj} --gseaplot {gseaplot_mark} --organism {organism}; '

    return(command)


def GSEA_markers(config, tool_name):
    """Get the counts of a number of bam files in a directory
    """
    logging.info(f'Starting {tool_name} process')

    # Extract the organism
    organism = config['options']['organism']

    # Extract the comparisons
    comparisons = config['comparisons']

    # Extract the infiles and the project name
    in_path = "/".join(config['tools_conf'][tool_name]['input']['DEtouched'].split('/')[0:-1])
    project = config['project']

    # Extract the out path and out markerfile
    out_path = "/".join(config['tools_conf'][tool_name]['output']['GSEAMtouched'].split('/')[0:-1])
    outmarker = config['tools_conf'][tool_name]['output']['GSEAMtouched']

    # Create the command to run the pca R script
    command = ""

    # Create outfolder if it doesnt exist
    if not os.path.exists(out_path):
        command += f"mkdir {out_path};"

    # Loop through the samples and controls
    for control in comparisons:
        samples = comparisons[control].split(",")
        for sample in samples:
            # Read the table with the whole data
            in_obj = f"{in_path}/{project}_{sample}_{control}.Rda"

            command = GSEA_markers_plots(in_obj, out_path, organism, command)

    command += f"touch {outmarker}; "

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
