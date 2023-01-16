import argparse
import logging
import os
import glob
import json
import subprocess
import pandas as pd
import python_scripts.python_functions as pf


def markers_plots(norm_counts, heatmap, organism, command, design, dims):
    """
    Do plots for markers
    """
    #<TODO>: Establish sets of markers for specific tissues
    gene_markers = {
        "mouse": {
            "Mitochondrial": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/Mitochondrial_ensembl.txt",
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
        },
        "human": {
            "Mitochondrial": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/Mitochondrial_human_ensembl.txt",
            "EnterocyteDist": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/EnterocyteDist_human_ensembl.txt",
            "EnterocyteProx": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/EnterocyteProx_human_ensembl.txt",
            "Enteroendocrine": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/Enteroendocrine_human_ensembl.txt",
            "Goblet": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/Goblet_human_ensembl.txt",
            "M_cells": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/M_cells_human_ensembl.txt",
            "Paneth": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/Paneth_human_ensembl.txt",
            "StemProg": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/StemProg_human_ensembl.txt",
            "TAProg": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/TAProg_human_ensembl.txt",
            "Tuft": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/Tuft_human_ensembl.txt",
            "Fibroblasts": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/Fibroblasts_human_ensembl.txt",
            "MO_DC": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/MO_DC_human_ensembl.txt",
            "Plasma_cells": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/Plasma_cells_human_ensembl.txt",
            "T_cells": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/T_cells_human_ensembl.txt",
            "B_cells": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/B_cells_human_ensembl.txt",
            "Mast_cells": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/Mast_cells_human_ensembl.txt",
            "NK_ILC1": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/NK_ILC1_human_ensembl.txt",
            "Endothelial": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/Endothelial_human_ensembl.txt",
            "Neutrophils": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/Neutrophils_human_ensembl.txt",
            "Smooth_muscle": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/Smooth_muscle_human_ensembl.txt",
            "EntericGlial": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/EntericGlial_human_ensembl.txt",
            "EntericNeuron": "/DATA/Thesis_proj/Requests_n_stuff/20210224_Christoph_cell_markers/MarkersV2/human_markers/EntericNeuron_human_ensembl.txt"
        }
    }

    for marker in gene_markers[organism]:
        marker_path = gene_markers[organism][marker]
        heatmap_mark = heatmap.replace(".png",f"_{marker}_marker.png")
        command += f'Rscript Rscripts/clustering.r --heatmap {heatmap_mark} --counts {norm_counts} --genelist {marker_path} --organism {organism} --dims {dims} --design {design}; '

    return(command)


def clustering_markers(config, tool_name):
    """Get each clustering plot for the established markers
    """
    logging.info(f'Starting {tool_name} process')

    out_dir = "/".join(config['tools_conf'][tool_name]['output']['markerstouched'].split('/')[0:-1])

    project = config['project']
    norm_counts = config['tools_conf'][tool_name]['input']['norm_counts']
    heatmap = config['tools_conf'][tool_name]['output']['markerstouched'].replace("markerstouched.txt",f"{project}_clustering_markers.png")
    organism = config['options']['organism']
    design = config['tools_conf'][tool_name]['input']['design']
    dims = config['tools_conf'][tool_name]['tool_conf']['dimensions']

    # Create the command to run the pca R script
    command = ""
    if not os.path.exists(out_dir):
        command += f"mkdir {out_dir};"

    # Check if this is a marker run
    if config['pipeline'] == 'RNAseq':
        command = markers_plots(norm_counts, heatmap, organism, command, design, dims)

    command += f"touch {config['tools_conf'][tool_name]['output']['markerstouched']}; "

    logging.info(f'Running command: {command}')
    pf.run_command(command)

def get_arguments():
    """
    Function that parse arguments given by the user, returning a dictionary
    that contains all the values.
    """

    # Create the top-level parser
    parser = argparse.ArgumentParser()

    # Mandatory variables
    parser.add_argument('--counts', required=True, help='Table with the counts of the assay')
    parser.add_argument('--design', required=True, help='Table with the design of the experiment')
    parser.add_argument('--heatmap', required=True, help='PNG with the heatmap result')
    parser.add_argument('--project', required=True, help='Project name')
    parser.add_argument('--organism', required=True, help='Organism')
    parser.add_argument('--dims', default="2000,2000", help='Dimensions for the plot')

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

    config = {
      "DEBUG": args.debug,
      "TESTING": args.test,
      "DRY_RUN": args.dry_run,
      "log_files": ["/tmp/full.log"],
      "project": args.project,
      "options": {
        "organism": args.organism,
        "sql_load": "True"
      },
      "tools_conf": {
        "clustering_markers": {
          "input": {
            "norm_counts": args.counts,
            "design": args.design
            },
          "output": {
            "heatmap": args.heatmap
            },
          "tool_conf": {
            "dimensions": args.dims
            }
          }
        }
      }

    # Startup the logger format
    logger = pf.create_logger(config['log_files'][0])

    clustering_markers(config, 'clustering_markers')

    logging.info(f'Finished clustering')

if __name__ == "__main__":
    main()
