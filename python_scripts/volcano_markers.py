import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import argparse
import logging
import os
import json
import subprocess
import math
import pandas as pd
import python_scripts.python_functions as pf

# Disable warning
pd.options.mode.chained_assignment = None  # default='warn'


def get_gene_markers(organism):
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

    return(gene_markers[organism])


def volcano_markers(config, tool_name):
    """Get the scatter plots for the markers
    """

    # Start the log for the process
    logging.info(f'Starting {tool_name} process')

    # Extract the organism
    organism = config['options']['organism']

    # Select the list of markers
    gene_markers = get_gene_markers(organism)

    # Extract the comparisons
    comparisons = config['comparisons']

    # Extract the infiles and the project name
    in_path = "/".join(config['tools_conf'][tool_name]['input']['DEtouched'].split('/')[0:-1])
    project = config['project']

    # Extract the out path and out markerfile
    out_path = "/".join(config['tools_conf'][tool_name]['output']['MVtouched'].split('/')[0:-1])
    outmarker = config['tools_conf'][tool_name]['output']['MVtouched']

    # Create out folder and marker
    command = ""
    if not os.path.exists(out_path):
        command += f"mkdir {out_path};"

    command += f"touch {outmarker}; "

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

    # Loop through the samples and controls
    for control in comparisons:
        samples = comparisons[control].split(",")
        for sample in samples:
            # Read the table with the whole data
            tab_name = f"{in_path}/{project}_{sample}_{control}.tsv"
            df = pd.read_csv(tab_name, sep='\t')

            # Iterate throug the genelists
            for glist in gene_markers:
                # Read the gene list into a list
                with open(gene_markers[glist]) as f:
                    genelist = f.read().splitlines()

                # Filter dataframe by the genelist
                wdf = df[df['EnsGenes'].isin(genelist)]

                if wdf.empty:
                    logging.info(f'Could not find genes expressed for marker {glist}')
                    continue

                # Modify the pvalue column
                wdf['Mod_pvalue'] = [-math.log10(i+1e-148) for i in wdf['pvalue'].tolist()]

                # Create the figure
                fig, ax = plt.subplots(figsize=(7,7))

                # Establish the limits for the plot
                ## Y top limit
                if max(wdf['Mod_pvalue'].tolist()) < 10:
                    ytoplim = 10
                else:
                    ytoplim = max(wdf['Mod_pvalue'].tolist())

                ## X left limit
                if min(wdf['log2FoldChange'].tolist()) > -5:
                    xleftlim = -5
                else:
                    xleftlim = min(wdf['log2FoldChange'].tolist())*1.1

                ## Y left limit
                if max(wdf['log2FoldChange'].tolist()) < 5:
                    xrightlim = 5
                else:
                    xrightlim = max(wdf['log2FoldChange'].tolist())*1.1

                # Plot the rectangles
                boxDR = plt.Rectangle(((xleftlim*1.1)-0.5, -math.log10(0.05)),
                                     abs(xleftlim*1.1), ytoplim*2, color="blue", alpha=0.2)
                ax.add_patch(boxDR)

                boxUR = plt.Rectangle((0.5, -math.log10(0.05)),
                                     xrightlim*2, ytoplim*2, color="red", alpha=0.2)
                ax.add_patch(boxUR)

                # Add the initial scatterplot
                ax.scatter(wdf['log2FoldChange'].tolist(), wdf['Mod_pvalue'].tolist(), s=100, c="grey", alpha=0.5)

                # Add the up and down regulated scatters
                wdf_filt1 = wdf[(wdf['log2FoldChange'] < -0.5) & (wdf['pvalue'] < 0.05)]
                ax.scatter(wdf_filt1['log2FoldChange'].tolist(), wdf_filt1['Mod_pvalue'].tolist(),
                           s=200, c="blue")

                wdf_filt2 = wdf[(wdf['log2FoldChange'] > 0.5) & (wdf['pvalue'] < 0.05)]
                ax.scatter(wdf_filt2['log2FoldChange'].tolist(), wdf_filt2['Mod_pvalue'].tolist(),
                           s=200, c="red")

                # With the subsets, calculate and plot percentages of genes up/down
                perc_down = (len(wdf_filt1)/len(wdf))*100
                ax.text(xleftlim*0.9, ytoplim*1.35, f"{perc_down:.1f}%", fontsize=30, c="blue")

                perc_up = (len(wdf_filt2)/len(wdf))*100
                ax.text(xrightlim*0.4, ytoplim*1.35, f"{perc_up:.1f}%", fontsize=30, c="red")

                # Set axis limits
                ax.set_ylim(bottom=0, top=ytoplim*1.5)
                ax.set_xlim(xleftlim, xrightlim)

                # Set axis labels
                ax.set_xlabel("Fold Change (log2FoldChange)", {'fontsize': 20})
                ax.set_ylabel("P value (-log10(pvalue))", {'fontsize': 20})

                # Set plot title
                ax.set_title(f"{glist} - {sample} - {control}", {'fontsize': 25})

                fig.savefig(f'{out_path}/{glist}_{sample}_{control}_scattermarkers.png')
                plt.close()

def get_arguments():
    """
    Function that parse arguments given by the user, returning a dictionary
    that contains all the values.
    """

    # Create the top-level parser
    parser = argparse.ArgumentParser()

    # Mandatory variables
    parser.add_argument('--DEtouched', required=True, help='Marker file for the differential expression')
    parser.add_argument('--MVtouched', required=True, help='Marker file for the end of the process')
    parser.add_argument('--project', required=True, help='Project name')
    parser.add_argument('--organism', required=True, help='Organism')

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

    out_dir = "/".join(args.DEtouched.split('/')[0:-1])
    project = args.project

    config = {
      "DEBUG": args.debug,
      "TESTING": args.test,
      "DRY_RUN": args.dry_run,
      "log_files": ["/tmp/full.log"],
      "project": project,
      "options": {
        "organism": args.organism,
        "sql_load": "True"
      },
      "comparisons": {
		"KO_H2O": "KO_IL22,KO_IL_A",
        "WT_H2O": "WT_IL22,WT_IL_A,KO_H2O",
        "WT_IL22": "WT_IL_A,KO_IL22",
        "WT_IL_A": "KO_IL_A",
        "KO_IL22": "KO_IL_A"
	  },
      "tools_conf": {
        "volcano_markers": {
          "input": {
            "DEtouched": args.DEtouched
            },
          "output": {
            "MVtouched": args.MVtouched
            },
          "tool_conf": {
            }
          }
        }
      }

    # Startup the logger format
    logger = pf.create_logger(config['log_files'][0])

    volcano_markers(config, 'volcano_markers')

    logging.info(f'Finished volcano_markers')

if __name__ == "__main__":
    main()
