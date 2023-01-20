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
    #<TODO>: This could be read from the mysql db (here it makes sense)
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

    return(gene_markers[organism])

def volcano_marker_plot(tab_name, gene_markers, out_path, title):
    """ Function to produce the plot
    """
    # Read the provided table
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
        ax.set_title(title.replace('GENEMARKER', glist), {'fontsize': 25})

        fig.savefig(out_path.replace('GENEMARKER', glist))
        plt.close()



def volcano_markers(config, tool_name):
    """Get the scatter plots for the markers
    """

    # Start the log for the process
    logging.info(f'Starting {tool_name} process')

    # Extract the organism
    organism = config['options']['organism']

    # Select the list of markers
    gene_markers = get_gene_markers(organism)

    if 'MVtouched' in config['tools_conf'][tool_name]['output']:
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

        pf.run_command(command)

        # Loop through the samples and controls
        for control in comparisons:
            samples = comparisons[control].split(",")
            for sample in samples:
                # Get the paths for in table, out plot format. Use wildcard /GENENAME/ in the name of the plot
                tab_name = f"{in_path}/{project}_{sample}_{control}.tsv"
                out_path = f'{out_path}/GENEMARKER_{sample}_{control}_scattermarkers.png'

                # Establish the plot title
                plot_title = f"GENEMARKER - {sample} - {control}"

                # Run the plots
                volcano_marker_plot(tab_name, gene_markers, out_path, plot_title)
    
    else:
        # get the path for the table
        in_path = config['tools_conf'][tool_name]['input']['RData']

        # get the path for the folder that's going to contain the plots.
        out_dir = config['tools_conf'][tool_name]['output']['out_dir']

        #  Make it if it doesn't exist
        if not os.path.exists(out_dir):
            pf.run_command(f"mkdir {out_dir};")

        # Make out a path with a wildcard for the plots
        out_path = f'{out_dir}/GENEMARKER_scattermarkers.png'

        # Make out a title with the wildcard for the same purpose
        plot_title = f"Marker genes for GENEMARKER"

        volcano_marker_plot(tab_name, gene_markers, out_path, plot_title)
          

def get_arguments():
    """
    Function that parse arguments given by the user, returning a dictionary
    that contains all the values.
    """

    # Create the top-level parser
    parser = argparse.ArgumentParser()

    # Mandatory variables
    parser.add_argument('--RData', required=True, help='R object with the result of the differential expression')
    parser.add_argument('--out_dir', required=True, help="Folder that will contain the plots. Will create it if it doesn't exist")
    parser.add_argument('--organism', required=True, help='organism of the input data')

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
        "options" : {
            "organism": args.organism
        },
        "tools_conf": {
            "volcano_plot": {
            "input": {
                "RData": args.RData
                },
            "output": {
                "out_dir": args.out_dir,
                },
            "tool_conf": {}
          }
        }
      }

    # Startup the logger format
    logger = pf.create_logger(config['log_files'][0])

    volcano_markers(config, 'volcano_markers')

    logging.info(f'Finished volcano_markers')

if __name__ == "__main__":
    main()
