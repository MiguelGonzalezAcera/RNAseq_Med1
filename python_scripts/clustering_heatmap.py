import argparse
import logging
import os
import pandas as pd
import python_scripts.python_functions as pf


def clustering_heatmap(config, tool_name):
    # Report the start in the log
    logging.info(f'Starting {tool_name} process')

    # Extract the inputs to variables
    # R object with the counts
    norm_counts = config['tools_conf'][tool_name]['input']['norm_counts']
    # Table with the counts
    norm_counts_tab = norm_counts.replace(".Rda",".tsv")
    # Design file
    design = config['tools_conf'][tool_name]['input']['design']

    # Extract the outputs
    # Result heatmap path
    heatmap = config['tools_conf'][tool_name]['output']['heatmap']
    # Result path for accessory table
    norm_counts_res = heatmap.replace('_heatmap.png','_NormCounts.tsv')
    # Out directory
    out_dir = "/".join(heatmap.split('/')[0:-1])

    # Extract the other parameters
    # Organism
    organism = config['options']['organism']
    # Dimensions of the plot
    dimensions = config['tools_conf'][tool_name]['tool_conf']['dimensions']

    # Create the command to run the pca R script
    command = ""
    # Make the outdir if it doesn't exist
    if not os.path.exists(out_dir):
        command += f"mkdir {out_dir};"

    # Get the genelist, if provided, or extract one using simple filters from the DE files
    # Get the DE path from inputs
    out_dir_DE = "/".join(config['tools_conf'][tool_name]['input']['DEtouched'].split('/')[0:-1])
    # Extract the design dictionary
    samples = config['comparisons']

    # Init the list of genes
    genes_list = []

    # Iter through the samples and controls
    for control in samples:
        sample_ids = samples[control].split(",")
        for sample in sample_ids:
            # Construct the path for the DE table using the path and the project name
            id_sample = out_dir_DE + "/" + config['project'] + "_" + f"{sample}_{control}.tsv"

            # Load and select the genes that are significant
            tdf = pd.read_csv(id_sample, sep='\t', index_col=None)
            genes = tdf[tdf['padj'] < 0.05]["EnsGenes"].tolist()

            # Add the selection to the list of genes
            genes_list += genes

    # Remove duplicated genes
    genes_list = list(set(genes_list))

    # Save the list of genes to a file, to grep it later
    genelist_path = f"{out_dir_DE}/significant_genes.txt"
    gene_file = open(genelist_path,"w")
    for gene in genes_list:
        gene_file.write(f"{gene}\n")
    gene_file.close()

    # Make the command for the clustering
    command += f'Rscript Rscripts/clustering.r --heatmap {heatmap} --counts {norm_counts} --genelist {genelist_path} --organism {organism} --design {design} --dims {dimensions}; '

    # Slice the table with the selected genes
    command += f'head -n +1 {norm_counts_tab} | awk \'{{print \"EnsemblID\\t\" $0}}\' > {norm_counts_res}; grep -f {genelist_path} {norm_counts_tab} >> {norm_counts_res}; '

    # log and run the command
    logging.info(f'Running command: {command}')
    pf.run_command(command)