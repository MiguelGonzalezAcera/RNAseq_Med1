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
    #<TODO>: I think the first option only makes sense if it is a standalone run, but I think this was to be removed, so...
    if 'genelist' in config['tools_conf'][tool_name]['input']:
        # Get the genelist file from inputs
        genelist_path = config['tools_conf'][tool_name]['input']['genelist']
    else:
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

def get_arguments():
    """
    Function that parse arguments given by the user, returning a dictionary
    that contains all the values.
    """

    # Create the top-level parser
    parser = argparse.ArgumentParser()

    # Mandatory variables
    parser.add_argument('--genelist', required=True, help='List of genes to represent in the heatmap')
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
        "clustering_heatmap": {
          "input": {
            "genelist": args.genelist,
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

    clustering_heatmap(config, 'clustering_heatmap')

    logging.info(f'Finished clustering')

if __name__ == "__main__":
    main()
