import argparse
import logging
import os
import glob
import pandas as pd
import python_scripts.python_functions as pf


def clustering_heatmap(config, tool_name):
    """Get the
    """

    out_dir_DE = "/".join(config['tools_conf'][tool_name]['input']['DEtouched'].split('/')[0:-1])
    out_dir = "/".join(config['tools_conf'][tool_name]['output']['heatmap'].split('/')[0:-1])

    norm_counts = config['tools_conf'][tool_name]['input']['norm_counts']
    heatmap = config['tools_conf'][tool_name]['output']['heatmap']
    samples = config['comparisons']
    organism = config['options']['organism']

    # Create the command to run the pca R script
    command = ""
    if not os.path.exists(out_dir):
        command += f"mkdir {out_dir};"

    genes_list = []

    for control in samples:
        sample_ids = samples[control].split(",")
        for sample in sample_ids:
            id_sample = out_dir_DE + "/" + config['project'] + "_" + f"{sample}_{control}.tsv"

            tdf = pd.read_csv(id_sample, sep='\t', index_col=None)
            genes = tdf[tdf['padj'] < 0.05]["EnsGenes"].tolist()

            genes_list += genes

    genes_list = list(set(genes_list))

    gene_file = open(f"{out_dir_DE}/significant_genes.txt","w")
    for gene in genes_list:
        gene_file.write(f"{gene}\n")
    gene_file.close()

    command += f'Rscript /DATA/RNAseq_test/Scripts/Rscripts/clustering.r --heatmap {heatmap} --counts {norm_counts} --genelist {out_dir_DE}/significant_genes.txt --organism {organism}; '

    print(command)

    pf.run_command(command)


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

    config = {
      "DEBUG": args.debug,
      "TESTING": args.test,
      "DRY_RUN": args.dry_run,
      "log_files": ["/tmp/full.log"],
      "tools_conf": {
        "pca": {
          "input": {
            "counts": args.counts,
            "design": args.design
            },
          "output": {
            "out_dir": args.out_dir
            },
          "tool_conf": {
            }
          }
        }
      }

    # Startup the logger format
    logger = pf.create_logger(config['log_files'][0])

    pca(config, 'pca')


if __name__ == "__main__":
    main()
