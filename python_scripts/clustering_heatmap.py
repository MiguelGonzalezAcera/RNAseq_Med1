import argparse
import logging
import os
import glob
import json
import pandas as pd


def clustering_heatmap(config, tool_name):
    """Get the
    """

    out_dir = "/".join(config['tools_conf'][tool_name]['output']['heatmap'].split('/')[0:-1])

    norm_counts = config['tools_conf'][tool_name]['input']['norm_counts']
    norm_counts_tab = norm_counts.replace(".Rda",".tsv")
    heatmap = config['tools_conf'][tool_name]['output']['heatmap']
    norm_counts_res = heatmap.replace('_heatmap.png','_NormCounts.tsv')
    organism = config['options']['organism']

    # Create the command to run the pca R script
    command = ""
    if not os.path.exists(out_dir):
        command += f"mkdir {out_dir};"

    if 'genelist' in config['tools_conf'][tool_name]['input']:
        genelist_path = config['tools_conf'][tool_name]['input']['genelist']
    else:
        out_dir_DE = "/".join(config['tools_conf'][tool_name]['input']['DEtouched'].split('/')[0:-1])
        samples = config['comparisons']

        genes_list = []

        for control in samples:
            sample_ids = samples[control].split(",")
            for sample in sample_ids:
                id_sample = out_dir_DE + "/" + config['project'] + "_" + f"{sample}_{control}.tsv"

                tdf = pd.read_csv(id_sample, sep='\t', index_col=None)
                genes = tdf[tdf['padj'] < 0.05]["EnsGenes"].tolist()

                genes_list += genes

        genes_list = list(set(genes_list))

        genelist_path = f"{out_dir_DE}/significant_genes.txt"
        gene_file = open(genelist_path,"w")
        for gene in genes_list:
            gene_file.write(f"{gene}\n")
        gene_file.close()

    command += f'Rscript /DATA/RNAseq_test/Scripts/Rscripts/clustering.r --heatmap {heatmap} --counts {norm_counts} --genelist {genelist_path} --organism {organism}; '

    command += f'head -n +1 {norm_counts_tab} | awk \'{{print \"EnsemblID\\t\" $0}}\' > {norm_counts_res}; grep -f {genelist_path} {norm_counts_tab} >> {norm_counts_res}; '


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

    config = {'tools_conf': {'clustering_heatmap': config_dict}}
    config['options'] = config['tools_conf']['clustering_heatmap']['options']

    clustering_heatmap(config, 'clustering_heatmap')


if __name__ == "__main__":
    main()
