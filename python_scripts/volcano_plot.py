import argparse
import logging
import os
import math
import pandas as pd
import python_scripts.python_functions as pf
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

def volcano(df_path, plot_path, genelist=[], dims=["6000","3600"], labels=False):
  """
  Make the actual plot and save it
  """

  # Read the table file with the result of the differential expression
  df = pd.read_csv(df_path, sep='\t')

  # Create the figure
  fig, ax = plt.subplots(figsize=(int(dims[0])/600,int(dims[1])/600))

  # Fix and plot the p values
  ax.scatter(df['log2FoldChange'].tolist(),[-math.log10(float(i)+1e-30) for i in df['padj'].tolist()], c='grey', s=2, alpha=0.05)

  # Plot threshold lines
  ax.axvline(x=1, c='orange', alpha=0.25)
  ax.axvline(x=-1, c='orange', alpha=0.25)
  ax.axhline(y=-math.log10(0.05), c='orange', alpha=0.25)

  # Put labels on axises
  ax.set_xlabel("Fold Change (log2FoldChange)")
  ax.set_ylabel("P value (-log10(pvalue))")

  # Draw the dots of 
  if genelist:
    # Show marked genes
    wdf = df[df['Genes'].isin(genelist)]

    # Plot the selection of genes
    ax.scatter(wdf['log2FoldChange'].tolist(),[-math.log10(float(i)+1e-30) for i in wdf['padj'].tolist()], c='blue', s=100)

    # Label the dots
    if labels:
      for i, label in enumerate(df['Genes'].tolist()):
        plt.annotate(label, (wdf['log2FoldChange'].tolist()[i] + 0.15, [-math.log10(float(i)+1e-30) for i in wdf['padj'].tolist()][i] + 0.15), fontsize=20)
  
  # Plot colors if there is no genelist
  else:
    # filter the dataframe by the conditions of the sections
    # Non-significant
    nsdf = df[(df['padj'] > 0.05) & ((df['log2FoldChange'] > 1) | (df['log2FoldChange'] < -1))]
    # low FC
    lfcdf = df[(df['padj'] < 0.05) & ((df['log2FoldChange'] < 1) | (df['log2FoldChange'] > -1))]
    # significant n high fc
    signdf = df[(df['padj'] < 0.05) & ((df['log2FoldChange'] > 1) | (df['log2FoldChange'] < -1))]

    # Draw selected sections
    ax.scatter(nsdf['log2FoldChange'].tolist(),[-math.log10(float(i)+1e-30) for i in nsdf['padj'].tolist()], c='orange', s=50, alpha=0.3)
    ax.scatter(lfcdf['log2FoldChange'].tolist(),[-math.log10(float(i)+1e-30) for i in lfcdf['padj'].tolist()], c='red', s=50, alpha=0.3)
    ax.scatter(signdf['log2FoldChange'].tolist(),[-math.log10(float(i)+1e-30) for i in signdf['padj'].tolist()], c='green', s=50)

  fig.tight_layout()

  plt.savefig(plot_path, dpi=600)

def volcano_plot(config, tool_name):
  """Framework to make the plot
  """

  # Start the logger
  logging.info(f'Starting {tool_name} process')

  # Get common parameters from the config
  # Get the dimension char
  dimensions = config['tools_conf'][tool_name]['tool_conf']['dimensions']
  #<TODO>: Raise an error if the dimensions introduced are not numbers
  dimensions = dimensions.split(',')

  # check if this run belongs in a pipeline or runs alone
  if 'volcanotouched' in config['tools_conf'][tool_name]['output']:
    # Get the directory of the DE files
    out_dir_DE = "/".join(config['tools_conf'][tool_name]['input']['DEtouched'].split('/')[0:-1])

    # Get the control file
    volcanotouched = config['tools_conf'][tool_name]['output']['volcanotouched']

    # Get the directory for the result plot
    out_dir = "/".join(volcanotouched.split('/')[0:-1])

    # Obtain the dictionary with the comparisons
    samples = config['comparisons']

    # Iter through the sample - control dictionary to generate the plots
    for control in samples:
        sample_ids = samples[control].split(",")
        for sample in sample_ids:
            # Make the name of the plot
            out_plot = out_dir + "/" + f"{sample}_{control}_volcano.png"

            # Get the name of the DE file (standardized)
            id_sample = out_dir_DE + "/" + config['project'] + "_" + f"{sample}_{control}.tsv"

            volcano(id_sample, out_plot, dims=dimensions)

    # Create the command to run the control in the pipeline
    # Make the outfolder if it doesn't exists
    if not os.path.exists(out_dir):
      command = f"mkdir {out_dir};"
    else:
      command = ""

    # Make the control file
    command += f'touch {volcanotouched};'

    # Run the command
    pf.run_command(command)

  else:
    # Get the path for the outfile 
    out_plot = config['tools_conf'][tool_name]['output']['volcano']

    # Get the tab for tie DE table
    id_tab = config['tools_conf'][tool_name]['input']['RData']

    # Get the out directory
    out_dir = "/".join(out_plot.split('/')[0:-1])

    # Read the list of genes, if any
    if config['tools_conf'][tool_name]['input']['genelist'] != "":
      # Read the lines of the file
      with open(config['tools_conf'][tool_name]['input']['genelist']) as f:
        gene_list = f.read().splitlines()
    else:
      gene_list = []

    # get the labels option
    labels = config['tools_conf'][tool_name]['tool_conf']['labels']

    # Do the plot
    volcano(id_tab, out_plot, genelist=gene_list, dims=dimensions, labels=labels)

def get_arguments():
    """
    Function that parse arguments given by the user, returning a dictionary
    that contains all the values.
    """

    # Create the top-level parser
    parser = argparse.ArgumentParser()

    # Mandatory variables
    parser.add_argument('--RData', required=True, help='R object with the result of the differential expression')
    parser.add_argument('--genelist', default="", help='Particular genelist to highlight in the plot')
    parser.add_argument('--labels', default=False, choices=[True,False], help='Show names of the genes. Only works with specific genelist')
    parser.add_argument('--dims', default="6000,3600", help='Dimensions for the plot')
    parser.add_argument('--volcano', required=True, help='Result file where to store the plot')

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
        "volcano_plot": {
          "input": {
            "RData": args.RData,
            "genelist": args.genelist
            },
          "output": {
            "volcano": args.volcano,
            },
          "tool_conf": {
            "labels": args.labels,
            "dimensions": args.dims
            }
          }
        }
      }

    # Startup the logger format
    logger = pf.create_logger(config['log_files'][0])

    volcano_plot(config, 'volcano_plot')

if __name__ == "__main__":
    main()
