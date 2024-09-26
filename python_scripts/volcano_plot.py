import logging
import os
import math
import pandas as pd
import python_scripts.python_functions as pf
import matplotlib.pyplot as plt


def volcano(df_path, plot_path, dims=["10","6"]):
  """
  Make the actual plot and save it
  """

  # Read the table file with the result of the differential expression
  df = pd.read_csv(df_path, sep='\t')

  # Create the figure
  fig, ax = plt.subplots(figsize=(int(dims[0]),int(dims[1])))

  # Plot the values
  ax.scatter(df['log2FoldChange'].tolist(),[-math.log10(float(i)+1e-320) for i in df['padj'].tolist()], c='grey', s=2, alpha=0.05)

  # Plot the up n down values
  df_up = df[(df['padj'] < 0.05) & (df['log2FoldChange'] > 0)]
  ax.scatter(
      df_up['log2FoldChange'].tolist(),
      [-math.log10(float(i)+1e-320) for i in df_up['padj'].tolist()], 
      c='orange', 
      s=[-math.log10(float(i)+1e-320) for i in df_up['padj'].tolist()], 
      alpha=0.5
  )

  df_dw = df[(df['padj'] < 0.05) & (df['log2FoldChange'] < 0)]
  ax.scatter(
      df_dw['log2FoldChange'].tolist(),
      [-math.log10(float(i)+1e-320) for i in df_dw['padj'].tolist()], 
      c='purple', 
      s=[-math.log10(float(i)+1e-320) for i in df_dw['padj'].tolist()], 
      alpha=0.5
  )

  # Plot threshold lines
  ax.axvline(x=1, c='orange', alpha=0.25)
  ax.axvline(x=-1, c='orange', alpha=0.25)
  ax.axhline(y=-math.log10(0.05), c='orange', alpha=0.25)

  # Put labels on axises
  ax.set_xlabel("Fold Change (log2FoldChange)")
  ax.set_ylabel("P value (-log10(padj))")
  
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
