import argparse
import logging
import os
import pandas as pd
import python_scripts.python_functions as pf


def GSVA(config, tool_name):
    """Get the
    """
    # Start log
    logging.info(f'Starting {tool_name} process')

    # Get input parameters
    norm_counts = config['tools_conf'][tool_name]['input']['norm_counts']
    design = config['tools_conf'][tool_name]['input']['design']

    # define the out directory
    out_dir = "/".join(config['tools_conf'][tool_name]['output']['heatmap'].split('/')[0:-1])

    # Define the rest of the output parameters
    heatmap = config['tools_conf'][tool_name]['output']['heatmap']

    # Get options
    organism = config['options']['organism']
    project = config['project']
    samples = config['comparisons']

    # Create the command to run the pca R script
    command = ""

    # Make the outfolder if it doesnt wxist
    if not os.path.exists(out_dir):
        command += f"mkdir {out_dir};"

    # Loop through controls
    for control in samples:
        # Add the main command
        command += f'Rscript Rscripts/GSVA.r --heatmap {heatmap} --counts {norm_counts} --design {design} --organism {organism} --out_obj {out_dir}/{project}_GSVA.Rda --control {control} --comparisons {samples[control]}; '

    # Run the command and log it
    pf.run_command(command)