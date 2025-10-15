import logging
import os
import pandas as pd
import python_scripts.python_functions as pf

def dexseq(config, tool_name):
    """Get the counts of a number of bam files in a directory
    """
    # Note in logger
    logging.info(f'Starting {tool_name} process')

    # Extract information
    # Inputs
    # counts file
    counts_dir = "/".join(config['tools_conf'][tool_name]['input']['counts'].split('/')[0:-1])
    # design file
    design = config['tools_conf'][tool_name]['input']['design']

    # Outputs
    # marker file
    DEXtouched = config['tools_conf'][tool_name]['output']['DEXtouched']
    # outside directory
    out_dir = "/".join(DEXtouched.split('/')[0:-1])

    # Other information
    # comparison set
    samples = config['comparisons']
    # project name
    project = config['project']
    # organism used
    organism = config['options']['organism']
    # tx2gene
    annotation = config['tools_conf'][tool_name]['input']['annotation']

    # Create the command to run the dexseq R script
    command = ""

    # create the main directory if it doesnt exist
    if not os.path.exists(out_dir):
        command += f"mkdir {out_dir};"

    # Create the deseq2 command for each control and sample
    for control in samples:
        sample_IDs = samples[control].split(',')
        for sample_name in sample_IDs:

            out_dir_sample = "/".join([out_dir, f"DTU_{sample_name}_{control}"])

            if not os.path.exists(out_dir_sample):
                command += f"mkdir {out_dir_sample};"

            # out object name
            out_obj = out_dir_sample + "/" + config['project'] +".Rda"

            command += f'Rscript Rscripts/dexseq.r --salmon_counts {counts_dir} --annotation {annotation} --design {design} --out_obj {out_obj} --organism {organism} --control {control} --comparison {sample_name}; '

    # Touch the markerfile
    command += f'touch {DEXtouched}; '

    # Run the commans
    pf.run_command(command)