import logging
import os
import pandas as pd
import python_scripts.python_functions as pf

def GO_enrichment(config, tool_name):
    """Get the counts of a number of bam files in a directory
    """

    # Record in logger
    logging.info(f'Starting {tool_name} process')

    organism = config['options']['organism']

    # Create the command to run the R script
    command = ""

    # Make all the paths and names for the pipeline run
    # Get the directory with the differential expression files
    out_dir_DE = "/".join(config['tools_conf'][tool_name]['input']['DEtouched'].split('/')[0:-1])

    # Extract the markerfile path
    GOtouched = config['tools_conf'][tool_name]['output']['gotouched']

    # Extract the dictionary with the comparisons
    samples = config['comparisons']

    # Get the directory to store the GO results from the markerfile. Make it if it does not exist
    out_dir = "/".join(GOtouched.split('/')[0:-1])

    if not os.path.exists(out_dir):
        command += f"mkdir {out_dir};"

    # Iter through the samples and controls to make a GO analysis for each comparison
    for control in samples:
        sample_ids = samples[control].split(",")
        for sample in sample_ids:
            # Get the name of the out dir
            id_dir = out_dir + "/" + f"{sample}_{control}"
            # Make the name of the outfile. The script witll generate the ontology variants
            id_tab = id_dir + "/" + f"{sample}_{control}_GO.tsv"
            # Get the DE file
            id_sample = out_dir_DE + "/" + config['project'] + "_" + f"{sample}_{control}.Rda"
            # Get the general universe file. Required for the GO function in R
            id_universe = id_sample + "_universe.Rda"

            # Create the out firectory and the GO R command
            command += f'Rscript Rscripts/GO_enrichment.r --out_tab {id_tab} --obj {id_sample} --universe {id_universe} --organism {organism}; '
    
    # create the markerfile
    command += f'touch {GOtouched}'

    # Run the command
    pf.run_command(command)