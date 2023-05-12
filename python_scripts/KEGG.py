import logging
import os
import pandas as pd
import python_scripts.python_functions as pf


def KEGG_enrichment(config, tool_name):
    """Get the counts of a number of bam files in a directory
    """

    logging.info(f'Starting {tool_name} process')

    # Get the organism
    organism = config['options']['organism']

    # Create the command to run the pca R script
    command = ""

    # Get directory of the DE files
    out_dir_DE = "/".join(config['tools_conf'][tool_name]['input']['DEtouched'].split('/')[0:-1])

    # Define the directory where the results will be stored
    out_dir = "/".join(config['tools_conf'][tool_name]['output']['keggtouched'].split('/')[0:-1])

    # Get the marker for end
    keggtouched = config['tools_conf'][tool_name]['output']['keggtouched']

    # Get the comparisons to use
    samples = config['comparisons']

    # Make the out dir if it doesnt exist
    if not os.path.exists(out_dir):
        command += f"mkdir {out_dir}; "

    # Iter through the samples to make the analysis
    for control in samples:
        sample_ids = samples[control].split(",")
        for sample in sample_ids:
            # Get an identifyer for the object
            id_obj = f"{sample}_{control}"

            # define a separate folder for each analysis
            id_dir = out_dir + "/" + id_obj

            # Create the folder in question
            if not os.path.exists(id_dir):
                command += f"mkdir {id_dir}; "

            # Define the name of the out table
            id_tab = id_dir + "/" + f"{id_obj}_KEGG.tsv"

            # Obtain the path for the R object containing the DE result
            id_sample = out_dir_DE + "/" + config['project'] + "_" + f"{id_obj}.Rda"

            # Make the commands
            command += f'Rscript Rscripts/KEGG_enrichment.r --out_tab {id_tab} --in_obj {id_sample} --id {id_obj} --organism {organism}; '
    
    # make the command for the ending control file
    command += f'touch {keggtouched}'

    pf.run_command(command)