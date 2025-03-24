import logging
import os
import json
import pandas as pd
import mysql.connector
import python_scripts.python_functions as pf

def fix_genelists(gene_markers, outpath, organism, mycursor):
    """"""
    # Create the dictionary with each marker
    marker_dict = {}

    # Iter through the markers to generate the dictionary
    for marker in gene_markers:
        # Get gene reference table
        if organism == 'mouse':
            marker_command = f'select * from markers_{marker};'
        elif organism == 'human':
            marker_command = f'select * from markers_{marker}_human;'

        # Execute the command to the database and retrieve the table into a dataframe
        mycursor.execute(marker_command)

        df_set = []
        for row in mycursor:
            df_set.append(row)
        resdf = pd.DataFrame(df_set)
        resdf.columns = ['ensembl','entrez','genename']

        # Add the ensembl ids as list to the dictionary
        marker_dict[marker] = resdf['ensembl'].tolist()

    # Generate the file name
    genegroups_path = f"{outpath}/genegroups.json"

    # dump onto json file
    with open(genegroups_path, 'w') as f:
        json.dump(marker_dict, f)

    # Return the path to the dictionary
    return genegroups_path

def GSEA_markers_plots(in_obj, outpath, organism, command, dims):
    """
    Do plots for markers
    """
    # Get the list of markers
    gene_markers = ["Mitochondrial", "EnterocyteDist", "EnterocyteProx", "Enteroendocrine", "Goblet", "Mcells", "Paneth",
        "Stem", "TAprog", "Tuft", "Fibroblasts", "MODC", "Plasma", "Tcells", "Bcells", "Mast", "NK",
        "Endothelial", "Neutrophils", "SmoothMuscle", "EntericGlia", "EntericNeuron"]
    
    # Establish a connection to the database for getting the ref lists
    mydb = mysql.connector.connect(
        host="localhost",
        user="root",
        passwd="Pl4ter!a",
        database="Refs"
    )

    mycursor = mydb.cursor()

    # Correct the markers into a json file
    genegroup = fix_genelists(gene_markers, outpath, organism, mycursor)

    # Make the model outfile
    gseaplot_mark = outpath + "/" + in_obj.split('/')[-1].replace(".Rda",f".svg")

    # Add the command
    command += f'Rscript Rscripts/GSEA.r --pathways {genegroup} --in_obj {in_obj} --gseaplot {gseaplot_mark} --organism {organism} --dims {dims}; '

    return(command)


def GSEA_markers(config, tool_name):
    """Get the counts of a number of bam files in a directory
    """
    logging.info(f'Starting {tool_name} process')

    # Extract the organism
    organism = config['options']['organism']

    # Extract the comparisons
    comparisons = config['comparisons']

    # Extract the infiles and the project name
    in_path = "/".join(config['tools_conf'][tool_name]['input']['DEtouched'].split('/')[0:-1])
    project = config['project']

    # Extract the out path and out markerfile
    outmarker = config['tools_conf'][tool_name]['output']['GSEAMtouched']
    out_path = "/".join(outmarker.split('/')[0:-1])

    # Get the dimensions
    dimensions = config['tools_conf'][tool_name]['tool_conf']['dimensions']

    # Create the command to run the pca R script
    command = ""

    # Create outfolder if it doesnt exist
    if not os.path.exists(out_path):
        command += f"mkdir {out_path};"

    # Loop through the samples and controls
    for control in comparisons:
        samples = comparisons[control].split(",")
        for sample in samples:
            # Read the table with the whole data
            in_obj = f"{in_path}/{project}_{sample}_{control}.Rda"

            # Create the corresponding command
            command = GSEA_markers_plots(in_obj, out_path, organism, command, dimensions)

    # touch the marker file
    command += f"touch {outmarker}; "

    # Run the finished command
    pf.run_command(command)