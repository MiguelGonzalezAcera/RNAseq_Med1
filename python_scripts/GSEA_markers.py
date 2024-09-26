import logging
import os
import pandas as pd
import mysql.connector
import python_scripts.python_functions as pf

def fix_genelists(marker, outpath, organism, mycursor):
    """"""
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
    resdf.columns = ["ensembl','entrez','genename"]

    # Add a column with the marker name
    resdf['group'] = [marker]*len(resdf['ensembl'].tolist())

    # Change column names and data types
    resdf.columns = ['group','entrez']

    # Remove rows with Nan
    resdf = resdf.dropna()
    resdf['entrez'] = resdf['entrez'].astype('int')

    # Dump onto file
    genegroups_path = f"{outpath}/{marker}_genegroups_test.txt"
    resdf.to_csv(genegroups_path, sep='\t', index=False, header=False)

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

    # Loop through the markers
    for marker in gene_markers:
        # Make the marker outfile
        gseaplot_mark = outpath + "/" + in_obj.split('/')[-1].replace(".Rda",f"_{marker}_GSEA.png")

        # Fix the markers into a gene group
        # Explanation: The GSEA function takes a named list of the Entrez IDs, so we must transform the gene markers to something more akin
        genegroup = fix_genelists(marker, outpath, organism, mycursor)

        # Add command
        command += f'Rscript Rscripts/GSEA.r --genegroup {genegroup} --in_obj {in_obj} --gseaplot {gseaplot_mark} --organism {organism} --dims {dims}; '

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