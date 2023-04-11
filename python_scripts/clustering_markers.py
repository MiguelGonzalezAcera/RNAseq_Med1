import logging
import os
import mysql.connector
import pandas as pd
import python_scripts.python_functions as pf


def markers_plots(norm_counts, heatmap, organism, command, design, dims):
    """
    Do plots for markers
    """
    # Set a connection for the database to retrieve the markers
    mydb = mysql.connector.connect(
        host="localhost",
        user="root",
        passwd="Pl4ter!a",
        database="Refs"
        )

    mycursor = mydb.cursor()

    # Get the paths for the markers
    #<TODO>: Establish sets of markers for specific tissues
    gene_markers = ["Mitochondrial", "EnterocyteDist", "EnterocyteProx", "Enteroendocrine", "Goblet", "Mcells", "Paneth",
            "Stem", "TAprog", "Tuft", "Fibroblasts", "MODC", "Plasma", "Tcells", "Bcells", "Mast", "NK",
            "Endothelial", "Neutrophils", "SmoothMuscle", "EntericGlial", "EntericNeuron"]

    # Generate a clustering R command for each marker and genelist
    for marker in gene_markers:
        # Ask the database for the marker in question
        if organism == 'mouse':
            marker_command = f'select * from markers_{marker};'
        elif organism == 'human':
            marker_command = f'select * from markers_{marker}_human;'

        # Execute the command to the database and retrieve the table into a dataframe
        mycursor.execute(marker_command)

        df_set = []
        for row in mycursor:
            df_set.append(row)

        df = pd.DataFrame(df_set)
        df.columns = ["EnsGenes", "EntrezID", "Genes"]

        # Select the ensembl column
        df_ens = df['EnsGenes'].tolist()

        # Dump the dataframe into a temporal file
        marker_path = heatmap.replace(".png",f"_{marker}_genelist.txt")
        with open(marker_path, 'w') as f:
            for item in df_ens:
                f.write("%s\n" % item)

        # Heatmap path
        heatmap_mark = heatmap.replace(".png",f"_{marker}_marker.png")
        # R command
        command += f'Rscript Rscripts/clustering.r --heatmap {heatmap_mark} --counts {norm_counts} --genelist {marker_path} --organism {organism} --dims {dims} --design {design}; '

    # Commit the changes and close the database
    mydb.commit()

    mycursor.close()
    mydb.close()
    
    return(command)


def clustering_markers(config, tool_name):
    """Get each clustering plot for the established markers
    """
    logging.info(f'Starting {tool_name} process')

    # Get the needed parameters
    # Inputs
    # Normalized counts file
    norm_counts = config['tools_conf'][tool_name]['input']['norm_counts']
    # Design file
    design = config['tools_conf'][tool_name]['input']['design']

    # Outputs
    # Marker file
    markerstouched = config['tools_conf'][tool_name]['output']['markerstouched']
    # Heatmap file
    heatmap = markerstouched.replace("markerstouched.txt",f"{project}_clustering_markers.png")
    # Results directory
    out_dir = "/".join(markerstouched.split('/')[0:-1])
    
    # Other
    # Project name
    project = config['project']
    # Organism
    organism = config['options']['organism']
    # Dimensions
    dims = config['tools_conf'][tool_name]['tool_conf']['dimensions']

    # Create the command to run the pca R script
    command = ""
    if not os.path.exists(out_dir):
        command += f"mkdir {out_dir};"

    # Add the instances for each marker
    command += markers_plots(norm_counts, heatmap, organism, command, design, dims)

    # Touch the marker file
    command += f"touch {markerstouched}; "

    # Run the command
    logging.info(f'Running command: {command}')
    pf.run_command(command)