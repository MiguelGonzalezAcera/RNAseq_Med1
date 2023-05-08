import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import logging
import os
import math
import pandas as pd
import mysql.connector
import python_scripts.python_functions as pf

# Disable warning
pd.options.mode.chained_assignment = None  # default='warn'

def get_gene_markers(organism, tmppath):
    #<TODO>: This could be read from the mysql db (here it makes sense)
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
            "Endothelial", "Neutrophils", "SmoothMuscle", "EntericGlia", "EntericNeuron"]

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
        marker_path = f"{tmppath}/volcano_{marker}_genelist.txt"
        with open(marker_path, 'w') as f:
            for item in df_ens:
                f.write("%s\n" % item)

    # Commit the changes and close the database
    mydb.commit()

    mycursor.close()
    mydb.close()

    return(gene_markers)

def volcano_marker_plot(tab_name, gene_markers, out_path, title, dims=["7","7"]):
    """ Function to produce the plot
    """
    # Read the provided table
    df = pd.read_csv(tab_name, sep='\t')

    # Get the path for the files
    out_path_folder = "/".join(out_path.split('/')[0:-1])

    # Iterate throug the genelists
    for glist in gene_markers:
        # Read the gene list into a list
        with open(f"{out_path_folder}/volcano_{glist}_genelist.txt") as f:
            genelist = f.read().splitlines()

        # Filter dataframe by the genelist
        wdf = df[df['EnsGenes'].isin(genelist)]

        if wdf.empty:
            logging.info(f'Could not find genes expressed for marker {glist}')
            continue

        # Modify the pvalue column
        wdf['Mod_pvalue'] = [-math.log10(i+1e-148) for i in wdf['pvalue'].tolist()]

        # Create the figure
        fig, ax = plt.subplots(figsize=(int(dims[0]),int(dims[1])))

        # Establish the limits for the plot
        ## Y top limit
        if max(wdf['Mod_pvalue'].tolist()) < 10:
            ytoplim = 10
        else:
            ytoplim = max(wdf['Mod_pvalue'].tolist())

        ## X left limit
        if min(wdf['log2FoldChange'].tolist()) > -5:
            xleftlim = -5
        else:
            xleftlim = min(wdf['log2FoldChange'].tolist())*1.1

        ## Y left limit
        if max(wdf['log2FoldChange'].tolist()) < 5:
            xrightlim = 5
        else:
            xrightlim = max(wdf['log2FoldChange'].tolist())*1.1

        # Plot the rectangles
        boxDR = plt.Rectangle(((xleftlim*1.1)-0.5, -math.log10(0.05)),
                                abs(xleftlim*1.1), ytoplim*2, color="blue", alpha=0.2)
        ax.add_patch(boxDR)

        boxUR = plt.Rectangle((0.5, -math.log10(0.05)),
                                xrightlim*2, ytoplim*2, color="red", alpha=0.2)
        ax.add_patch(boxUR)

        # Add the initial scatterplot
        ax.scatter(wdf['log2FoldChange'].tolist(), wdf['Mod_pvalue'].tolist(), s=100, c="grey", alpha=0.5)

        # Add the up and down regulated scatters
        wdf_filt1 = wdf[(wdf['log2FoldChange'] < -0.5) & (wdf['pvalue'] < 0.05)]
        ax.scatter(wdf_filt1['log2FoldChange'].tolist(), wdf_filt1['Mod_pvalue'].tolist(),
                    s=200, c="blue")

        wdf_filt2 = wdf[(wdf['log2FoldChange'] > 0.5) & (wdf['pvalue'] < 0.05)]
        ax.scatter(wdf_filt2['log2FoldChange'].tolist(), wdf_filt2['Mod_pvalue'].tolist(),
                    s=200, c="red")

        # With the subsets, calculate and plot percentages of genes up/down
        perc_down = (len(wdf_filt1)/len(wdf))*100
        ax.text(xleftlim*0.9, ytoplim*1.35, f"{perc_down:.1f}%", fontsize=30, c="blue")

        perc_up = (len(wdf_filt2)/len(wdf))*100
        ax.text(xrightlim*0.4, ytoplim*1.35, f"{perc_up:.1f}%", fontsize=30, c="red")

        # Set axis limits
        ax.set_ylim(bottom=0, top=ytoplim*1.5)
        ax.set_xlim(xleftlim, xrightlim)

        # Set axis labels
        ax.set_xlabel("Fold Change (log2FoldChange)", {'fontsize': 20})
        ax.set_ylabel("P value (-log10(pvalue))", {'fontsize': 20})

        # Set plot title
        ax.set_title(title.replace('GENEMARKER', glist), {'fontsize': 25})

        fig.savefig(out_path.replace('GENEMARKER', glist))
        plt.close()

def volcano_markers(config, tool_name):
    """Get the scatter plots for the markers
    """

    # Start the log for the process
    logging.info(f'Starting {tool_name} process')

    # Extract the organism
    organism = config['options']['organism']

    # Extract the comparisons
    comparisons = config['comparisons']

    # Extract the infiles and the project name
    in_path = "/".join(config['tools_conf'][tool_name]['input']['DEtouched'].split('/')[0:-1])
    project = config['project']

    # Extract the out path and out markerfile
    outmarker = config['tools_conf'][tool_name]['output']['MVtouched']
    out_path = "/".join(outmarker.split('/')[0:-1])

    # Select the list of markers
    gene_markers = get_gene_markers(organism, out_path)

    # Get the dimension char
    dimensions = config['tools_conf'][tool_name]['tool_conf']['dimensions']
    dimensions = dimensions.split(',')

    # Create out folder and marker
    command = ""
    if not os.path.exists(out_path):
        command += f"mkdir {out_path};"

    command += f"touch {outmarker}; "

    pf.run_command(command)

    # Loop through the samples and controls
    for control in comparisons:
        samples = comparisons[control].split(",")
        for sample in samples:
            # Get the paths for in table, out plot format. Use wildcard /GENENAME/ in the name of the plot
            tab_name = f"{in_path}/{project}_{sample}_{control}.tsv"
            out_path = f'{out_path}/GENEMARKER_{sample}_{control}_scattermarkers.png'

            # Establish the plot title
            plot_title = f"GENEMARKER - {sample} - {control}"

            # Run the plots
            volcano_marker_plot(tab_name, gene_markers, out_path, plot_title, dimensions)