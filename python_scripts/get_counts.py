import argparse
import logging
import pandas as pd
import python_scripts.python_functions as pf

def transformToRanges(counts):
    """
    Function to transform the segment of featureCounts to ranges for R
    """
    # Read the file
    df = pd.read_csv(counts, sep='\t', index_col=None)

    # Transform the pseudo-dataframe (cells with more than one field of data) into a normal dataframe.
    new_list = []
    for index, row in df.iterrows():
        chr_list = row.Chr.split(";")
        start_list = row.Start.split(";")
        end_list = row.End.split(";")
        strand_list = row.Strand.split(";")

        for i in range(0,len(chr_list)):
            newline = [row.Geneid, chr_list[i], start_list[i], end_list[i], strand_list[i]]
            new_list.append(newline)

    df2 = pd.DataFrame(new_list)
    df2.columns = ['Geneid','Chr','Start','End','Strand']
    out_counts = counts.replace(".tmpranges.tsv",".ranges.tsv")
    df2.to_csv(out_counts, sep='\t', index=False)

def fixFormat(counts):
    """
    Load the counts table and remove duplicated gene names
    """
    df = pd.read_csv(counts, sep='\t', index_col=None)
    df = df.drop_duplicates(subset='Geneid')
    df.to_csv(counts, sep='\t', index=False)

def counts(config, tool_name):
    """Get the counts of a number of bam files in a directory
    """
    logging.info(f'Starting {tool_name} process')

    # Get information for the scripts
    # Input
    bamdir = "/".join(config['tools_conf'][tool_name]['input']['bamdir'].split('/')[0:-1])

    # Output
    output = config['tools_conf'][tool_name]['output']['counts']

    # Other
    annot = config['tools_conf'][tool_name]['input']['annot']

    # NOTE: On hold
    # rangestable = output.replace(".tsv",".tmpranges.tsv")

    # List all the bam files in the directory
    filelist = pf.list_files_dir(bamdir, ext = '*.bam')

    # Create the featurecounts command
    tmpoutput = output.replace('.tsv','.tmp.tsv')
    command = f'featureCounts -p -a {annot} -o {tmpoutput} {" ".join(filelist)}; '

    # Reformat the output of featurecounts into a readable table
    command += f"cat {tmpoutput} | tail -n +2 | sed -r 's/\t([^\t]+)\//\t/g' | sed 's/.bam//g' | cut --complement -f 2,3,4,5,6 | perl -pe 's|(\.).*?\t|\t|' > {output};"

    # NOTE: On hold
    # command += f"cat {tmpoutput} | tail -n +2 | sed -r 's/\t([^\t]+)\//\t/g' | sed 's/.bam//g' | cut -f 1,2,3,4,5 > {rangestable}"

    # Run the command(s)
    pf.run_command(command)

    # Human annotation has an annoying detail (some gene names and ensemblIDs are duplicated and scripts downstream hate that) and needs to be corrected
    if config['options']['organism'] == "human":
        fixFormat(output)

    # Create the ranges file. Useful later for the fpkm
    # NOTE: I don't need this at the moment, but I have the feeling that I will need this for other stuff in the future,
    # so I'm gonna leave this here, but commented, along with the other chunks of code neccesary for this.
    # transformToRanges(rangestable)