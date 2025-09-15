import logging
import os
import pandas as pd
import python_scripts.python_functions as pf

def counts_sal(config, tool_name):
    """Get the counts of a number of bam files in a directory
    """
    logging.info(f'Starting {tool_name} process')

    # Get information for the scripts
    # Input
    bamdir = "/".join(config['tools_conf'][tool_name]['input']['bamdir'].split('/')[0:-1])

    # Output
    # Control file
    counts_sal_touched = config['tools_conf'][tool_name]['output']['counts_sal_touched']
    # Out folder
    salmondir = "/".join(counts_sal_touched.split('/')[0:-1])

    # Other
    annot = config['tools_conf'][tool_name]['input']['annot']
    gentr = config['tools_conf'][tool_name]['input']['gentr']

    # Make the bam directory if it does not exist
    if not os.path.exists(salmondir):
        command += f"mkdir {salmondir}; "

    # Create the salmon command
    command = f"for file in {bamdir}/*unsorted.bam; do echo $file; salmon quant --geneMap {annot} --libType A -t {gentr} -a $file -o {salmondir} -q; mv -v {salmondir}/quant.sf ${{file%.unsorted.bam}}.sf; mv -v {salmondir}/quant.genes.sf ${{file%.unsorted.bam}}.genes.sf; done; mv -v {bamdir}/*.sf {salmondir}; mkdir {salmondir}/genes_unsorted; mv -v {salmondir}/*.genes.sf {salmondir}/genes_unsorted;"

    # Touch the marker file
    command += f"touch {counts_sal_touched}"

    # Run the command(s)
    pf.run_command(command)
