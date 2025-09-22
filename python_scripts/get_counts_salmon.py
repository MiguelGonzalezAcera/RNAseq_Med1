import logging
import os
import pandas as pd
import glob
import python_scripts.python_functions as pf

def counts_sal(config, tool_name):
    """Get the counts of a number of bam files in a directory
    """
    logging.info(f'Starting {tool_name} process')

    # Get information for the scripts
    # Input
    bamdir = "/".join(config['tools_conf'][tool_name]['input']['bamdir'].split('/')[0:-1])

    # List all the bam files in the directory
    filelist = pf.list_files_dir(bamdir, ext = '*.unsorted.bam')

    # Output
    # Control file
    counts_sal_touched = config['tools_conf'][tool_name]['output']['counts_sal_touched']
    # Out folder
    salmondir = "/".join(counts_sal_touched.split('/')[0:-1])

    # Other
    annot = config['tools_conf'][tool_name]['input']['annot']
    gentr = config['tools_conf'][tool_name]['input']['gentr']

    # Init the command
    command = ""

    # Make the bam directory if it does not exist
    if not os.path.exists(salmondir):
        command += f"mkdir {salmondir}; "

    g_uns_path = f"{salmondir}/genes_unsorted"

    if not os.path.exists(g_uns_path):
        command += f"mkdir {g_uns_path}; "

    # Create the salmon command
    for file in filelist:
        # Some names
        file_sf = file.replace('.unsorted.bam','.sf')
        file_genes_sf = file.replace('.unsorted.bam','.genes.sf')

        command += f"salmon quant --geneMap {annot} -p 10 --libType A -t {gentr} -a {file} -o {salmondir} -q; mv -v {salmondir}/quant.sf {file_sf}; mv -v {salmondir}/quant.genes.sf {file_genes_sf}; mv -v {bamdir}/*.sf {salmondir}; mv -v {salmondir}/*.genes.sf {salmondir}/genes_unsorted; "

    # Touch the marker file
    command += f"touch {counts_sal_touched}"

    # Run the command(s)
    pf.run_command(command)

    #------------------------------------------------------------------

    # Get the counts by transcript for future analysis using isoforms
    salmon_files = glob.glob(f'{salmondir}/*.sf')

    resdf_counts = pd.DataFrame()
    resdf_tpm = pd.DataFrame()

    for file in salmon_files:
        filename = file.split("/")[-1].replace(".sf","")
        
        df = pd.read_csv(file, sep='\t')

        df_reads = df[['Name','NumReads']]
        df_reads.columns = ['TransID',filename]

        if resdf_counts.empty:
            resdf_counts = df_reads
        else:
            resdf_counts = pd.merge(resdf_counts, df_reads, how='outer', on='TransID')

        df_tpm = df[['Name','TPM']]
        df_tpm.columns = ['TransID',filename]

        if resdf_tpm.empty:
            resdf_tpm = df_tpm
        else:
            resdf_tpm = pd.merge(resdf_tpm, df_tpm, how='outer', on='TransID')

    resdf_counts.to_csv(f'{salmondir}/transcript_counts.tsv', sep='\t', index=False)
    resdf_tpm.to_csv(f'{salmondir}/transcript_countsTPM.tsv', sep='\t', index=False)
