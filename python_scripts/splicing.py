import logging
import os
import pandas as pd
import python_scripts.python_functions as pf

def splicing(config, tool_name):
    # Report the start in the log
    logging.info(f'Starting {tool_name} process')

    # Get the paths to the files that are going to be input/output
    # Inputs
    # bam file of files
    bamfof = config['tools_conf'][tool_name]['input']['bamfof']
    # Path from the bamfof, for purposes
    bamfofpath = "/".join(bamfof.split('/')[0:-1])
    # design file
    design = config['tools_conf'][tool_name]['input']['design']

    # Outputs
    # Control file
    splicetouched = config['tools_conf'][tool_name]['output']['splicetouched']
    # Out folder
    splicedir = "/".join(splicetouched.split('/')[0:-1])

    # Other information
    annot = config['tools_conf'][tool_name]['input']['annot']
    threads = config['tools_conf'][tool_name]['tool_conf']['threads']
    # comparison set
    samples = config['comparisons']
    # project name
    project = config['project']
    # organism used
    organism = config['options']['organism']

    # Read the fof to make them into a string to pass to the command
    with open(bamfof) as file:
        bamfoflines = [line.rstrip() for line in file]
    
    bamfiles = ",".join(bamfoflines)

    # Part 1. Run the build step for the experiment (steps from https://spladder.readthedocs.io/en/latest/spladder_cohort.html)

    # Create the initial command object
    command = ""

    # Step 1 (parallel). Single graphs.
    command += f"parallel --joblog splicing_logfile.log -j {threads} \"spladder build -o {splicedir} -a {annot} -b {{}} --merge-strat single --no-extract-ase --ignore-mismatches\" :::: {bamfof}; "

    # Step 2 (single. use bam fof). Merged graph.
    command += f"spladder build -o {splicedir} -a {annot} -b {bamfiles} --merge-strat merge_graphs --no-extract-ase --ignore-mismatches; "

    # Step 3 (first parallel second single). Quantification.
    command += f"parallel --joblog splicing_logfile.log -j {threads} \"spladder build -o {splicedir} -a {annot} -b {{}} --merge-strat merge_graphs --no-extract-ase --quantify-graph --qmode single --ignore-mismatches\" :::: {bamfof}; "

    command += f"spladder build -o {splicedir} -a {annot} -b {bamfiles} --merge-strat merge_graphs --no-extract-ase --quantify-graph --qmode collect --ignore-mismatches; "

    # Step 4 (single). Calling events.
    command += f"spladder build -o {splicedir} -a {annot} -b {bamfiles}; "

    # Run the commands
    pf.run_command(command)

    # Part 2. Run the differential splicing tests.

    # Get the design to find the sample names associated to each condition
    designdf = pd.read_csv(design, index_col=None, sep = '\t')
    designdf.columns = ['Sample_name', 'Condition']

    # Iter through the comparisons
    for control in samples:
        sample_ids = samples[control].split(",")
        for sample in sample_ids:
            # Restart the command
            command = ""

            # Get the bam files from the control samples
            controlbam = ",".join([f"{bamfofpath}/{bam}.bam" for bam in designdf[designdf['Condition'] == control]['Sample_name'].tolist()])

            # Get the bam files from the other samples
            samplebam = ",".join([f"{bamfofpath}/{bam}.bam" for bam in designdf[designdf['Condition'] == sample]['Sample_name'].tolist()])

            # Make the command for testing
            command = f"spladder test --conditionA {controlbam} --conditionB {samplebam} --outdir {splicedir}; "

            # Run the command
            pf.run_command(command)

            # Rename the resulting folder
            os.rename(f"{splicedir}/testing", f"{splicedir}/{sample}_{control}_testing")

    # Last part: Touch the marker file for ending

    command = f"touch {splicetouched}; "
    pf.run_command(command)


    