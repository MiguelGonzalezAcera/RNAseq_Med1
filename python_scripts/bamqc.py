import logging
import os
import python_scripts.python_functions as pf

def bamqc(config, tool_name):
    # Report the start in the log
    logging.info(f'Starting {tool_name} process')

    # Get the paths to the files that are going to be input/output
    # Inputs
    # path to the bamfof file
    bamfof = config['tools_conf'][tool_name]['input']['bamfof']

    # Outputs
    # Control file
    bamqctouched = config['tools_conf'][tool_name]['output']['bamqctouched']
    # Out folder
    bamqcfolder = "/".join(bamqctouched.split('/')[0:-1])

    # Read the fof to make them into a string to pass to the command
    with open(bamfof) as file:
        bamfoflines = [line.rstrip() for line in file]

    # Init the command
    command = ''

    # Create the out folder
    if not os.path.exists(bamqcfolder):
        command += f"mkdir {bamqcfolder}; "

    # iter through the bamfiles to generate the commands
    for bamfile in bamfoflines:
        # Get name of the bam file
        bamname = bamfile.split('/')[-1]

        # Generate the picard command
        command += f"java -jar $PICARD CollectAlignmentSummaryMetrics VALIDATION_STRINGENCY=SILENT  I= {bamfile} O= {bamqcfolder}/{bamname}qc; "
    
    # Generate the touched command
    command += f"touch {bamqctouched}; "
    pf.run_command(command)
