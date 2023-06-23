import logging
import os
import glob
import python_scripts.python_functions as pf

def fastqc(config, tool_name):
    # Report the start in the log
    logging.info(f'Starting {tool_name} process')

    # Get the paths to the files that are going to be input/output
    # Inputs
    # path to (some) fastq files
    fastq_r1 = config['tools_conf'][tool_name]['input']['fastq_r1']
    # Gett all the fastq files
    fastq_path = "/".join(fastq_r1[0].split('/')[0:-1])
    fastq_files = glob.glob(f'{fastq_path}/*.fastq.gz')

    # Outputs
    # Control file
    fastqctouched = config['tools_conf'][tool_name]['output']['fastqctouched']
    # Out folder
    fastqcfolder = "/".join(fastqctouched.split('/')[0:-1])

    # Other information
    threads = config['tools_conf'][tool_name]['tool_conf']['threads']

    # Init the command
    command = ''

    # Create the out folder
    if not os.path.exists(fastqcfolder):
        command += f"mkdir {fastqcfolder}; "

    # Make the fastqc command
    command += f'fastqc -o {fastqcfolder} -t {threads}'

    for fastq in fastq_files:
        command += f' {fastq}'
    command += '; '
    
    command += f"touch {fastqctouched}; "
    pf.run_command(command)


    