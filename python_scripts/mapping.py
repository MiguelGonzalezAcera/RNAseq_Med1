"""
This script functions as a main program to map sequences from a given set of
fastq files (R1 and R2).

The user can select which the output would be and which one would be the
reference genome to use (defaults to mm10).
"""
import argparse
import logging
import os
import python_scripts.python_functions as pf

def mapping(config, tool_name):
    """Method to build local parameters for the tool to work

    The method sets the local parameters of the class that are going to be
    needed for the process

    """
    # Start the logger
    logging.info(f'Starting {tool_name} process')

    # Get the paths to the files that are going to be input/output
    # Inputs
    # R1 files
    R1_FILES = config['tools_conf'][tool_name]['input']['fastq_r1']

    # Outputs
    # Control file
    mappingtouched = config['tools_conf'][tool_name]['output']['mappingtouched']
    # Out folder
    bamdir = "/".join(mappingtouched.split('/')[0:-1])

    # Other information
    genomePath = config['tools_conf']['genome']
    threads = config['tools_conf'][tool_name]['tool_conf']['threads']

    # Init the command
    command = ""

    # Pre-load the genome, as it would save time and processing power for each sample
    command += f"STAR --genomeLoad LoadAndExit --genomeDir {genomePath}; "

    # Make the bam directory if it does not exist
    if not os.path.exists(bamdir):
        command += f"mkdir {bamdir}; "
    
    # Iter through the fastq filenames
    for filer1 in R1_FILES:
        # Make a different command when the run is with paired or single end
        if config['options']['reads'] == 'paired':
            # Make the name of the bam file from the fastq file
            bamfile = bamdir + "/" + filer1.split("/")[-1].replace('_1.fastq.gz','.bam')
            # Replace extension for the R2 file
            filer2 = filer1.replace('_1.fastq.gz','_2.fastq.gz')

            # Make the star mapper command, with the samtools indexing of the bam file
            command += f'STAR --runThreadN {threads} --readFilesCommand gzip -cd --genomeDir {genomePath} --readFilesIn {filer1} {filer2} --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate > {bamfile}; samtools index {bamfile}; '
        else:
            # Make the name of the bam file from the fastq file
            bamfile = bamdir + "/" + filer1.split("/")[-1].replace('.fastq.gz','.bam')

            # Make the star mapper command, with the samtools indexing of the bam file
            command += f'STAR --runThreadN {threads} --readFilesCommand gzip -cd --genomeDir {genomePath} --readFilesIn {filer1} --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate > {bamfile}; samtools index {bamfile}; '
    
    # Remove the loaded genome from memory
    command += f"STAR --genomeLoad Remove --genomeDir {genomePath}; "

    # Create tracking file
    command += f"touch {mappingtouched}"

    # Run the commands
    pf.run_command(command)

def get_arguments():
    """
    Function that parse arguments given by the user, returning a dictionary
    that contains all the values.
    """

    # Create the top-level parser
    parser = argparse.ArgumentParser()

    # Mandatory variables
    parser.add_argument('--fastq', '-fq', nargs='*', required=True, help='fastq file/s, comma separated, no space. If the run is single end, extension should be .fastq.gz. If the run is paired, only the R1 files should be stated, with extension _1.fastq.gz, and the R2 files should be named the same with extension _2.fastq.gz.')
    parser.add_argument('--bamdir', required=True, help='Folder for bam files')

    # Default parameters
    parser.add_argument('--genome', default='/DATA/references/star_genomes/mmu39/star_indices_overhang150/', help='Genome folder for STAR mapping')
    parser.add_argument('--run_type', default='paired', choices = ['paired', 'single'], help='Specify if the run is done with single end or paired end reads. Defaults to paired.')

    # Test and debug variables
    parser.add_argument('--dry_run', action='store_true', default=False, help='debug')
    parser.add_argument('--debug', '-d', action='store_true', default=False, help='dry_run')
    parser.add_argument('--test', '-t', action='store_true', default=False, help='test')

    # parse some argument lists
    args = parser.parse_args()

    return args


def main():
    """
    Main function of the script. Launches the rest of the process
    """

    # Get arguments from user input
    args = get_arguments()

    # Make the config dictionary
    config = {
      "DEBUG": args.debug,
      "TESTING": args.test,
      "DRY_RUN": args.dry_run,
      "log_files": ["/tmp/full.log"],
      "tools_conf": {
        "mapping": {
          "input": {
            "fastq_r1": args.fastq.split(',')
          },
          "output": {
            "mappingtouched": args.bamdir + "/mappingtouched.txt"
            },
          "tool_conf": {
            "threads": "2"
            }
          },
        'genome': args.genome,
        },
      'options': {
        'reads': args.run_type
      }
      }

    # Startup the logger format

    mapping(config, 'mapping')


if __name__ == "__main__":
    main()
