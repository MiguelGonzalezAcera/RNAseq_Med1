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
    logging.info(f'Starting {tool_name} process')

    # Get the paths to the files that are going to be input/output
    R1_FILES = config['tools_conf'][tool_name]['input']['fastq_r1']
    if 'fastq_r2' in config['tools_conf'][tool_name]['input']:
        R2_FILES = config['tools_conf'][tool_name]['input']['fastq_r2']
    else:
        R2_FILES = ''
    bamdir = "/".join(config['tools_conf'][tool_name]['output']['mappingtouched'].split('/')[0:-1])
    mappingtouched = config['tools_conf'][tool_name]['output']['mappingtouched']
    genomePath = config['tools_conf']['genome']
    threads = config['tools_conf'][tool_name]['tool_conf']['threads']

    ## Possible Command style
    # STAR --genomeLoad LoadAndExit --genomeDir /DATA/references/star_genomes/mmu38/star_indices_overhang150/; for file in $(cat fastq.test.txt); do echo $file STAR --runThreadN 10 --readFilesCommand gzip -cd --genomeDir /DATA/references/star_genomes/mmu38/star_indices_overhang150/ --readFilesIn $file ${file%_1.fastq.gz}_2.fastq.gz --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate > ${file%_1.fastq.gz}.bam; samtools-1.9 index ${file%_1.fastq.gz}.bam; mv ${file%_1.fastq.gz}.bam* BAM/; done; STAR --genomeLoad Remove --genomeDir /DATA/references/star_genomes/mmu38/star_indices_overhang150/
    command = ""

    command += f"STAR --genomeLoad LoadAndExit --genomeDir {genomePath}; "
    if not os.path.exists(bamdir):
        command += f"mkdir {bamdir}; "
    for filer1 in R1_FILES:
        if R2_FILES:
            bamfile = bamdir + "/" + filer1.split("/")[-1].replace('_1.fastq.gz','.bam')
            filer2 = filer1.replace('_1.fastq.gz','_2.fastq.gz')
            command += f'STAR --runThreadN {threads} --readFilesCommand gzip -cd --genomeDir {genomePath} --readFilesIn {filer1} {filer2} --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate > {bamfile}; samtools index {bamfile}; '
        else:
            bamfile = bamdir + "/" + filer1.split("/")[-1].replace('.fastq.gz','.bam')
            command += f'STAR --runThreadN {threads} --readFilesCommand gzip -cd --genomeDir {genomePath} --readFilesIn {filer1} --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate > {bamfile}; samtools index {bamfile}; '
    command += f"STAR --genomeLoad Remove --genomeDir {genomePath}; "
    command += f"touch {mappingtouched}"

    pf.run_command(command)

def get_arguments():
    """
    Function that parse arguments given by the user, returning a dictionary
    that contains all the values.
    """

    # Create the top-level parser
    parser = argparse.ArgumentParser()

    # Mandatory variables
    parser.add_argument('--fastq_r1', '-f1', nargs='*', required=True, help='fastq_r1 file/s')
    parser.add_argument('--fastq_r2', '-f2', nargs='*', help='fastq_r2 file/s. Only when paired end')
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

    # Select if the run is going to have two sets of fastq or one
    if args.run_type == 'single':
      fastq_dict = {
        "fastq_r1": args.fastq_r1.split(',')
      }
    else:
      fastq_dict = {
        "fastq_r1": args.fastq_r1.split(','),
        "fastq_r2": args.fastq_r2.split(',')
      }

    config = {
      "DEBUG": args.debug,
      "TESTING": args.test,
      "DRY_RUN": args.dry_run,
      "log_files": ["/tmp/full.log"],
      "tools_conf": {
        "mapping": {
          "input": fastq_dict,
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
