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

    # Get the paths to the files that are going to be input/output
    R1_FILES = config['tools_conf'][tool_name]['input']['fastq_r1']
    R2_FILES = config['tools_conf'][tool_name]['input']['fastq_r2']
    bamdir = "/".join(config['tools_conf'][tool_name]['output']['mappingtouched'].split('/')[0:-1])
    mappingtouched = config['tools_conf'][tool_name]['output']['mappingtouched']
    genomePath = config['tools_conf'][tool_name]['tool_conf']['genome']
    threads = config['tools_conf'][tool_name]['tool_conf']['threads']

    ## Possible Command style
    # STAR --genomeLoad LoadAndExit --genomeDir /DATA/references/star_genomes/mmu38/star_indices_overhang150/; for file in $(cat fastq.test.txt); do echo $file STAR --runThreadN 10 --readFilesCommand gzip -cd --genomeDir /DATA/references/star_genomes/mmu38/star_indices_overhang150/ --readFilesIn $file ${file%_1.fastq.gz}_2.fastq.gz --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate > ${file%_1.fastq.gz}.bam; samtools-1.9 index ${file%_1.fastq.gz}.bam; mv ${file%_1.fastq.gz}.bam* BAM/; done; STAR --genomeLoad Remove --genomeDir /DATA/references/star_genomes/mmu38/star_indices_overhang150/
    command = ""

    command += f"STAR --genomeLoad LoadAndExit --genomeDir {genomePath}; "
    if not os.path.exists(bamdir):
        command += f"mkdir {bamdir}; "
    for filer1 in R1_FILES:
        filer2 = filer1.replace('_1.fastq.gz','_2.fastq.gz')
        bamfile = bamdir + "/" + filer1.split("/")[-1].replace('_1.fastq.gz','.bam')
        command += f'STAR --runThreadN {threads} --readFilesCommand gzip -cd --genomeDir {genomePath} --readFilesIn {filer1} {filer2} --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate > {bamfile}; samtools-1.9 index {bamfile}; '
    command += f"STAR --genomeLoad Remove --genomeDir {genomePath}; "
    command += f"touch {mappingtouched}"

    print(command)
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
    parser.add_argument('--fastq_r2', '-f2', nargs='*', required=True, help='fastq_r2 file/s')
    parser.add_argument('--bamdir', required=True, help='Folder for bam files')

    # Default parameters
    parser.add_argument('--genome', default='/DATA/references/star_genomes/mmu38/star_indices_overhang150/', help='Genome folder for STAR mapping')

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

    config = {
      "DEBUG": args.debug,
      "TESTING": args.test,
      "DRY_RUN": args.dry_run,
      "log_files": ["/tmp/full.log"],
      "tools_conf": {
        "mapping": {
          "input": {
            "fastq_r1": args.fastq_r1.split(','),
            "fastq_r2": args.fastq_r2.split(','),
            "samplelist": a
            },
          "output": {
            "bam_dir": args.bamdir
            },
          "tool_conf": {
            "genome": args.genome,
            "threads": "2"
            }
          }
        }
      }

    # Startup the logger format

    mapping(config, 'mapping')


if __name__ == "__main__":
    main()
