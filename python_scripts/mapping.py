"""

This script functions as a main program to map sequences from a given set of
fastq files (R1 and R2).

The user can select which the output would be and which one would be the
reference genome to use (defaults to mm10).

"""
import argparse
import logging
import os
import python_functions as pf

def mapping(config, tool_name, logger):
    """Method to build local parameters for the tool to work

    The method sets the local parameters of the class that are going to be
    needed for the process

    """

    # Get the paths to the files that are going to be input/output
    R1_FILES = config['tools_conf'][tool_name]['input']['fastq_r1']
    R2_FILES = config['tools_conf'][tool_name]['input']['fastq_r2']
    bamdir = config['tools_conf'][tool_name]['output']['bam_dir']
    genomePath = config['tools_conf'][tool_name]['tool_conf']['genome']
    threads = config['tools_conf'][tool_name]['tool_conf']['threads']
    path = "/".join(R1_FILES[0].split("/")[:-1])
    print(path)

    ## Possible Command style
    # STAR --genomeLoad LoadAndExit --genomeDir /DATA/references/star_genomes/mmu38/star_indices_overhang150/; for file in $(cat fastq.test.txt); do echo $file STAR --runThreadN 10 --readFilesCommand gzip -cd --genomeDir /DATA/references/star_genomes/mmu38/star_indices_overhang150/ --readFilesIn $file ${file%_1.fastq.gz}_2.fastq.gz --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate > ${file%_1.fastq.gz}.bam; samtools-1.9 index ${file%_1.fastq.gz}.bam; mv ${file%_1.fastq.gz}.bam* BAM/; done; STAR --genomeLoad Remove --genomeDir /DATA/references/star_genomes/mmu38/star_indices_overhang150/
    command = ""

    command += f"STAR --genomeLoad LoadAndExit --genomeDir {genomePath}; "
    command += f'for file in {path + "/*1.fastq.gz"}; do STAR --runThreadN {threads} --readFilesCommand gzip -cd \
    --genomeDir {genomePath} --readFilesIn $file ${{file%_1.fastq.gz}}_2.fastq.gz --outSAMtype BAM SortedByCoordinate --outStd \
    BAM_SortedByCoordinate > ${{file%_1.fastq.gz}}.bam; samtools-1.9 index ${{file%_1.fastq.gz}}.bam; done; '
    command += f"STAR --genomeLoad Remove --genomeDir {genomePath}; "
    if not os.path.exists(bamdir):
        command += f"mkdir {bamdir}; "
    command += f'mv {path + "/*.bam*"} {bamdir}; '

    logger.info(command)

    pf.run_command(command, logger)
    print(command)

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
            "fastq_r1": args.fastq_r1,
            "fastq_r2": args.fastq_r2
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
    logger = pf.create_logger(config['log_files'][0])

    mapping(config, 'mapping', logger)


if __name__ == "__main__":
    main()
