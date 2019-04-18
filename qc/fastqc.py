import argparse
import os
import logging
import python_functions as pf

def fastqc(config, tool_name, logger):
    """Method to build local parameters for the tool to work

    The method sets the local parameters of the class that are going to be
    needed for the process

    """

    # Get the paths to the files that are going to be input/output
    R1_FILES = config['tools_conf'][tool_name]['input']['fastq_r1']
    R2_FILES = config['tools_conf'][tool_name]['input']['fastq_r2']
    fastqcdir = config['tools_conf'][tool_name]['output']['fastqcdir']
    threads = config['tools_conf'][tool_name]['tool_conf']['threads']

    # Add the fastq files in a single list
    FQfiles = R1_FILES + R2_FILES

    # Construct the command
    full_cmd = ""

    if not os.path.exists(fastqcdir):
        full_cmd += f"mkdir {fastqcdir};"

    full_cmd += f"parallel -j {threads} fastqc -o {fastqcdir} {{}} :::: <(ls {' '.join(FQfiles)})"
    full_cmd += f"touch {fastqcdir}/fastqc.done.txt"

    logger.info(command)

    pf.run_command(command, logger)


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
    parser.add_argument('--fastqcdir', required=True, help='Folder for fastqc files')

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
      "outfolder": args.outpath,
      "log_files": ["/tmp/full.log"],
      "tools_conf": {
        "fastqc": {
          "input": {
            "fastq_r1": args.fastq_r1.split(','),
            "fastq_r2": args.fastq_r2.split(',')
            },
          "output": {
            "fastqcdir": args.fastqcdir
            },
          "tool_conf": {
            "threads": "2"
            }
          }
        }
      }

    # Startup the logger format
    logger = pf.create_logger(config['log_files'][0])

    fastqc(config, 'fastqc', logger)


if __name__ == "__main__":
    main()
