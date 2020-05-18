import argparse
import os
import python_scripts.python_functions as pf

def coverage_per_base(config, tool_name):
    """Get the counts of a number of bam files in a directory
    """

    logging.info(f'Starting {tool_name} process')

    bamdir = "/".join(config['tools_conf'][tool_name]['input']['bamdir'].split('/')[0:-1])
    bed = config['tools_conf'][tool_name]['input']['bed']
    cpbdir = "/".join(config['tools_conf'][tool_name]['output']['cpbtouched'].split('/')[0:-1])
    cpbtouched = config['tools_conf'][tool_name]['output']['cpbtouched']
    genome = config['tools_conf']['genometxt']

    # Lsit all the bam files in the directory
    filelist = pf.list_files_dir(bamdir, ext = '*.bam')

    # Get path
    path = "/".join(filelist[0].split("/")[0:-1])

    # Save the filenames under a single file, for use in parallel fashion
    cpb_file = open(f'{path}/bam.fof', 'w')
    for cpb in filelist:
        cpb_file.write(f"{cpb}\n")
    cpb_file.close()

    # Create the cpb command
    command = ""
    # Option 1: sequential command
    # for file in filelist
    #     command += f"bedtools coverage -g {genome} -sorted -d -a {annot} -b {file} > {file.replace('.bam','.cpb')};"

    # Option 2: Parallel (make sure that we have enough ram memory)
    command += f'parallel --joblog split_logfile.log -j 10 "bedtools coverage -g {genome} -sorted -d -a {bed} -b {{}} > {{= s:\.bam:\.cpb: =}}" :::: {path}/bam.fof; '

    if not os.path.exists(cpbdir):
        command += f"mkdir {cpbdir}; "
    command += f'mv {path}/*.cpb {cpbdir}; '

    command += f'touch {cpbtouched}'

    pf.run_command(command)

def get_arguments():
    """
    Function that parse arguments given by the user, returning a dictionary
    that contains all the values.
    """

    # Create the top-level parser
    parser = argparse.ArgumentParser()

    # Mandatory variables
    parser.add_argument('--bamdir', required=True, help='Folder with bam files')
    parser.add_argument('--bed', required=True, help='Bedfile with the regiosn dor the cpb')
    parser.add_argument('--genome', required=True, help='Genome to assist with the process')
    parser.add_argument('--cpbdir', required=True, help='Folder for the cppb files')

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
        "coverage_per_base": {
          "input": {
            "bamdir": args.bamdir,
            "bed": args.bed
            },
          "output": {
            "cpbdir": args.cpbdir
            },
          "tool_conf": {
            "genome": args.genome
            }
          }
        }
      }

    coverage_per_base(config, 'coverage_per_base')


if __name__ == "__main__":
    main()
