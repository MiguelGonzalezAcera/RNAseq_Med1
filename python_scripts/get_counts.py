import argparse
import logging
import python_functions as pf

def counts(config, tool_name, logger):
    """Get the counts of a number of bam files in a directory
    """

    bamdir = config['tools_conf'][tool_name]['input']['bamdir']
    annot = config['tools_conf'][tool_name]['input']['annot']
    output = config['tools_conf'][tool_name]['output']['counts']

    # Lsit all the bam files in the directory
    filelist = pf.list_files_dir(bamdir, ext = '.bam')

    # Create the featurecounts command
    command = f"featureCounts -a {annot} -o {output} {" ".join(filelist)}"
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
    parser.add_argument('--bamdir', nargs='*', required=True, help='Folder with bam files')
    parser.add_argument('--counts', required=True, help='Table with the counts')
    parser.add_argument('--annot', required=True, help='Annotation file (same than used in mapping)')

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
        "counts": {
          "input": {
            "bamdir": args.bamdir,
            "annot": args.annot
            },
          "output": {
            "counts": args.counts
            },
          "tool_conf": {
            }
          }
        }
      }

    # Startup the logger format
    logger = pf.create_logger(config['log_files'][0])

    counts(config, 'get_counts', logger)


if __name__ == "__main__":
    main()
