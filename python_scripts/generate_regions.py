import argparse
import logging
import pandas as pd
import python_scripts.python_functions as pf


def regions(config, tool_name):
    """Get the counts of a number of bam files in a directory
    """

    bedfile = config['tools_conf'][tool_name]['input']['bedfile']
    genome = config['tools_conf'][tool_name]['tool_conf']['genome']
    list = config['tools_conf'][tool_name]['output']['list']

    # Create the picard command
    command = ""
    command += f"java -jar /SOFTWARE/bin/picard-2.19.0.jar BedToIntervalList I={bedfile}" + \
           f" O={list} SD={genome}"

    pf.run_command(command)

def get_arguments():
    """
    Function that parse arguments given by the user, returning a dictionary
    that contains all the values.
    """

    # Create the top-level parser
    parser = argparse.ArgumentParser()

    # Mandatory variables
    parser.add_argument('--bedfile', required=True, help='Regions file in bed format')
    parser.add_argument('--genome_dict', required=True, help='Genome dictionary file to generate the proper Regions')
    parser.add_argument('--list', required=True, help='File generated for the mapping. Used later in picard')

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
        "generate_regions": {
          "input": {
            "bedfile": args.bedfile
            },
          "output": {
            "list": args.list
            },
          "tool_conf": {
            "genome": args.genome_dict
            }
          }
        }
      }

    # Startup the logger format

    regions(config, 'get_counts')


if __name__ == "__main__":
    main()
