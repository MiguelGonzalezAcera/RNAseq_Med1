import argparse
import logging
import os
import pandas as pd
import python_scripts.python_functions as pf


def GSVA(config, tool_name):
    """Get the
    """
    # Start log
    logging.info(f'Starting {tool_name} process')

    # Get input parameters
    norm_counts = config['tools_conf'][tool_name]['input']['norm_counts']
    design = config['tools_conf'][tool_name]['input']['design']

    # define the out directory
    out_dir = "/".join(config['tools_conf'][tool_name]['output']['heatmap'].split('/')[0:-1])

    # Define the rest of the output parameters
    heatmap = config['tools_conf'][tool_name]['output']['heatmap']

    # Get options
    organism = config['options']['organism']
    project = config['project']
    samples = config['comparisons']

    # Create the command to run the pca R script
    command = ""

    # Make the outfolder if it doesnt wxist
    if not os.path.exists(out_dir):
        command += f"mkdir {out_dir};"

    # Loop through controls
    for control in samples:
        # Add the main command
        command += f'Rscript Rscripts/GSVA.r --heatmap {heatmap} --counts {norm_counts} --design {design} --organism {organism} --out_obj {out_dir}/{project}_GSVA.Rda --control {control} --comparisons {samples[control]}; '

    # Run the command and log it
    pf.run_command(command)


def get_arguments():
    """
    Function that parse arguments given by the user, returning a dictionary
    that contains all the values.
    """

    # Create the top-level parser
    parser = argparse.ArgumentParser()

    # Mandatory variables
    parser.add_argument('--counts', required=True, help='Table with the counts of the assay')
    parser.add_argument('--design', required=True, help='Table with the design of the experiment')
    parser.add_argument('--project', required=True, help='Project name')
    parser.add_argument('--heatmap', required=True, help='Heatmap to save')
    parser.add_argument('--organism', required=True, help='Organism')
    parser.add_argument('--control', required=True, help='Category in the design table to serve as control in the differential expression')
    parser.add_argument('--samples', required=True, help='Samples to be compared against the control, comma separated, no spaces')

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
      "project": args.project,
      "options": {
        "organism": args.organism
      },
      "comparisons": {
        args.control:args.samples
	  },
      "tools_conf": {
        "GSVA": {
          "input": {
            "counts": args.counts,
            "design": args.design
            },
          "output": {
            "heatmap": args.heatmap
            },
          "tool_conf": {
            }
          }
        }
      }

    # Startup the logger format
    logger = pf.create_logger(config['log_files'][0])

    GSVA(config, 'GSVA')

if __name__ == "__main__":
    main()
