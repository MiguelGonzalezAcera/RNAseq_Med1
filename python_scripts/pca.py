import argparse
import logging
import os
import imageio
import python_scripts.python_functions as pf

def gif(filelist, out_dir):
    """Create a gif file from a list of images
    """
    # result name
    result = f"{out_dir}/pca_3d.gif"

    # Make the gif using the list of files provided
    with imageio.get_writer(result, mode='I', fps=25) as writer:
        for filename in filelist:
            image = imageio.imread(filename)
            writer.append_data(image)

def pca(config, tool_name):
    """Get the counts of a number of bam files in a directory
    """
    logging.info(f'Starting {tool_name} process')

    # Get parameters
    # Inputs
    # Counts files
    tr_counts = config['tools_conf'][tool_name]['input']['tr_counts'].replace(".Rda", ".tsv")
    # Design file
    design = config['tools_conf'][tool_name]['input']['design']

    # Outputs
    # Control files
    pcatouched = config['tools_conf'][tool_name]['output']['pcatouched']
    # Out directory
    out_dir = "/".join(pcatouched.split('/')[0:-1])

    # Create the command to run the pca R script
    command = ""
    
    # Make directory if it soes not exist
    if not os.path.exists(out_dir):
        command += f"mkdir {out_dir};"

    # Command
    command += f'Rscript Rscripts/pca.r --counts {tr_counts} --design {design} --out_dir {out_dir}; '
    command += f'touch {pcatouched}'

    # Run command
    pf.run_command(command)

    # List files to make gif if needed
    try:
        filelist = pf.list_files_dir(out_dir, ext = "*3d*")
    except:
        filelist = []

    # Run the creation of a gif if filelist is not empty
    if len(filelist) != 0:
        gif(sorted(filelist), out_dir)

def get_arguments():
    """
    Function that parse arguments given by the user, returning a dictionary
    that contains all the values.
    """

    # Create the top-level parser
    parser = argparse.ArgumentParser()

    # Mandatory variables
    parser.add_argument('--counts', required=True, help='Table with the counts of the assay, straight from featurecounts (so far)')
    parser.add_argument('--design', required=True, help='Table with the design of the experiment')
    parser.add_argument('--out_dir', required=True, help='Directory for all of the plots)')

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
        "pca": {
          "input": {
            "counts": args.counts,
            "design": args.design
            },
          "output": {
            "out_dir": args.out_dir
            },
          "tool_conf": {
            }
          }
        }
      }

    # Startup the logger format
    logger = pf.create_logger(config['log_files'][0])

    pca(config, 'pca')


if __name__ == "__main__":
    main()
