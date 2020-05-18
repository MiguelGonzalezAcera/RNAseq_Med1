import argparse
import logging
import os
import glob
import pandas as pd
import python_scripts.python_functions as pf


def volcano_plot(config, tool_name):
    """Get the
    """

    logging.info(f'Starting {tool_name} process')

    out_dir_DE = "/".join(config['tools_conf'][tool_name]['input']['DEtouched'].split('/')[0:-1])

    out_dir = "/".join(config['tools_conf'][tool_name]['output']['volcanotouched'].split('/')[0:-1])
    volcanotouched = config['tools_conf'][tool_name]['output']['volcanotouched']
    samples = config['comparisons']
    organism = config['options']['organism']

    # Create the command to run the pca R script
    command = ""
    if not os.path.exists(out_dir):
        command += f"mkdir {out_dir};"

    path_list = []

    for control in samples:
        sample_ids = samples[control].split(",")
        for sample in sample_ids:
            id_tab = out_dir + "/" + f"{sample}_{control}_volcano.png"
            id_sample = out_dir_DE + "/" + config['project'] + "_" + f"{sample}_{control}.Rda"
            id_obj = f"{sample}_{control}"

            command += f'Rscript /DATA/RNAseq_test/Scripts/Rscripts/volcano_plot.r --out_plot {id_tab} --res {id_sample} --organism {organism}; '

            path_list.append([id_sample, id_tab])

    command += f'touch {volcanotouched}'

    # Add volcano plots paths to the design table
    path_df = pd.DataFrame(path_list, columns = ['Robj_path', 'Volcano_path'])
    design_df = pd.read_csv(config['tools_conf'][tool_name]['input']['design_tab'], sep='\t', index_col=None)
    design_df = pd.merge(design_df, path_df, on=['Robj_path'])
    design_df.to_csv(config['tools_conf'][tool_name]['output']['design_tab'], sep='\t', index=False)

    pf.run_command(command)


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
