import argparse
import logging
import os
import glob
import mysql.connector
import pandas as pd
import python_scripts.python_functions as pf

def update_design(config, tool_name):
    """Get the counts of a number of bam files in a directory
    """

    design_tab = config['tools_conf'][tool_name]['input']['design_tab']
    project = config['project']
    organism = config['options']['organism']
    prupdtouched = config['tools_conf'][tool_name]['output']['prupdtouched']

    # Establish the connection to the datbase
    mydb = mysql.connector.connect(
      host="localhost",
      user="root",
      passwd="Plater1a",
      database="Projects"
    )

    mycursor = mydb.cursor()

    df = pd.read_csv(design_tab, sep='\t', index_col=None)

    for index, row in df.iterrows():
        insert_command = f"""insert into Projects.{project}(Comparison,Control,Sample,Table_path,Robj_path,Volcano_path) value('{row[0]}','{row[1]}','{row[2]}','{row[3]}','{row[4]}','{row[5]}');"""
        logging.info(insert_command)
        mycursor.execute(insert_command)

    mydb.commit()

    mycursor.close()
    mydb.close()

    command = f'touch {prupdtouched}'
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

    update_design(config, 'pca')


if __name__ == "__main__":
    main()
