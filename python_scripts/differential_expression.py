import argparse
import logging
import os
import glob
import mysql.connector
import pandas as pd
import python_scripts.python_functions as pf

def deseq2(config, tool_name):
    """Get the counts of a number of bam files in a directory
    """

    counts = config['tools_conf'][tool_name]['input']['counts']
    design = config['tools_conf'][tool_name]['input']['design']
    samples = config['comparisons']
    project = config['project']
    out_dir = "/".join(config['tools_conf'][tool_name]['output']['DEtouched'].split('/')[0:-1])
    out_obj = out_dir + "/" + config['project'] +".Rda"
    organism = config['options']['organism']
    DEtouched = config['tools_conf'][tool_name]['output']['DEtouched']
    design_name = config['tools_conf'][tool_name]['output']['design_tab']

    # Create the command to run the deseq2 R script
    command = ""
    if not os.path.exists(out_dir):
        command += f"mkdir {out_dir};"

    for control in samples:
        command += f'Rscript /DATA/RNAseq_test/Scripts/Rscripts/deseq2.r --counts {counts} --design {design} --out_obj {out_obj} --organism {organism} --control {control} --comparisons {samples[control]}; '
    command += f'touch {DEtouched}'

    print(command)

    pf.run_command(command)

    # Once all of the tables have been generated, we must load them into the mysql server
    # Establish the connection to the datbase
    mydb = mysql.connector.connect(
      host="localhost",
      user="root",
      passwd="Plater1a",
      database="RNAseq"
    )

    mycursor = mydb.cursor()

    # Create the table for the design
    design_list = []

    # Create and execute the commands to generate and insert the table into the server per each table
    for control in samples:
        sample_ids = samples[control].split(",")
        for sample in sample_ids:
            table_name = out_dir + "/" + config['project'] + "_" + f"{sample}_{control}.tsv"

            if  config['options']['sql_load'] == "True":
                create_command = f"""create table RNAseq.{project}_{sample}_{control}(id INT(11) NOT NULL AUTO_INCREMENT, baseMean FLOAT(10,5) NOT NULL, log2FoldChange FLOAT(10,5) NOT NULL, lcfSE FLOAT(10,5) NOT NULL, stat FLOAT(10,5) NOT NULL, pvalue FLOAT(10,5) NOT NULL, padj FLOAT(10,5) NOT NULL, EnsGenes VARCHAR(255) NOT NULL, Genename VARCHAR(255), primary key(id), foreign key(EnsGenes) references Refs.mouse_genes(EnsGenes));"""
                print(create_command)
                mycursor.execute(create_command)
                insert_command = f"""load data local infile '{table_name}' into table RNAseq.{project}_{sample}_{control} fields terminated by '\\t' enclosed by '"' lines terminated by '\\n' ignore 1 rows (baseMean,log2FoldChange,lcfSE,stat,pvalue,padj,EnsGenes,Genename);"""
                print(insert_command)
                mycursor.execute(insert_command)

            # Add the line to the list:
            line = [f"{project}_{sample}_{control}",f"{control}",f"{sample}",table_name,table_name.replace(".tsv",".Rda")]
            design_list.append(line)

    # Transform the project into table and save it
    df = pd.DataFrame(design_list, columns=['Comparison','Control','Sample','Table_path','Robj_path'])
    df.to_csv(design_name, sep='\t', index=False)

    mydb.commit()

    mycursor.close()
    mydb.close()


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