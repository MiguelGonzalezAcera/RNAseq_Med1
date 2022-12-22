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
    logging.info(f'Starting {tool_name} process')

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
        command += f'Rscript Rscripts/deseq2.r --counts {counts} --design {design} --out_obj {out_obj} --organism {organism} --control {control} --comparisons {samples[control]}; '
    command += f'touch {DEtouched}'

    pf.run_command(command)

    # Once all of the tables have been generated, we must load them into the mysql server
    # Establish the connection to the datbase
    mydb = mysql.connector.connect(
      host="localhost",
      user="root",
      passwd="Plater1a",
      database="RNAseq",
      allow_local_infile=True
    )

    mycursor = mydb.cursor()

    # Create the table for the design
    design_list = []

    # Create and execute the commands to generate and insert the table into the server per each table
    for control in samples:
        sample_ids = samples[control].split(",")
        for sample in sample_ids:
            table_name = out_dir + "/" + config['project'] + "_" + f"{sample}_{control}_expanded.tsv"

            if  config['options']['sql_load'] == "True":
                # Check if table exists
                create_command = f"""show tables like '{project}_{sample}_{control}';"""

                mycursor.execute(create_command)

                result = mycursor.fetchone()

                if result:
                    logging.info(f'Table {project}_{sample}_{control} already exists')
                else:
                    with open(table_name) as f:
                        header = f.readline()
                    header = header.replace("\"", "").rstrip().split('\t')
                    j_header = ",".join(header)

                    create_command = f"create table RNAseq.{project}_{sample}_{control}(id INT(11) NOT NULL AUTO_INCREMENT, "
                    for field in header:
                        if field == 'EnsGenes':
                            create_command += 'EnsGenes VARCHAR(255) NOT NULL, '
                        elif field == 'Genes':
                            create_command += 'Genes VARCHAR(255), '
                        elif field == 'FLAG':
                            create_command += 'FLAG VARCHAR(255), '
                        elif field == 'pvalue':
                            create_command += 'pvalue DECIMAL(30,30), '
                        elif field == 'padj':
                            create_command += 'padj DECIMAL(30,30), '
                        else:
                            create_command += f'{field} FLOAT(10,5), '

                    create_command += f'primary key(id), foreign key(EnsGenes) references Refs.{organism}_genes(EnsGenes));'

                    logging.info(create_command)
                    mycursor.execute(create_command)
                    insert_command = f"""load data local infile '{table_name}' into table RNAseq.{project}_{sample}_{control} fields terminated by '\\t' enclosed by '"' lines terminated by '\\n' ignore 1 rows ({j_header});"""
                    logging.info(insert_command)
                    mycursor.execute(insert_command)

            # Add the line to the list:
            line = [f"{project}_{sample}_{control}",f"{control}",f"{sample}",table_name,table_name.replace("_expanded.tsv",".Rda")]
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
    parser.add_argument('--project', required=True, help='Project name')
    parser.add_argument('--out_dir', required=True, help='Directory for all of the tables)')
    parser.add_argument('--organism', required=True, help='Organism')

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

    out_dir = "/".join(args.out_dir.split('/')[0:-1])
    project = args.project

    config = {
      "DEBUG": args.debug,
      "TESTING": args.test,
      "DRY_RUN": args.dry_run,
      "log_files": ["/tmp/full.log"],
      "project": project,
      "options": {
        "organism": args.organism,
        "sql_load": "True"
      },
      "comparisons": {
		"fl_mock":"fl_AA,ko_mock",
        "fl_AA":"ko_AA",
        "ko_mock":"ko_AA"
	  },
      "tools_conf": {
        "differential_expression": {
          "input": {
            "counts": args.counts,
            "design": args.design
            },
          "output": {
            "DEtouched": args.out_dir,
            "design_tab": f"{out_dir}/{project}_design.tmp",
            "norm_counts": f"{out_dir}/{project}_norm_counts.Rda"
            },
          "tool_conf": {
            }
          }
        }
      }

    # Startup the logger format
    logger = pf.create_logger(config['log_files'][0])

    deseq2(config, 'differential_expression')


if __name__ == "__main__":
    main()
