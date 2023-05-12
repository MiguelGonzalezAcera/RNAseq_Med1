import logging
import os
import mysql.connector
import pandas as pd
import python_scripts.python_functions as pf

def deseq2(config, tool_name):
    """Get the counts of a number of bam files in a directory
    """
    # Note in logger
    logging.info(f'Starting {tool_name} process')

    # Extract information
    # Inputs
    # counts file
    counts = config['tools_conf'][tool_name]['input']['counts']
    # design file
    design = config['tools_conf'][tool_name]['input']['design']

    # Outputs
    # marker file
    DEtouched = config['tools_conf'][tool_name]['output']['DEtouched']
    # outside directory
    out_dir = "/".join(DEtouched.split('/')[0:-1])
    # out object name
    out_obj = out_dir + "/" + config['project'] +".Rda"

    # Other information
    # comparison set
    samples = config['comparisons']
    # project name
    project = config['project']
    # organism used
    organism = config['options']['organism']

    # Create the command to run the deseq2 R script
    command = ""
    if not os.path.exists(out_dir):
        command += f"mkdir {out_dir};"

    # Create the deseq2 command for each control
    for control in samples:
        command += f'Rscript Rscripts/deseq2.r --counts {counts} --design {design} --out_obj {out_obj} --organism {organism} --control {control} --comparisons {samples[control]}; '
    # Touch the markerfile
    command += f'touch {DEtouched}'

    # Run the commans
    pf.run_command(command)

    # -----------------------------------------------------------------

    # Once all of the tables have been generated, we must load them into the mysql server
    # Establish the connection to the datbase and make a cursor
    mydb = mysql.connector.connect(
      host="localhost",
      user="root",
      passwd="Pl4ter!a",
      database="RNAseq",
      allow_local_infile=True
    )
    mycursor = mydb.cursor()

    # Create and execute the commands to generate and insert the table into the server per each table
    for control in samples:
        sample_ids = samples[control].split(",")
        for sample in sample_ids:
            # Get the name of the table to be uploaded to the database
            table_name = out_dir + "/" + config['project'] + "_" + f"{sample}_{control}_expanded.tsv"

            # Sometimes the pipeline is asked not to upload the things, so we must check beforehand
            if  config['options']['sql_load'] == "True":
                # Ask the database if table exists
                create_command = f"""show tables like '{project}_{sample}_{control}';"""

                # Run the command and fetch a result, which will be empty if no table is in the database by that name
                mycursor.execute(create_command)
                result = mycursor.fetchone()

                # If the table exists, just log it and be done with it
                if result:
                    logging.info(f'Table {project}_{sample}_{control} already exists')
                else:
                    # Get the header from the table to upload
                    with open(table_name) as f:
                        header = f.readline()
                    header = header.replace("\"", "").rstrip().split('\t')

                    # Get the header in a single string comma separated to upload later
                    j_header = ",".join(header)

                    # Make the command to upload the table
                    # Star by naming the table and defining the ID number field. Then iter through the heade defining each field
                    create_command = f"create table RNAseq.{project}_{sample}_{control}(id INT(11) NOT NULL AUTO_INCREMENT, "
                    for field in header:
                        if field == 'EnsGenes':
                            # Ensemble Gene ID is a character and must always exist
                            create_command += 'EnsGenes VARCHAR(255) NOT NULL, '
                        elif field == 'Genes':
                            # Gene name. character
                            create_command += 'Genes VARCHAR(255), '
                        elif field == 'FLAG':
                            # QC flag. character
                            create_command += 'FLAG VARCHAR(255), '
                        elif field == 'pvalue':
                            # p value. decimal
                            create_command += 'pvalue DECIMAL(30,30), '
                        elif field == 'padj':
                            # p adjusted value. decimal
                            create_command += 'padj DECIMAL(30,30), '
                        else:
                            # All other fields are the counts per sample. float
                            create_command += f'{field} FLOAT(10,5), '

                    # Close the command naming ID as the primary key and linking the ensemblID to our annotation
                    create_command += f'primary key(id), foreign key(EnsGenes) references Refs.{organism}_genes(EnsGenes));'

                    # log and run the command to create the table
                    logging.info(create_command)
                    mycursor.execute(create_command)

                    # create the command to insert the data in the newly created table
                    insert_command = f"""load data local infile '{table_name}' into table RNAseq.{project}_{sample}_{control} fields terminated by '\\t' enclosed by '"' lines terminated by '\\n' ignore 1 rows ({j_header});"""
                    
                    # log and insert the data
                    logging.info(insert_command)
                    mycursor.execute(insert_command)

    # Commit the changes and close the database
    mydb.commit()

    mycursor.close()
    mydb.close()