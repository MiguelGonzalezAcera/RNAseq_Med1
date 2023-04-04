import logging
import mysql.connector
import python_scripts.python_functions as pf

def load_design(config, tool_name):
    """Get the counts of a number of bam files in a directory
    """
    # Record in the log
    logging.info(f'Starting {tool_name} process')

    # Get information
    # Design file
    design = config['tools_conf'][tool_name]['input']['design']
    # Marker file
    prloadtouched = config['tools_conf'][tool_name]['output']['prloadtouched']
    # Project name
    project = config['project']

    # Establish the connection to the DESIGN datbase
    mydb = mysql.connector.connect(
      host="localhost",
      user="root",
      passwd="Pl4ter!a",
      database="Designs",
      allow_local_infile=True
    )

    mycursor = mydb.cursor()

    # Check if table already exists
    create_command = f"""show tables like '{project}';"""

    mycursor.execute(create_command)
    result = mycursor.fetchone()

    # If the table does not exist, load the design into the database
    if result:
        logging.info(f'Table {project} already exists in Designs')
    else:
        create_command = f"""create table Designs.{project}(ID INT unique auto_increment, Sample VARCHAR(255) NOT NULL, Treatment VARCHAR(255) NOT NULL, primary key(ID));"""
        logging.info(create_command)
        mycursor.execute(create_command)

        insert_command = f"""load data local infile '{design}' into table Designs.{project} fields terminated by '\\t' enclosed by '"' lines terminated by '\\n' ignore 1 rows (Sample,Treatment);"""
        logging.info(insert_command)
        mycursor.execute(insert_command)

    mydb.commit()

    mycursor.close()
    mydb.close()

    # Create the marker file for snakemake
    command = f'touch {prloadtouched}'
    pf.run_command(command)
