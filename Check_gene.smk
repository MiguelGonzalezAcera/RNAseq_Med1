import os
import json
import argparse
import datetime
import Gene_query
import qc
import glob
import pandas as pd
import logging

# Get initial data
# config = {}
# config['param'] = '/VAULT/20200325_Cytrobacter_GSE71734/config.json'

with open(config['param'], 'r') as f:
    config_dict = json.load(f)

# Get the out folder
outfolder = config_dict['outfolder']
jobid = config_dict['jobid']
genename = config_dict['genename']

if not os.path.exists(f"{outfolder}/static/{jobid}"):
    command = f"mkdir {outfolder}/static/{jobid}; mkdir {outfolder}/templates/{jobid}"
    output = os.system(command)

# Set logger
logging.basicConfig(filename=f'{outfolder}/Check_gene.log', level=logging.DEBUG, format='#[%(levelname)s]: - %(asctime)s - %(message)s')
logging.info(f'Starting consult for {genename}')

# Rules
rule gene_consult:
    output:
        FC_models_table = f"{outfolder}/static/{jobid}/{genename}_FC_models_table.tsv",
        FC_course_table = f"{outfolder}/static/{jobid}/{genename}_FC_course_table.tsv",
        barplot_counts = f"{outfolder}/static/{jobid}/{genename}_barplot_counts.png",
        barplot_FC = f"{outfolder}/static/{jobid}/{genename}_barplot_FC.png",
        course_plot = f"{outfolder}/static/{jobid}/{genename}_course_plot.png"
    run:
        tool_name = 'gene_consult'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        Gene_query.gene_consult.gene_consult(config_dict, tool_name)

rule report:
    input:
        FC_models_table = rules.gene_consult.output.FC_models_table,
        FC_course_table = rules.gene_consult.output.FC_course_table,
        barplot_counts = rules.gene_consult.output.barplot_counts,
        barplot_FC = rules.gene_consult.output.barplot_FC,
        course_plot = rules.gene_consult.output.course_plot,
    output:
        report = f"{outfolder}/static/{jobid}/{genename}_report.pdf"
    run:
        tool_name = 'report'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        Gene_query.report.report(config_dict, tool_name)

rule html:
    input:
        FC_models_table = rules.gene_consult.output.FC_models_table,
        FC_course_table = rules.gene_consult.output.FC_course_table,
        barplot_counts = rules.gene_consult.output.barplot_counts,
        barplot_FC = rules.gene_consult.output.barplot_FC,
        course_plot = rules.gene_consult.output.course_plot,
        report = rules.report.output.report
    output:
        html = f"{outfolder}/templates/{jobid}/{genename}_report.html"
    run:
        tool_name = 'html'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        Gene_query.single_gene_html.single_html(config_dict, tool_name)

rule all:
    input:
        report = rules.report.output.report,
        html = rules.html.output.html
    run:
        tool_name = 'all'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }

        config_dict['results'] = {"results": [
        {
            "name": "REPORT",
            "value": rules.report.output.report
        }]
        }

        print(config_dict['results'])

        with open(config['param'], 'w') as f:
            json.dump(config_dict, f)

        logging.info(f'Finished gene query for {genename}')
