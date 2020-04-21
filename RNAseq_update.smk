import os
import json
import argparse
import datetime
import python_scripts
import qc

# Get initial data
with open(config['param'], 'r') as f:
    config_dict = json.load(f)

# Get the out folder
outfolder = config_dict['outfolder']

design = config_dict['design']
project = config_dict['project']

reanalisis_numm = config_dict['reanalisis_numm']

# Rules
rule deseq2:
    input:
        counts = f"{outfolder}/counts.tsv",
        design = design
    output:
        DEtouched = f"{outfolder}/detables/DEtouched_r{reanalisis_numm}.txt",
        design_tab = f"{outfolder}/detables/{project}_design_r{reanalisis_numm}.tmp"
    run:
        tool_name = 'differential_expression'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        python_scripts.differential_expression.deseq2(config_dict, tool_name)

rule volcano_plot:
    input:
        DEtouched = rules.deseq2.output.DEtouched,
        design_tab = rules.deseq2.output.design_tab
    output:
        volcanotouched = f"{outfolder}/plots/volcanotouched_r{reanalisis_numm}.txt",
        design_tab = f"{outfolder}/detables/{project}_design_r{reanalisis_numm}.txt"
    run:
        tool_name = 'volcano_plot   '
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        python_scripts.volcano_plot.volcano_plot(config_dict, tool_name)

rule update_project:
    input:
        design_tab = rules.volcano_plot.output.design_tab
    output:
        prupdtouched = f"{outfolder}/detables/loadedtouched_r{reanalisis_numm}.txt"
    run:
        tool_name = 'update_project_table'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        python_scripts.update_design.update_design(config_dict, tool_name)

rule KEGG:
    input:
        DEtouched = rules.deseq2.output.DEtouched,
    output:
        keggtouched = f"{outfolder}/KEEG_enrichment/keggtouched_r{reanalisis_numm}.txt"
    run:
        tool_name = 'KEEG_enrichment'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        python_scripts.KEGG.KEGG_enrichment(config_dict, tool_name)

rule GO:
    input:
        DEtouched = rules.deseq2.output.DEtouched,
    output:
        gotouched = f"{outfolder}/GO_enrichment/gotouched_r{reanalisis_numm}.txt"
    run:
        tool_name = 'GO_enrichment'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        python_scripts.GO.GO_enrichment(config_dict, tool_name)

rule all:
    input:
        keggtouched = rules.KEGG.output.keggtouched,
        gotouched = rules.GO.output.gotouched,
        volcanotouched = rules.volcano_plot.output.volcanotouched,
        prloadtouched = rules.update_project.output.prupdtouched
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
        }]
        }

        print(config_dict['results'])
