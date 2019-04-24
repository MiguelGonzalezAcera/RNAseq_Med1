import os
import json
import argparse
import datetime
import python_functions as pf

logger = pf.create_logger(config_dict['log_files'][0])

rule Mapping
    input:
        fastq_r1 = config_dict['r1_files'].split(','),
        fastq_r2 = config_dict['r2_files'].split(','),
        samplelist = config_dict['samples'].split(',')
    output:
        bamdir = f"{outfolder}/bamfiles"
    run:
        tool_name = 'mapping'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        python_scripts.mapping.mapping(config_dict, tool_name)

rule cpb
    input:
        bamdir = rule.Mapping.output.bamdir
    output:
        cpbdir = f"{outfolder}/cpbfiles"
    run:
        tool_name = 'get_cpb'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        python_scripts.get_cpb.cpb(config_dict, tool_name)

rule Counts
    input:
        bamdir = rule.Mapping.output.bamdir
    output:
        counts = f"{outfolder}/counts.tsv"
    run:
        tool_name = 'get_counts'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        python_scripts.get_counts.counts(config_dict, tool_name)
