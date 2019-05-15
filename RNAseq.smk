import os
import json
import argparse
import datetime
import python_scripts
import qc

# Get initial data
config_dict_path = "/DATA/RNAseq_test/Test_project/TEST001_config.json"
with open(config_dict_path, 'r') as f:
    config_dict = json.load(f)

# Get the out folder
outfolder = config_dict['outfolder']

# Fastq files
fastq_r1 = config_dict['r1_files'].split(',')
fastq_r2 = config_dict['r2_files'].split(',')
design = config_dict['design']

# Rules
rule Mapping:
    input:
        fastq_r1 = config_dict['r1_files'].split(','),
        fastq_r2 = config_dict['r2_files'].split(',')
    output:
        bam_dir = f"{outfolder}/bamfiles"
    run:
        tool_name = 'mapping'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {
                "genome": "/DATA/references/star_genomes/mmu38/star_indices_overhang150/",
                "threads": "2"
            }
        }
        python_scripts.mapping.mapping(config_dict, tool_name)

rule Counts:
    input:
        bamdir = rules.Mapping.output.bam_dir,
        annot = "/DATA/references/star_genomes/mmu38/annotation/Mus_musculus.GRCm38.96.gtf"
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

rule Fastqc:
    input:
        fastq_r1 = config_dict['r1_files'].split(','),
        fastq_r2 = config_dict['r2_files'].split(',')
    output:
        fastqcdir = f"{outfolder}/fastqc",
        fastqceval = f"{outfolder}/fastqc/fastqc.results.txt"
    run:
        tool_name = 'fastqc'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {
                'threads': "2"
            }
        }
        qc.fastqc.fastqc(config_dict, tool_name)

rule PCA:
    input:
        counts = rules.Counts.output.counts,
        design = design
    output:
        out_dir = f"{outfolder}/pca"
    run:
        tool_name = 'pca'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        python_scripts.pca.pca(config_dict, tool_name)

rule all:
    input:
        fastqceval = rules.Fastqc.output.fastqceval,
        bamdir = rules.Mapping.output.bam_dir,
        counts = rules.Counts.output.counts,
        pca = rules.PCA.output.out_dir
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
            "name": "BAMDIR",
            "value": rules.Mapping.output.bam_dir
        },
        {
            "name": "Counts",
            "value": rules.Counts.output.counts
        },
        {
            "name": "Fastqc",
            "value": rules.Fastqc.output.fastqceval
        }]
        }

        print(config_dict['results'])

# rule cpb:
#     input:
#         bamdir = rule.Mapping.output.bamdir
#     output:
#         cpbdir = f"{outfolder}/cpbfiles"
#     run:
#         tool_name = 'coverage_per_base'
#         config_dict['tools_conf'][tool_name] = {
#             'input': {i[0]: i[1] for i in input.allitems()},
#             'output': {i[0]: i[1] for i in output.allitems()},
#             'software': {},
#             'tool_conf': {}
#         }
#         python_scripts.get_cpb.cpb(config_dict, tool_name)
