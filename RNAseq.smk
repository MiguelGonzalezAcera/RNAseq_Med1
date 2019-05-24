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
        mappingtouched = f"{outfolder}/bamfiles/mappingtouched.txt"
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

rule GenerateRegions:
    input:
        bedfile = "/DATA/references/star_genomes/mmu38/annotation/Mus_musculus.GRCm38.96.merged.sorted.bed"
    output:
        list = "/DATA/references/star_genomes/mmu38/annotation/Mus_musculus.GRCm38.96.merged.sorted.list"
    run:
        tool_name = 'generate_regions'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {
                "genome": "/DATA/references/star_genomes/mmu38/sequence/Mus_musculus.GRCm38.dna.toplevel.dict"
            }
        }
        python_scripts.generate_regions.regions(config_dict, tool_name)

rule Counts:
    input:
        bamdir = rules.Mapping.output.mappingtouched,
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

rule CoveragePerBase:
    input:
        bamdir = rules.Mapping.output.mappingtouched,
        bed = "/DATA/references/star_genomes/mmu38/annotation/Mus_musculus.GRCm38.96.merged.sorted.bed"
    output:
        cpbtouched = f"{outfolder}/cpbfiles/cpbtouched.txt"
    run:
        tool_name = 'coverage_per_base'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {
                "genome": "/DATA/references/star_genomes/mmu38/sequence/Mus_musculus.GRCm38.dna.toplevel.txt"
            }
        }
        python_scripts.get_cpb.coverage_per_base(config_dict, tool_name)

rule Fastqc:
    input:
        fastq_r1 = config_dict['r1_files'].split(','),
        fastq_r2 = config_dict['r2_files'].split(',')
    output:
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

rule Bamqc:
    input:
        bamdir = rules.Mapping.output.mappingtouched,
        cpbtouched = rules.CoveragePerBase.output.cpbtouched,
        list = rules.GenerateRegions.output.list,
        bedfile = "/DATA/references/star_genomes/mmu38/annotation/Mus_musculus.GRCm38.96.merged.sorted.bed"
    output:
        bamqctouched = f"{outfolder}/bamqc/bamqctouched.txt"
    run:
        tool_name = 'bamqc_eval'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        qc.bamqc.bamqc(config_dict, tool_name)

rule BamqcEval:
    input:
        bamqctouched = rules.Bamqc.output.bamqctouched,
    output:
        evaloutfile = f"{outfolder}/bamqc/bamqceval.json"
    run:
        tool_name = 'bamqc'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {
                "genome": "/DATA/references/star_genomes/mmu38/sequence/Mus_musculus.GRCm38.dna.toplevel.fa"
            }
        }
        qc.bamqc_eval.bamqcEval(config_dict, tool_name)

rule PCA:
    input:
        counts = rules.Counts.output.counts,
        design = design
    output:
        pcatouched = f"{outfolder}/pca/pcatouched.txt"
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
        bamdir = rules.Mapping.output.mappingtouched,
        counts = rules.Counts.output.counts,
        bamqc = rules.Bamqc.output.bamqctouched,
        pca = rules.PCA.output.pcatouched,
        bamqceval = rules.BamqcEval.output.evaloutfile
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
            "value": rules.Mapping.output.mappingtouched
        },
        {
            "name": "Counts",
            "value": rules.Counts.output.counts
        },
        {
            "name": "Fastqc",
            "value": rules.Fastqc.output.fastqceval
        },
        {
            "name": "Bamqc",
            "value": rules.Bamqc.output.bamqctouched
        },
        {
            "name": "BamqcEval",
            "value": rules.BamqcEval.output.evaloutfile
        }]
        }

        print(config_dict['results'])
