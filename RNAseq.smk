import os
import json
import argparse
import datetime
import python_scripts
import qc

# Get initial data
config_dict_path = "/VAULT/20200204_Timo_KO_mice/config.json"
with open(config_dict_path, 'r') as f:
    config_dict = json.load(f)

# Get the out folder
outfolder = config_dict['outfolder']

# Fastq files
fastq_r1 = config_dict['r1_files'].split(',')
fastq_r2 = config_dict['r2_files'].split(',')
fastq = fastq_r1 + fastq_r2
if config_dict['options']['reads'] == 'single':
    fastq.remove("")
design = config_dict['design']
project = config_dict['project']

bedfile_path = config_dict['tools_conf']['bedfile']
list_path = bedfile_path.replace('.bed','.list')
annot_path = config_dict['tools_conf']['annot']

# Rules
if config_dict['options']['reads'] == 'single':
    rule Mapping:
        input:
            fastq_r1 = config_dict['r1_files'].split(',')
        output:
            mappingtouched = f"{outfolder}/bamfiles/mappingtouched.txt"
        run:
            tool_name = 'mapping'
            config_dict['tools_conf'][tool_name] = {
                'input': {i[0]: i[1] for i in input.allitems()},
                'output': {i[0]: i[1] for i in output.allitems()},
                'software': {},
                'tool_conf': {
                    "threads": "2"
                }
            }
            python_scripts.mapping.mapping(config_dict, tool_name)
else:
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
                    "threads": "2"
                }
            }
            python_scripts.mapping.mapping(config_dict, tool_name)

rule GenerateRegions:
    input:
        bedfile = bedfile_path
    output:
        list = list_path
    run:
        tool_name = 'generate_regions'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        python_scripts.generate_regions.regions(config_dict, tool_name)

rule Counts:
    input:
        bamdir = rules.Mapping.output.mappingtouched,
        annot = annot_path
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
        bed = bedfile_path
    output:
        cpbtouched = f"{outfolder}/cpbfiles/cpbtouched.txt"
    run:
        tool_name = 'coverage_per_base'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        python_scripts.get_cpb.coverage_per_base(config_dict, tool_name)

rule Fastqc:
    input:
        fastq = fastq
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
        bedfile = bedfile_path
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
            'tool_conf': {}
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

rule deseq2:
    input:
        counts = rules.Counts.output.counts,
        design = design
    output:
        DEtouched = f"{outfolder}/detables/DEtouched.txt",
        design_tab = f"{outfolder}/detables/{project}_design.txt"
    run:
        tool_name = 'differential_expression'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        python_scripts.differential_expression.deseq2(config_dict, tool_name)

rule load_project:
    input:
        design_tab = rules.deseq2.output.design_tab
    output:
        prloadtouched = f"{outfolder}/detables/loadedtouched.txt"
    run:
        tool_name = 'load_project_table'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        python_scripts.load_design.load_design(config_dict, tool_name)

rule KEGG:
    input:
        DEtouched = rules.deseq2.output.DEtouched,
    output:
        keggtouched = f"{outfolder}/KEEG_enrichment/keggtouched.txt"
    run:
        tool_name = 'KEEG_enrichment'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        python_scripts.KEGG.KEGG_enrichment(config_dict, tool_name)

rule volcano_plot:
    input:
        DEtouched = rules.deseq2.output.DEtouched,
    output:
        volcanotouched = f"{outfolder}/plots/volcanotouched.txt"
    run:
        tool_name = 'volcano_plot   '
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        python_scripts.volcano_plot.volcano_plot(config_dict, tool_name)

rule GO:
    input:
        DEtouched = rules.deseq2.output.DEtouched,
    output:
        gotouched = f"{outfolder}/GO_enrichment/gotouched.txt"
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
        fastqceval = rules.Fastqc.output.fastqceval,
        bamdir = rules.Mapping.output.mappingtouched,
        counts = rules.Counts.output.counts,
        bamqc = rules.Bamqc.output.bamqctouched,
        pca = rules.PCA.output.pcatouched,
        bamqceval = rules.BamqcEval.output.evaloutfile,
        keggtouched = rules.KEGG.output.keggtouched,
        gotouched = rules.GO.output.gotouched,
        volcanotouched = rules.volcano_plot.output.volcanotouched,
        prloadtouched = rules.load_project.output.prloadtouched
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
