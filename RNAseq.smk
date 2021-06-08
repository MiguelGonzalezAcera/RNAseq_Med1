import os
import json
import argparse
import datetime
import python_scripts
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
design = config_dict['design']
project = config_dict['project']

# Set logger
logging.basicConfig(filename=f'{outfolder}/RNAseq.log', level=logging.DEBUG, format='#[%(levelname)s]: - %(asctime)s - %(message)s')
logging.info(f'Starting RNAseq {project}')

# Fastq files
design_file = pd.read_csv(design, sep='\t', index_col=0).reset_index()
design_file.columns = ['sample','tr']

fastq_path = config_dict['fastq_path']

if config_dict['options']['reads'] == 'single':
    fastq_r1 = []

    for name in design_file['sample'].tolist():
        fastq_r1.append(glob.glob(f'{fastq_path}/{name}.fastq.gz')[0])
else:
    fastq_r1 = []
    fastq_r2 = []

    for name in design_file['sample'].tolist():
        fastq_r1.append(glob.glob(f'{fastq_path}/{name}_1.fastq.gz')[0])
        fastq_r2.append(glob.glob(f'{fastq_path}/{name}_2.fastq.gz')[0])

if not fastq_r1:
    logger.error(f'FASTQ files not found in {fastq_path}')
    raise ValueError(f'FASTQ files not found in {fastq_path}')

bedfile_path = config_dict['tools_conf']['bedfile']
list_path = bedfile_path.replace('.bed','.list')
annot_path = config_dict['tools_conf']['annot']


# Rules
if config_dict['options']['reads'] == 'single':
    rule Mapping:
        input:
            fastq_r1 = fastq_r1
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
            fastq_r1 = fastq_r1,
            fastq_r2 = fastq_r2
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

if config_dict['options']['splicing'] == 'True':
    include: './subworkflows/splicing.smk'
    splicetouched = rules.Splicing.output.splicetouched
else:
    splicetouched = config['param']

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

if config_dict['options']['qc'] == 'True':
    include: './subworkflows/qc.smk'
    cpbtouched = rules.CoveragePerBase.output.cpbtouched
    fastqceval = rules.Fastqc.output.fastqceval
    bamqceval = rules.BamqcEval.output.evaloutfile
else:
    cpbtouched = config['param']
    fastqceval = config['param']
    bamqceval = config['param']

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
        design_tab = f"{outfolder}/detables/{project}_design.tmp",
        norm_counts = f"{outfolder}/detables/{project}_norm_counts.Rda"
    run:
        tool_name = 'differential_expression'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        python_scripts.differential_expression.deseq2(config_dict, tool_name)

rule clustering_heatmap:
    input:
        DEtouched = rules.deseq2.output.DEtouched,
        norm_counts = rules.deseq2.output.norm_counts
    output:
        heatmap = f"{outfolder}/plots/{project}_clustering_heatmap.png"
    run:
        tool_name = 'clustering_heatmap'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        python_scripts.clustering_heatmap.clustering_heatmap(config_dict, tool_name)

rule clustering_markers:
    input:
        norm_counts = rules.deseq2.output.norm_counts
    output:
        markerstouched = f"{outfolder}/markers/markerstouched.txt"
    run:
        tool_name = 'clustering_markers'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        python_scripts.clustering_markers.clustering_markers(config_dict, tool_name)

rule volcano_markers:
    input:
        DEtouched = rules.deseq2.output.DEtouched
    output:
        MVtouched = f"{outfolder}/markers/MVtouched.txt"
    run:
        tool_name = 'volcano_markers'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        python_scripts.volcano_markers.volcano_markers(config_dict, tool_name)

rule GSEA_markers:
    input:
        DEtouched = rules.deseq2.output.DEtouched
    output:
        GSEAMtouched = f"{outfolder}/markers/GSEAMtouched.txt"
    run:
        tool_name = 'GSEA_markers'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        python_scripts.GSEA_markers.GSEA_markers(config_dict, tool_name)

rule volcano_plot:
    input:
        DEtouched = rules.deseq2.output.DEtouched,
        design_tab = rules.deseq2.output.design_tab
    output:
        volcanotouched = f"{outfolder}/plots/volcanotouched.txt",
        design_tab = f"{outfolder}/detables/{project}_design.txt"
    run:
        tool_name = 'volcano_plot'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        python_scripts.volcano_plot.volcano_plot(config_dict, tool_name)

rule load_project:
    input:
        design_tab = rules.volcano_plot.output.design_tab
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
        keggtouched = f"{outfolder}/KEGG_enrichment/keggtouched.txt"
    run:
        tool_name = 'KEGG_enrichment'
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

rule report:
    input:
        markerstouched = rules.clustering_markers.output.markerstouched,
        MVtouched = rules.volcano_markers.output.MVtouched,
        GSEAMtouched = rules.GSEA_markers.output.GSEAMtouched
    output:
        report = f"{outfolder}/report.pdf"
    run:
        tool_name = 'report'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input.allitems()},
            'output': {i[0]: i[1] for i in output.allitems()},
            'software': {},
            'tool_conf': {}
        }
        python_scripts.markers_report.report(config_dict, tool_name)

rule all:
    input:
        fastqceval = fastqceval,
        cpbtouched = cpbtouched,
        splicetouched = splicetouched,
        pca = rules.PCA.output.pcatouched,
        bamqceval = bamqceval,
        keggtouched = rules.KEGG.output.keggtouched,
        gotouched = rules.GO.output.gotouched,
        volcanotouched = rules.volcano_plot.output.volcanotouched,
        heatmap = rules.clustering_heatmap.output.heatmap,
        prloadtouched = rules.load_project.output.prloadtouched,
        report = rules.report.output.report
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
            "name": "Counts",
            "value": rules.Counts.output.counts
        }]
        }

        print(config_dict['results'])

        with open(config['param'], 'w') as f:
            json.dump(config_dict, f)

        logging.info(f'Finished RNAseq {project}')
