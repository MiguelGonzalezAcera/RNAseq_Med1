import os
import json
import argparse
import datetime
import python_scripts
import glob
import pandas as pd
import logging

# Get initial data
# Open the configuration json
with open(config['param'], 'r') as f:
    config_dict = json.load(f)

# Get the out folder, design path and project name
outfolder = config_dict['outfolder']
design = config_dict['design']
project = config_dict['project']

# Set and start logger
logging.basicConfig(filename=f'{outfolder}/RNAseq.log', level=logging.DEBUG, format='#[%(levelname)s]: - %(asctime)s - %(message)s')
logging.info(f'Starting RNAseq {project}')

# Get the design into a dataframe
design_file = pd.read_csv(design, sep='\t', index_col=0).reset_index()
design_file.columns = ['sample','tr']

# Read the fastq files
fastq_path = config_dict['fastq_path']

# Select file extension naming convention between single end and paired end
if config_dict['options']['reads'] == 'single':
    fastq_r1 = []

    for name in design_file['sample'].tolist():
        fastq_r1.append(glob.glob(f'{fastq_path}/{name}.fastq.gz')[0])
else:
    fastq_r1 = []

    for name in design_file['sample'].tolist():
        fastq_r1.append(glob.glob(f'{fastq_path}/{name}_1.fastq.gz')[0])

# raise error when no files are found in the selected path
if not fastq_r1:
    logger.error(f'FASTQ files not found in {fastq_path}')
    raise ValueError(f'FASTQ files not found in {fastq_path}')

# get the annotation files
annot_path = config_dict['tools_conf']['annot']

# ------------------Snakemake pipeline------------------
# Rules
rule Mapping:
    input:
        fastq_r1 = fastq_r1
    output:
        mappingtouched = f"{outfolder}/bamfiles/mappingtouched.txt",
        bamfof = f"{outfolder}/bamfiles/bam.fof"
    run:
        tool_name = 'mapping'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input._allitems()},
            'output': {i[0]: i[1] for i in output._allitems()},
            'software': {},
            'tool_conf': {
                "threads": "2"
            }
        }
        python_scripts.mapping.mapping(config_dict, tool_name)

rule FastQC:
    input:
        fastq_r1 = fastq_r1,
        bamdir = rules.Mapping.output.mappingtouched
    output:
        fastqctouched = f"{outfolder}/fastqc/fastqctouched.txt"
    run:
        tool_name = 'fastqc'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input._allitems()},
            'output': {i[0]: i[1] for i in output._allitems()},
            'software': {},
            'tool_conf': {
                "threads": "5"
            }
        }
        python_scripts.fastqc.fastqc(config_dict, tool_name)

rule BamQC:
    input:
        bamfof = rules.Mapping.output.bamfof
    output:
        bamqctouched = f"{outfolder}/bamqc/bamqctouched.txt"
    run:
        tool_name = 'bamqc'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input._allitems()},
            'output': {i[0]: i[1] for i in output._allitems()},
            'software': {},
            'tool_conf': {
            }
        }
        python_scripts.bamqc.bamqc(config_dict, tool_name)

rule Splicing:
    input:
        bamfof = rules.Mapping.output.bamfof,
        annot = annot_path,
        design = design
    output:
        splicetouched = f"{outfolder}/splicing/splicetouched.txt"
    run:
        tool_name = 'splicing'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input._allitems()},
            'output': {i[0]: i[1] for i in output._allitems()},
            'software': {},
            'tool_conf': {
                "threads": "15"
            }
        }
        python_scripts.splicing.splicing(config_dict, tool_name)

rule Counts:
    input:
        bamdir = rules.Mapping.output.mappingtouched,
        annot = annot_path
    output:
        counts = f"{outfolder}/counts.tsv"
    run:
        tool_name = 'get_counts'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input._allitems()},
            'output': {i[0]: i[1] for i in output._allitems()},
            'software': {},
            'tool_conf': {}
        }
        python_scripts.get_counts.counts(config_dict, tool_name)

rule PCA:
    input:
        counts = rules.Counts.output.counts,
        design = design
    output:
        pcatouched = f"{outfolder}/pca/pcatouched.txt"
    run:
        tool_name = 'pca'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input._allitems()},
            'output': {i[0]: i[1] for i in output._allitems()},
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
        norm_counts = f"{outfolder}/detables/{project}_norm_counts.Rda"
    run:
        tool_name = 'differential_expression'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input._allitems()},
            'output': {i[0]: i[1] for i in output._allitems()},
            'software': {},
            'tool_conf': {}
        }
        python_scripts.differential_expression.deseq2(config_dict, tool_name)

rule GSVA:
    input:
        norm_counts = rules.deseq2.output.norm_counts,
        design = design
    output:
        heatmap = f"{outfolder}/GSVA/{project}_GSVA_heatmap.png"
    run:
        tool_name = 'GSVA'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input._allitems()},
            'output': {i[0]: i[1] for i in output._allitems()},
            'software': {},
            'tool_conf': {}
        }
        python_scripts.GSVA.GSVA(config_dict, tool_name)

rule clustering_heatmap:
    input:
        DEtouched = rules.deseq2.output.DEtouched,
        norm_counts = rules.deseq2.output.norm_counts,
        design = design
    output:
        heatmap = f"{outfolder}/plots/{project}_clustering_heatmap.png"
    run:
        tool_name = 'clustering_heatmap'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input._allitems()},
            'output': {i[0]: i[1] for i in output._allitems()},
            'software': {},
            'tool_conf': {
                "dimensions": "2000,2000"
            }
        }
        python_scripts.clustering_heatmap.clustering_heatmap(config_dict, tool_name)

rule clustering_markers:
    input:
        norm_counts = rules.deseq2.output.norm_counts,
        design = design
    output:
        markerstouched = f"{outfolder}/markers/markerstouched.txt"
    run:
        tool_name = 'clustering_markers'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input._allitems()},
            'output': {i[0]: i[1] for i in output._allitems()},
            'software': {},
            'tool_conf': {
                "dimensions": "2000,2000"
            }
        }
        python_scripts.clustering_markers.clustering_markers(config_dict, tool_name)

rule volcano_plot:
    input:
        DEtouched = rules.deseq2.output.DEtouched
    output:
        volcanotouched = f"{outfolder}/plots/volcanotouched.txt"
    run:
        tool_name = 'volcano_plot'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input._allitems()},
            'output': {i[0]: i[1] for i in output._allitems()},
            'software': {},
            'tool_conf': {
                "dimensions": "10,6"
            }
        }
        python_scripts.volcano_plot.volcano_plot(config_dict, tool_name)

rule volcano_markers:
    input:
        DEtouched = rules.deseq2.output.DEtouched
    output:
        MVtouched = f"{outfolder}/markers/MVtouched.txt"
    run:
        tool_name = 'volcano_markers'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input._allitems()},
            'output': {i[0]: i[1] for i in output._allitems()},
            'software': {},
            'tool_conf': {
                "dimensions": "7,7"
            }
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
            'input': {i[0]: i[1] for i in input._allitems()},
            'output': {i[0]: i[1] for i in output._allitems()},
            'software': {},
            'tool_conf': {
                "dimensions": "2000,2000"
            }
        }
        python_scripts.GSEA_markers.GSEA_markers(config_dict, tool_name)

# NOTE: this one does not load the tables with the results or the counts, that's done in the DE step.
# This one just loads the design in a different database.
# The DEtouched thing is only there in case there is an error in the naming and the DE fails.
# It would run the DE and this at the same time and I don't want it to load something mistaken if the input's wrong
rule load_project:
    input:
        design = design,
        DEtouched = rules.deseq2.output.DEtouched
    output:
        prloadtouched = f"{outfolder}/detables/loadedtouched.txt"
    run:
        tool_name = 'load_project_table'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input._allitems()},
            'output': {i[0]: i[1] for i in output._allitems()},
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
            'input': {i[0]: i[1] for i in input._allitems()},
            'output': {i[0]: i[1] for i in output._allitems()},
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
            'input': {i[0]: i[1] for i in input._allitems()},
            'output': {i[0]: i[1] for i in output._allitems()},
            'software': {},
            'tool_conf': {}
        }
        python_scripts.GO.GO_enrichment(config_dict, tool_name)

rule report:
    input:
        markerstouched = rules.clustering_markers.output.markerstouched,
        MVtouched = rules.volcano_markers.output.MVtouched,
        GSEAMtouched = rules.GSEA_markers.output.GSEAMtouched,
        GSVAhmap = rules.GSVA.output.heatmap
    output:
        report = f"{outfolder}/report.pdf"
    run:
        tool_name = 'report'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input._allitems()},
            'output': {i[0]: i[1] for i in output._allitems()},
            'software': {},
            'tool_conf': {}
        }
        python_scripts.markers_report.report(config_dict, tool_name)

rule all:
    input:
        pca = rules.PCA.output.pcatouched,
        keggtouched = rules.KEGG.output.keggtouched,
        gotouched = rules.GO.output.gotouched,
        fastqctouched = rules.FastQC.output.fastqctouched,
        bamtouched = rules.BamQC.output.bamqctouched,
        volcanotouched = rules.volcano_plot.output.volcanotouched,
        heatmap = rules.clustering_heatmap.output.heatmap,
        prloadtouched = rules.load_project.output.prloadtouched,
        splicetouched = rules.Splicing.output.splicetouched,
        report = rules.report.output.report
    run:
        tool_name = 'all'
        config_dict['tools_conf'][tool_name] = {
            'input': {i[0]: i[1] for i in input._allitems()},
            'output': {i[0]: i[1] for i in output._allitems()},
            'software': {},
            'tool_conf': {}
        }

        # Construct a dictionary with the main results
        config_dict['results'] = {"results": [
            {
                "name": "Counts",
                "value": rules.Counts.output.counts
            },
            {
                "name": "Splicing",
                "value": "/".join(rules.Splicing.output.splicetouched.split('/')[0:-1])
            },
            {
                "name": "PCA",
                "value": "/".join(rules.PCA.output.pcatouched.split('/')[0:-1])
            },
            {
                "name": "Differential expression",
                "value": "/".join(rules.deseq2.output.DEtouched.split('/')[0:-1])
            },
            {
                "name": "Plots",
                "value": "/".join(rules.volcano_plot.output.volcanotouched.split('/')[0:-1])
            },
            {
                "name": "GO",
                "value": "/".join(rules.GO.output.gotouched.split('/')[0:-1])
            },
            {
                "name": "KEGG",
                "value": "/".join(rules.KEGG.output.keggtouched.split('/')[0:-1])
            },
            {
                "name": "Markers report",
                "value": rules.report.output.report
            }
        ]}

        # Save the dictionary
        # NOTE: This dictionary only stores the information added in the 'all' rule.
        # Other rules do not add or remove anything from the original object
        with open(config['param'], 'w') as f:
            json.dump(config_dict, f)

        logging.info(f'Finished RNAseq {project}')
