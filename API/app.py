#!/usr/bin/env python3
from flask import Flask, jsonify, make_response, request
from snakemake import snakemake
import threading
from flask_cors import cross_origin
import json
import os
import sys
import subprocess
import datetime

app = Flask(__name__)

def get_script_names():
    data = {
        'RNAseq': {'path': '/DATA/RNAseq_test/Scripts/RNAseq.smk'}
    }

    return data

def get_references_names():
    data = {
        'human': {
            "tools_conf": {
                "genome": "/DATA/references/star_genomes/hs38/Gencode26/",
                "genomedict": "/DATA/references/star_genomes/hs38/Gencode26/GRCh38.primary_assembly.genome.dict",
                "bedfile":"/DATA/references/star_genomes/hs38/Gencode26/Homo_sapiens.GRCh38.96.corr.bed",
                "annot": "/DATA/references/star_genomes/hs38/Gencode26/gencode.v26.primary_assembly.annotation.gtf",
                "genometxt": "/DATA/references/star_genomes/hs38/Gencode26/GRCh38.primary_assembly.genome.txt",
                "genomefasta": "/DATA/references/star_genomes/hs38/Gencode26/GRCh38.primary_assembly.genome.fa"
                }
            },
        'mouse': {
            "tools_conf": {
                "genome": "/DATA/references/star_genomes/mmu38/star_indices_overhang150/",
                "genomedict": "/DATA/references/star_genomes/mmu38/sequence/Mus_musculus.GRCm38.dna.toplevel.dict",
                "bedfile":"/DATA/references/star_genomes/mmu38/annotation/Mus_musculus.GRCm38.96.merged.sorted.bed",
                "annot": "/DATA/references/star_genomes/mmu38/annotation/Mus_musculus.GRCm38.96.gtf",
                "annot_gff3": "/DATA/references/star_genomes/mmu38/annotation/Mus_musculus.GRCm38.96.gff3",
                "genometxt": "/DATA/references/star_genomes/mmu38/sequence/Mus_musculus.GRCm38.dna.toplevel.txt",
                "genomefasta": "/DATA/references/star_genomes/mmu38/sequence/Mus_musculus.GRCm38.dna.toplevel.fa"
                }
            }
    }

    return data

def generate_configuration(postdata, pipeline):
        # Get the configuration of the different pipelines
    refs_names = get_references_names()

    config_json = {**postdata, **refs_names[postdata['options']['organism']]}

    # Serialize class attributes into a configuration fileW
    config_json_path = config_json['outfolder'] + '/config.json'
    with open(config_json_path, 'w') as outfile:
        json.dump(config_json, outfile)

    # Config path must be return in order to send it to the pipeline
    return True, config_json_path


def generate_response(postdata, dag):

    response = jsonify({"date": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                        "dag": str(dag), "status": "ok"}), 200

    return make_response(response)


def generate_dag(snakemakepath, config_json_path, pipeline):

    # Hago una primera ejecucion en vacio para generar el DAG
    cmd = f"snakemake --snakefile {snakemakepath} " +\
          f'--config param={config_json_path} --dag'
    p = subprocess.Popen(cmd, shell=True, executable='/bin/bash',
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         universal_newlines=True)
    text = p.stdout.read()
    # retcode = p.wait()

    return text

def launch_process(config_json_path, postdata, pipeline):

    # Return name of the pipeline scripts
    config_names = get_script_names()

    print('#[INFO] - Launching process')

    # Create the class snakemake and launch the process
    # snk_status = snakemake(config_names[pipeline]['path'],
    #                        config={"param": config_json_path, "postdata": str(postdata)},
    #                        lock=False, targets=['all'], cores=1, force_incomplete=True,
    #                        use_conda=True, debug=True,
    #                        )
    #
    # # If job broke during execution, call endpoint to delete job
    # if snk_status is False:
    #     print('#[ERR] - Cancelling job...')

    cmd = f"snakemake -s {config_names[pipeline]['path']} all -j 10 " +\
        f"--config param={config_json_path} --use-conda"

    os.system(cmd)


def launch_job(postdata, config_json_path, pipeline, mode='RUN'):
    # Get path to the needed scripts
    config_names = get_script_names()

    if mode == 'TEST':
        dag = {}
        snakemake(config_names[pipeline]['path'],
                  config={"param": config_json_path}, lock=False,
                  targets=['all'], cores=8, use_conda=True)

    else:
        dag = generate_dag(config_names[pipeline]['path'], config_json_path, pipeline)

        # Launching as class would block the response until the pipeline has finished.
        # https://stackoverflow.com/questions/21284319/can-i-make-one-method-of-class-run-in-background
        thread = threading.Thread(target=launch_process, args=(config_json_path, postdata, pipeline,))
        thread.start()

    return True, dag


@app.errorhandler(404)
def not_found(error):
    return make_response(jsonify({'Description': f"There is no pipeline in that direction",
                                  "Raw_data": f"{error}"}), 404)


@app.route('/RNAseq/', methods=['POST'])
@cross_origin(origin="*")
def launch_nextgene():
    postdata = request.get_json()
    pipeline = 'RNAseq'

    # Create and save configuration
    status_config, config_json_path = generate_configuration(postdata, pipeline)

    # Initialize job
    status_job, dag = launch_job(postdata, config_json_path, pipeline)

    return generate_response(postdata, dag)

if __name__ == "__main__":
    app.run(host="0.0.0.0", use_reloader=True, debug=True)
