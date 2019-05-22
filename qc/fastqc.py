import argparse
import os
import logging
import zipfile
import glob
import json
import python_scripts.python_functions as pf

def fastqc(config, tool_name):
    """Method to build local parameters for the tool to work

    The method sets the local parameters of the class that are going to be
    needed for the process

    """

    # Get the paths to the files that are going to be input/output
    R1_FILES = config['tools_conf'][tool_name]['input']['fastq_r1']
    R2_FILES = config['tools_conf'][tool_name]['input']['fastq_r2']
    fastqcdir = "/".join(config['tools_conf'][tool_name]['output']['fastqceval'].split('/')[0:-1])
    threads = config['tools_conf'][tool_name]['tool_conf']['threads']

    # Add the fastq files in a file
    FQfilelist = R1_FILES + R2_FILES
    FQfiles = '\n'.join(FQfilelist)
    fof = f"{fastqcdir}/fastq_files.fof"

    f = open(fof, 'w')
    f.write(FQfiles)
    f.close()

    # Construct the command
    full_cmd = ""

    if not os.path.exists(fastqcdir):
        full_cmd += f"mkdir {fastqcdir};"

    full_cmd += f'parallel -j {threads} "fastqc -o {fastqcdir} {{}}" :::: {fof}; '

    pf.run_command(full_cmd)

    fastqc_eval(fastqcdir)

def fastqc_eval(fastqcdir):
    qc_data = {'files': {}}
    estados = []

    pattern = os.path.join(fastqcdir, '*.zip')
    files = glob.glob(pattern)

    for in_file in files:

        zip_ref = zipfile.ZipFile(in_file, 'r')
        zip_ref.extractall(fastqcdir)
        zip_ref.close()


        eval_file = os.path.join(in_file.strip('.zip'), 'summary.txt')

        estado = 'OK'

        for line in open(eval_file, 'r'):
            line = line.strip('\n').split('\t')

            # TODO: Change the ones that fail usually in nextgene
            if line[1] in ['Per base sequence content',
                           'Per sequence GC content',
                           'Sequence Length Distribution',
                           'Overrepresented sequences',
                           'Sequence Duplication Levels',
                           'Kmer Content']:
                pass

            else:
                if line[0] == 'FAIL':
                    estado = 'FAIL'

                if line[0] == 'WARN' and estado != 'FAIL':
                    estado = 'WARN'


        qc_data['files'][in_file.strip('.zip')] = estado
        estados.append(estado)

    qc_data['resume'] = 'PASS'

    if any(estados) == 'WARN':
        qc_data['resume'] = 'WARN'

    if any(estados) == 'FAIL':
        qc_data['resume'] = 'FAIL'

    resfile = f"{fastqcdir}/fastqc.results.txt"
    out_file = open(resfile, 'w')
    json.dump(qc_data, out_file)
    out_file.close()


def get_arguments():
    """
    Function that parse arguments given by the user, returning a dictionary
    that contains all the values.
    """

    # Create the top-level parser
    parser = argparse.ArgumentParser()

    # Mandatory variables
    parser.add_argument('--fastq_r1', '-f1', nargs='*', required=True, help='fastq_r1 file/s')
    parser.add_argument('--fastq_r2', '-f2', nargs='*', required=True, help='fastq_r2 file/s')
    parser.add_argument('--fastqcdir', required=True, help='Folder for fastqc files')

    # Test and debug variables
    parser.add_argument('--dry_run', action='store_true', default=False, help='debug')
    parser.add_argument('--debug', '-d', action='store_true', default=False, help='dry_run')
    parser.add_argument('--test', '-t', action='store_true', default=False, help='test')

    # parse some argument lists
    args = parser.parse_args()

    return args


def main():
    """
    Main function of the script. Launches the rest of the process
    """

    # Get arguments from user input
    args = get_arguments()

    config = {
      "DEBUG": args.debug,
      "TESTING": args.test,
      "DRY_RUN": args.dry_run,
      "outfolder": args.outpath,
      "log_files": ["/tmp/full.log"],
      "tools_conf": {
        "fastqc": {
          "input": {
            "fastq_r1": args.fastq_r1.split(','),
            "fastq_r2": args.fastq_r2.split(',')
            },
          "output": {
            "fastqcdir": args.fastqcdir
            },
          "tool_conf": {
            "threads": "2"
            }
          }
        }
      }

    fastqc(config, 'fastqc')


if __name__ == "__main__":
    main()
