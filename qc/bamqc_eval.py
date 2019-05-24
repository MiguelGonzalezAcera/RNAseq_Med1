import json
import pysam
import pandas as pd
import os
import python_scripts.python_functions as pf

def bamqcEval(config, tool_name):
    """
    """
    # Get paths of data
    bamqcpath = "/".join(config['tools_conf'][tool_name]['input']['bamqctouched'].split('/')[0:-1])

    # Get the list of files
    algnfilelist = sorted(pf.list_files_dir(bamqcpath, ext = '*.algn.txt'))
    mtrgfilelist = sorted(pf.list_files_dir(bamqcpath, ext = '*.mtrg.txt'))
    gresfilelist = sorted(pf.list_files_dir(bamqcpath, ext = '*.genome_results.txt'))
    cvrgfilelist = sorted(pf.list_files_dir(bamqcpath, ext = '*.coverages.json'))

    OUTEVAL = config['tools_conf'][tool_name]['output']['evaloutfile']

    # Init the dictionary
    qc_data = {}

    for i in range(0, len(algnfilelist)):
        # Assign all of the names of the files we are going to analyze
        ALIGNQCBAM = algnfilelist[i]
        print(ALIGNQCBAM)
        METRQCBAM = mtrgfilelist[i]
        QUALIMAP = gresfilelist[i]
        COVERAGE = cvrgfilelist[i]

        # Get sample name
        sampleid = algnfilelist[i].split("/")[-1]

        qc_data[f'data_{sampleid}'] = {}

        bam_qc = open(ALIGNQCBAM, 'r')

        qc_df = pd.DataFrame()
        c = 0

        for line in bam_qc:
            if line.startswith('#') or line == '\n':
                continue

            line = line.strip('\n').split('\t')
            qc_df = pd.concat([qc_df, pd.Series(line, name=c)], axis=1)
            c += 1

        bam_qc_2 = open(METRQCBAM, 'r')
        qc_df2 = pd.DataFrame()
        c = 0
        cont = True

        for line in bam_qc_2:
            if line.startswith('#') or line == '\n' or cont is False:
                continue

            line = line.strip('\n').split('\t')
            if line[0] == 'coverage_or_base_quality':
                cont = False
                continue

            qc_df2 = pd.concat([qc_df2, pd.Series(line, name=c)], axis=1)
            c += 1

        for idx, row in qc_df2.iterrows():
            qc_data[f'data_{sampleid}'][row[0]] = row[1]


        with open(COVERAGE) as json_file:
            covdict = json.load(json_file)

        qc_data[f'data_{sampleid}'].update(covdict)

        qc_data['resume'] = 'PASS'

        value = "0"
        for line in open(QUALIMAP, 'r'):
            if line.startswith('     median insert size'):
                value = line.strip('\n').replace('     median insert size = ', '')

        qc_data[f'data_{sampleid}']['MEDIAN_INSERT_SIZE'] = value

        bam_qc.close()
        bam_qc_2.close()

    # Save the dict to files
    outfile = open(OUTEVAL, 'w')
    json.dump(qc_data, outfile)
    outfile.close()

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
      "log_files": ["/tmp/full.log"],
      "tools_conf": {
        "bamqc": {
          "input": {
            "bamqctouched": args.bamdir,
            "cpbtouched": args.cpbdir
            },
          "output": {
            "bameval": args.outfile
            },
          "tool_conf": {
            }
          }
        }
      }

    bamqcEval(config, 'bamqceval')


if __name__ == "__main__":
    main()
