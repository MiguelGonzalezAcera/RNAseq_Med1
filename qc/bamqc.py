import json
import pysam
import pandas as pd
import os
import python_scripts.python_functions as pf

def bamqc(config, tool_name):
    """
    Method that generates the command to run the mapping of a list of sequences
    R1 and R2 against the reference genome. Because it is a command run by
    CLI, there's the possibility of adding different softwares to map
    the sequences.

    Parameters
    ----------
    None : None

    Returns
    -------
    cmd : str
        Command that indicates the process of mapping the sequences against
        the reference genome.

    """

    # Get the paths to the files that are going to be input/output
    BAMDIR = "/".join(config['tools_conf'][tool_name]['input']['bamdir'].split('/')[0:-1])
    CPBDIR = "/".join(config['tools_conf'][tool_name]['input']['cpbtouched'].split('/')[0:-1])
    LIST = config['tools_conf'][tool_name]['input']['list']

    GENOME = config['tools_conf'][tool_name]['tool_conf']['genome']

    OUTDIR = "/".join(config['tools_conf'][tool_name]['output']['bamqctouched'].split('/')[0:-1])
    bamqctouched = config['tools_conf'][tool_name]['output']['bamqctouched']

    BED = config['tools_conf'][tool_name]['input']['bedfile']

    filelist = sorted(pf.list_files_dir(BAMDIR, ext = '*.bam'))
    cpblist = sorted(pf.list_files_dir(CPBDIR, ext = '*.cpb'))

    for i in range(0,len(filelist)):
        # Set the names of the files
        BAM = filelist[i]
        CPB = cpblist[i]
        QCALNOUT = BAM.replace('.bam','.algn.txt')
        QCMTROUT = BAM.replace('.bam','.mtrg.txt')
        COVFILE = BAM.replace('.bam','.coverages.json')

        # Get sample name
        sample = BAM.split("/")[-1].replace(".bam", "")

        # Init command
        command = ''

        # Create outfolder, if it doesnt exist
        if not os.path.exists(OUTDIR):
            command += f"mkdir {OUTDIR};"
        # bedtools qc software to generate information file
        command += f'java -jar /SOFTWARE/bin/picard-2.19.0.jar CollectAlignmentSummaryMetrics VALIDATION_STRINGENCY=SILENT' +\
            f' R={GENOME} I={BAM} O={QCALNOUT}; '

        # BWA mapping software to generate sam file
        command += f'java -jar /SOFTWARE/bin/picard-2.19.0.jar CollectTargetedPcrMetrics VALIDATION_STRINGENCY=SILENT' +\
            f' I={BAM} O={QCMTROUT} TI={LIST} AI={LIST}; '

        command += f"qualimap bamqc --java-mem-size=20G -bam {BAM} -gff {BED} -outdir {OUTDIR}; "
        # rename qualimap result
        command += f"mv {OUTDIR}/genome_results.txt {OUTDIR}/{sample}.genome_results.txt; "
        command += f'touch {bamqctouched}; '

        value_dict = {}

        # Read file with coverage per base
        cov_df = pd.read_csv(CPB, header=None, sep='\t', low_memory=False)

        # cov_df.columns = ['#chr', 'start', 'end', 'gene', 'rank', 'refseq', 'relPos', 'cov']
        cov_df.columns = ['#chr', 'start', 'end', 'gene', 'relPos', 'cov']

        for value in [1, 2, 5, 10, 20, 40, 50, 80, 100, 200, 500, 1000]:
            value_dict[f'PCT_TARGET_BASES_{value}X'] = str(len(cov_df[cov_df['cov'] > value]) * 100 / len(cov_df['cov']))

        value_dict['MEAN_TARGET_COVERAGE'] = str(round(cov_df['cov'].mean(), 2))

        outfile = open(COVFILE, 'w')
        json.dump(value_dict, outfile)
        outfile.close()

        command += f'mv {QCALNOUT} {QCMTROUT} {COVFILE} {OUTDIR}'

        pf.run_command(command)

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
            "bamdir": args.bamdir,
            "cpbtouched": args.cpbdir,
            "list": args.list,
            "bed": args.bed
            },
          "output": {
            "outdir": args.outdir
            },
          "tool_conf": {
            "genome": args.genome
            }
          }
        }
      }

    bamqc(config, 'bamqc')


if __name__ == "__main__":
    main()
