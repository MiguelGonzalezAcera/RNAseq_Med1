import argparse
import os
import json
import glob
import pandas as pd
from statistics import mean

def filterByCell(ser):
    res = []
    for index, val in ser.items():
        res.append(mean(map(int,val.split(','))) > 20)

    return(pd.Series(res, index=ser.index))

def filter_rmats(rmats):
    """"""
    df = pd.read_csv(rmats, sep='\t', index_col=None)
    df = df[df['PValue']<0.05]
    df = df[(df["IncLevelDifference"] < -0.2) | (df["IncLevelDifference"] > 0.2)]

    df = df[(filterByCell(df['IJC_SAMPLE_1'])) | (filterByCell(df['SJC_SAMPLE_1']))]
    df = df[(filterByCell(df['IJC_SAMPLE_2'])) | (filterByCell(df['SJC_SAMPLE_2']))]
    df['chr'] = df['chr'].str.replace("chr","")

    df.to_csv(rmats.replace(".txt","_filt.txt"), sep='\t', index=False)


def rmats(config, mappingtouched, annot, splicetouched):
    """Get the counts of a number of bam files in a directory
    """

    samples = config['comparisons']
    organism = config['options']['organism']
    bam_path = "/".join(mappingtouched.split('/')[0:-1])
    out_dir = "/".join(splicetouched.split('/')[0:-1])

    if not os.path.exists(out_dir):
        out_dir_command = "mkdir {0}".format(out_dir)
        os.system(out_dir_command)

    # Get the files with the bamfiles for each treatment
    for control in samples:
        sample_ids = samples[control].split(",")
        for sample in sample_ids:
            design = pd.read_csv(config['design'], sep='\t', index_col=0)

            folder_command = "mkdir {0}/{1}_{2}_splice".format(out_dir,sample,control)
            os.system(folder_command)

            if not os.path.exists("{0}/{1}_bamfof.txt".format(out_dir,control)):
                control_list = design[design['Tr1'] == control].index.tolist()

                control_list_wpath = []
                for sampleID in control_list:
                    print("{0}/{1}*.bam".format(bam_path,sampleID))
                    control_list_wpath.append(glob.glob("{0}/{1}*.bam".format(bam_path,sampleID))[0])
                control_bam = ",".join(control_list_wpath)

                control_file = open("{0}/{1}_bamfof.txt".format(out_dir,control), "w")
                control_file.write(control_bam)
                control_file.close()
            else:
                control_bam = [line.rstrip('\n') for line in open("{0}/{1}_bamfof.txt".format(out_dir,control))][0]

            if not os.path.exists("{0}/{1}_bamfof.txt".format(out_dir,sample)):
                sample_list = design[design['Tr1'] == sample].index.tolist()

                sample_list_wpath = []
                for sampleID in sample_list:
                    sample_list_wpath.append(glob.glob("{0}/{1}*.bam".format(bam_path,sampleID))[0])
                sample_bam = ",".join(sample_list_wpath)

                sample_file = open("{0}/{1}_bamfof.txt".format(out_dir,sample), "w")
                sample_file.write(sample_bam)
                sample_file.close()
            else:
                sample_bam = [line.rstrip('\n') for line in open("{0}/{1}_bamfof.txt".format(out_dir,sample))][0]


            rmats_command = "python /SOFTWARE/tools/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py --b1 {0}/{1}_bamfof.txt --b2 {0}/{2}_bamfof.txt --gtf {3} --od {0}/{2}_{1}_splice -t paired --nthread 10 --readLength 150".format(out_dir,control,sample,annot)
            print(rmats_command)
            os.system(rmats_command)

            events = ['SE','A5SS','A3SS','MXE','RI']
            for event in events:
                filter_rmats("{0}/{2}_{1}_splice/{3}.MATS.JCEC.txt".format(out_dir,control,sample,event))
                rmatsplot_command = "rmats2sashimiplot --b1 {0} --b2 {1} -e {2}/{4}_{3}_splice/{5}.MATS.JCEC_filt.txt --l1 {3} --l2 {4} -t {5} --exon_s 1 --intron_s 5 -o {2}/{4}_{3}_splice/{5}_viz".format(control_bam,sample_bam,out_dir,control,sample,event)
                print(rmatsplot_command)
                os.system(rmatsplot_command)

    finish_command = "touch {0}".format(splicetouched)

    os.system(finish_command)


def get_arguments():
    """
    Function that parse arguments given by the user, returning a dictionary
    that contains all the values.
    """

    # Create the top-level parser
    parser = argparse.ArgumentParser()

    # Mandatory variables
    parser.add_argument('--annot', required=True, help='annotation file in gtf')
    parser.add_argument('--mappingtouched', required=True, help='Input control file, mapping paths')
    parser.add_argument('--splicetouched', required=True, help='Out control file')
    parser.add_argument('--config', required=True, help='Config dict file')

    # parse some argument lists
    args = parser.parse_args()

    return args


def main():
    """
    Main function of the script. Launches the rest of the process
    """

    # Get arguments from user input
    args = get_arguments()

    with open(args.config) as json_file:
        config = json.load(json_file)

    rmats(config, args.mappingtouched, args.annot, args.splicetouched)


if __name__ == "__main__":
    main()
