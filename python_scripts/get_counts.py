import argparse
import logging
import pandas as pd
import python_scripts.python_functions as pf

def transformToRanges(counts):
    """Function to transform the segment of featureCounts to ranges for R
    """
    # Read the file
    df = pd.read_csv(counts, sep='\t', index_col=None)

    # Transform the pseudo-dataframe (cells with more than one field of data) into a normal dataframe.
    new_list = []
    for index, row in df.iterrows():
        chr_list = row.Chr.split(";")
        start_list = row.Start.split(";")
        end_list = row.End.split(";")
        strand_list = row.Strand.split(";")

        for i in range(0,len(chr_list)):
            newline = [row.Geneid, chr_list[i], start_list[i], end_list[i], strand_list[i]]
            new_list.append(newline)

    df2 = pd.DataFrame(new_list)
    df2.columns = ['Geneid','Chr','Start','End','Strand']
    out_counts = counts.replace(".tmpranges.tsv",".ranges.tsv")
    df2.to_csv(out_counts, sep='\t', index=False)

def counts(config, tool_name):
    """Get the counts of a number of bam files in a directory
    """

    bamdir = "/".join(config['tools_conf'][tool_name]['input']['bamdir'].split('/')[0:-1])
    annot = config['tools_conf'][tool_name]['input']['annot']
    output = config['tools_conf'][tool_name]['output']['counts']
    rangestable = output.replace(".tsv",".tmpranges.tsv")

    # Lsit all the bam files in the directory
    filelist = pf.list_files_dir(bamdir, ext = '*.bam')

    # Create the featurecounts command
    tmpoutput = output.replace('.tsv','.tmp.tsv')
    command = f'featureCounts -a {annot} -o {tmpoutput} {" ".join(filelist)}; '
    command += f"cat {tmpoutput} | tail -n +2 | sed -r 's/\t([^\t]+)\//\t/g' | sed 's/.bam//g' | cut --complement -f 2,3,4,5,6 > {output};"
    command += f"cat {tmpoutput} | tail -n +2 | sed -r 's/\t([^\t]+)\//\t/g' | sed 's/.bam//g' | cut -f 1,2,3,4,5 > {rangestable}"

    pf.run_command(command)

    # Create the ranges file. Useful later for the fpkm
    transformToRanges(rangestable)

def get_arguments():
    """
    Function that parse arguments given by the user, returning a dictionary
    that contains all the values.
    """

    # Create the top-level parser
    parser = argparse.ArgumentParser()

    # Mandatory variables
    parser.add_argument('--bamdir', required=True, help='Folder with bam files')
    parser.add_argument('--counts', required=True, help='Table with the counts')
    parser.add_argument('--annot', required=True, help='Annotation file (same than used in mapping)')

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
      "log_files": ["/tmp/full.log"],
      "tools_conf": {
        "get_counts": {
          "input": {
            "bamdir": args.bamdir,
            "annot": args.annot
            },
          "output": {
            "counts": args.counts
            },
          "tool_conf": {
            }
          }
        }
      }

    # Startup the logger format

    counts(config, 'get_counts')


if __name__ == "__main__":
    main()
