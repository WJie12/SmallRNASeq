# Filename: unmapanalysis.py
#
# -i            Input file path
# -o            Output file path

import re
import sys
import sys
import pandas as pd
import numpy as np

def getopts(argv):
    opts = {}
    while argv:
        if argv[0][0] == '-':
            opts[argv[0]] = argv[1]
            argv = argv[2:]
        else:
            argv = argv[1:]
    return opts


def tofa(input,output):
    with open(output, 'w') as ofile:
        with open(input) as f:
            idx = 1
            for line in f:

                index1 = line.rfind(',')
                if index1 is not -1:
                    line = '>' + str(idx) + '_' + line[index1+1:]+line[0:index1]+line[-1]
                ofile.write(line)
                idx += 1


def analysis_unmap_reads(input_path, output_path):
    df = pd.read_csv(input_path + 'unmap.fastq', sep='_', header=None)
    df = pd.DataFrame(df.values.reshape(df.shape[0]//4, 8))[[1, 2]]

    df = df.drop_duplicates([1, 2], keep="last")

    df_group = df.groupby([2])
    num_reads = len(df_group)
    size_reads = df_group.size()
    sort_reads= size_reads.sort_values(ascending=False)
    sort_reads.to_csv(output_path + 'n_rnadb_n_genome_result.csv', index=True)


def main():
    args = sys.argv[1:]

    try:
        opts = getopts(args)
    except IndexError:
        print("Usage:")
        print(" -i        Input file")
        print(" -o        Output file")
        return 0

    output_path= opts.get("-o")
    if output_path is None:
        print("No output file specified.")
        return -1

    input_path = opts.get("-i")
    if input_path is None:
        print("No input file specified.")
        return -2

    analysis_unmap_reads(input_path, output_path)

    tofa(output_path + 'n_rnadb_n_genome_result.csv', output_path + 'unmap.fa')


if __name__ == "__main__":
    main()
