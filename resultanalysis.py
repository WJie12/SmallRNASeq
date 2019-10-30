# Filename: resultanalysis.py
#
# -i            Input file path
# -o            Output file path

import re
import sys
import pandas as pd
import matplotlib.pyplot as plt


def getopts(argv):
    opts = {}
    while argv:
        if argv[0][0] == '-':
            opts[argv[0]] = argv[1]
            argv = argv[2:]
        else:
            argv = argv[1:]
    return opts


def analysis_reads(input_path, output_path):
    # load smallRNA DB mapping result, merge
    rna_dbs = ['hg19-tRNAs', 'human_rRNA_5.8S', 'human_rRNA_5S', 'human_rRNA_12S', 'human_rRNA_16S', 'human_rRNA_18S',
               'human_rRNA_28S', 'human_rRNA_45S', 'human_rRNA_other', 'miRBase_21-hsa', 'piR_human', 'Rfam-12.3-human']
    genome = 'genome'
    for idx in range(0, len(rna_dbs)):
        if idx == 0:
            rs_rnadb = pd.read_table(input_path + rna_dbs[idx] + 'nohead.txt', header=None, sep=' ',
                                     names=['name', 'reads', rna_dbs[idx]])
            print(rna_dbs[idx], rs_rnadb.shape)
        else:
            df = pd.read_table(input_path + rna_dbs[idx] + 'nohead.txt', header=None, sep=' ',
                               names=['name', 'reads', rna_dbs[idx]])
            print(rna_dbs[idx], df.shape)
            rs_rnadb = pd.merge(rs_rnadb, df, how='outer', on=['name', 'reads'],
                                suffixes=[rna_dbs[idx - 1], rna_dbs[idx]]).fillna('*')
            print(rs_rnadb.shape)
    # load genome mapping result, merge
    rs_genome = pd.read_table(input_path + genome + 'nohead.txt', header=None, sep=' ',
                              names=['name', 'reads', 'genome'])
    print(rs_genome.shape)

    rs = pd.merge(rs_genome, rs_rnadb, how='outer', on=['name']).fillna('*')
    print(rs.shape)

    rs_group = rs.groupby(rna_dbs)
    rs_csv = rs_group.describe().reset_index()
    rs_csv.to_csv(output_path + 'y_rnadb_y_genome_result.csv', index=False)
    print("--------y_rnadb_y_genome_result.csv finished---------")

    rnadbunmap = rs_group.get_group(('*',) * len(rna_dbs))
    rnadbunmap_group = rnadbunmap.groupby('genome').describe().reset_index()
    rnadbunmap_group.to_csv(output_path + 'n_rnadb_y_genome_result.csv', index=False)
    print("--------n_rnadb_y_genome_result.csv finished---------")
    return 0


def main():
    args = sys.argv[1:]

    try:
        opts = getopts(args)
    except IndexError:
        print("Usage:")
        print(" -i        Input file")
        print(" -o        Output file")
        return 0

    output_path = opts.get("-o")
    if output_path is None:
        print("No output file specified.")
        return -1

    input_path = opts.get("-i")
    if input_path is None:
        print("No input file specified.")
        return -2

    analysis_reads(input_path, output_path)



if __name__ == "__main__":
    main()
