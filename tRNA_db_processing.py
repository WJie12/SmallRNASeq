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


def tRNA_db_processing(input,output):
    with open(output, 'w') as ofile:
        with open(input) as f:
            for line in f:
                if line.startswith('>') and line.rfind("His") != -1:
                    line = line+"G"
                elif line.startswith('>'):
                    line = line
                else:
                    line = line[:-1] + "CCA" + line[-1]
                ofile.write(line)


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

    tRNA_db_processing(input_path, output_path)


if __name__ == "__main__":
    main()
