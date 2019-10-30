# Filename: remae.py
#
# -i            Input file
# -o            Output file

import re
import sys

def getopts(argv):
    opts = {}
    while argv:
        if argv[0][0] == '-':
            opts[argv[0]] = argv[1]
            argv = argv[2:]
        else:
            argv = argv[1:]
    return opts


def main():
    args = sys.argv[1:]

    try:
        opts = getopts(args)
    except IndexError:
        print("Usage:")
        print(" -i        Input file")
        print(" -o        Output file")
        return 0

    outputfile = opts.get("-o")
    if outputfile is None:
        print("No output file specified.")
        return -1

    inputfile = opts.get("-i")
    if inputfile is None:
        print("No input file specified.")
        return -2

    # All inputs have been specified at this point, now validate.
    fileRegEx = re.compile(r"^[A-Za-z0-9./\-_\[\]]+$")

    if not fileRegEx.match(outputfile):
        print("Illegal output filename.")
        return -5
    if not fileRegEx.match(inputfile):
        print("Illegal input filename.")
        return -6

    with open(outputfile, 'w') as ofile:
        with open(inputfile) as f:
            idx = 1
            for line in f:
                if line.startswith('@'):
                    index1 = line.rfind('_')
                    if index1 is not -1:
                        line = line[0:index1] + line[index1+1:]
                    index2 = line.rfind('_')
                    if index2 is not -1:
                        #line = line[0] + line[index2+1:-1] + '_' + line[1:index2] + line[-1]
                        line = line[0] + str(idx) + line[index2:]
                    line =  line.split(' ')[0] + line[-1]
                    idx += 1
                    ofile.write(line)
                else:
                    ofile.write(line)
    
    return 0


if __name__ == "__main__":
    main()
