''' This script forms a .bedGraph file from the file with our method's results.

Input - file with our method's results; output - .bedGraph file.

Usage: python3 stats_to_bedgraph_th.py /path/to/the/method's/result/file '''

import sys

def main(input_path):
    ''' input_path - full path to the file with our method's results '''
    with open(input_path, 'r') as fi:
        # skipping the header
        next(fi)
        # constructing filename for the output file
        out = ''.join([input_path[:-3], 'bedGraph'])
        with open(out, 'w') as fo:
            for line in fi:
                l = (line.strip()).split('\t')
                e = int(l[2])
                e += 1
                if (int(l[3])+int(l[4])) == 0:
                    p = 0
                else:
                    p = round( 100 * (int(l[3]) / (int(l[3])+int(l[4]))) )
                lo = '\t'.join([l[0], l[1], str(e), l[3], l[4], str(p), '\n'])
                fo.write(lo)

if __name__=="__main__":
    # getting the arguments for main function from the command line
    input_path = sys.argv[1]
    main(input_path)
