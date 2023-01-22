''' This script performs the CtoT and GtoA transfromations of the reference genome sequences of
all chromosomes in one run.

To execute it this way, your input sequences files have to be placed in the same directory,
have the same name, altering only in the chromosome number,
and end with "N.fasta", where N is the number of chromosome (or X/Y).

Input - path to the chromosomes fasta sequences, up to "N.fasta" part; path to the output directory;
output - files of CtoT and GtoA transformations of each chromosome sequence in the given output directory.
Usage: python3 3_direct_transforms_th.py /path/to/input/sequences/before/N.fasta/part /path/to/the/output/directory
Can be changed to perform transformations for just selected chromosomes in a run. '''

import sys

def main(input_file, output_path):
    '''input_file - path to the chromosomes fasta sequences, up to "N.fasta" part,
    output_path - path to the output directory'''

    # list with all chromosomes' numbers
    chroms = ['Y', '22', '21', '20', '19', '18', '17', '16', '15', '14', '13', '12', '11', '10', '9', '8', 'X', '7', '6', '5', '4', '3', '2', '1']

    # to make the names of output files match the name of their respective original sequence files,
    # we are getting the filename of original sequence file from the input_file variable
    prts = input_file.split("/")
    file = prts[(len(prts)-1)]

    # for each chromosome
    for chrom in chroms:
        # constructing the input filename
        inf = ''.join([input_file, chrom, '.fasta'])
        # constructing the CtoT transformation's filename
        ctf = ''.join([output_path, '/' , file, chrom, '_CtoT.fasta'])
        # constructing the GtoA transformation's filename
        gaf = ''.join([output_path, '/', file, chrom, '_GtoA.fasta'])

        with open(inf, 'r') as fasta:
             with open(ctf, 'w') as ct:
                 with open(gaf, 'w') as ga:
                    line = fasta.readline()
                    # for each line of the input file
                    while line != '':
                        line_ct = line
                        line_ga = line
                        # if the current line is a line of a sequence
                        if line[0] != '>':
                            # replacing all Cs with Ts for CtoT
                            line_ct = line_ct.replace('C', 'T')
                            # replacing all Gs with As for GtoA
                            line_ga = line_ga.replace('G', 'A')
                        ct.write(line_ct)
                        ga.write(line_ga)
                        line = fasta.readline()
        # log message
        log = ''.join(["Files for chromosome ", chrom, " written.\n"])
        print(log) 


if __name__=="__main__":
    # getting the arguments for main function from the command line
    input_file = sys.argv[1]
    output_path = sys.argv[2]
    main(input_file, output_path)
