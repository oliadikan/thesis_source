''' This script separates the full reference genome file into separate chromosomes files 
and converts the sequence to the upper case at the same time

Usage: python3 1_chrom_separation_upper_th.py /full/path/to/the/input/file /path/to/the/output/directory '''

import sys

def main(input_file, output_path):
    '''input_file - full path to the file with the full reference genome sequence, 
    output_path - path to the folder, where the separate chromosomes files will be created'''
    # list of chromosomes' names
    chromo = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]
    # list of chromosomes' primary assembly ids
    nc_add = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24"]
    with open(input_file, 'r') as full:
        # for each chromosome
        for i in range(24):
            # constructing path to the output chromosome fasta
            out_file = ''.join([output_path, '/Chr', chromo[i], '.fasta'])
            # constructing the chromosome's primary assembly id
            nc = ''.join(['NC_0000', nc_add[i]])
            with open(out_file, 'w') as chrom: 
                line = full.readline()
                # chromosome fasta header
                chrm = ''
                # chromosome sequence
                seq = ''
                while line != '':
                    # if fasta header
                    if line[0] == '>':
                        # if it is the chromosome's primary assembly header 
                        if chrm != '' and nc in chrm: 
                            # convert all the sequence symbols to upper case
                            seq = seq.upper()
                            # writing chromosome's primary assembly header to the output file
                            chrom.write(chrm)
                            # writing chromosome's sequence to the output file
                            chrom.write(seq)
                        # changing the header
                        chrm = line
                        # emptying the sequence
                        seq = ''
                    # if part of the sequence
                    else:
                        # appending to the sequence
                        seq += line
                    line = full.readline()
                # processing the last part of the input file
                if chrm != '' and nc in chrm: 
                            seq = seq.upper()
                            chrom.write(chrm)
                            chrom.write(seq)
            # log message
            log = ''.join(["Chromosome ", chromo, " written.\n"])
            print(log)       

if __name__=="__main__":
    # getting the arguments for main function from the command line
    input_file = sys.argv[1]
    output_path = sys.argv[2]
    main(input_file, output_path)
