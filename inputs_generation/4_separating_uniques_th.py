''' This script separates the unique 29-mers, obtained from a k-mer counting tool KMC3, into
3 files: file with unique 29-mers, originating from CtoT transformations (contain only AGT),
file with unique 29-mers, originating from GtoA transformations (contain only ACT),
file with unique 29-mers, which we cannot be sure about their origin (contain only AT).

Input - path to the file with the unique 29-mers, obtained from a k-mer counting tool KMC3; 
path to the output directory;
output - files with CtoT, GtoA and AT 29-mers in the given output directory.
Usage: python3 4_separating_uniques_th.py /path/to/the/unique/29-mers/file /path/to/the/output/directory 
Can be changed to process output files from other k-mer counting tools. '''

import sys

def main(input_file, output_path):
    '''input_file - path to the file with the unique 29-mers, 
    obtained from a k-mer counting tool KMC3,
    output_path - path to the output directory'''

    # constructing filename of the CtoT file
    ctt_file = ''.join([output_path, '/unique_29mers_ctot'])
    # constructing filename of the GtoA file
    gta_file = ''.join([output_path, '/unique_29mers_gtoa'])
    # constructing filename of the AT file
    at_file = ''.join([output_path, '/unique_29mers_at'])

    with open(input_file, 'r') as fi: 
        with open(ctt_file, 'w') as ctot:  
            with open(gta_file, 'w') as gtoa: 
                with open(at_file, 'w') as at: 
                    # constructing CtoT alphabet
                    ctot_ab = set('AGT')
                    # constructing GtoA alphabet
                    gtoa_ab = set('ACT')
                    # constructing AT alphabet
                    at_ab = set('AT')
                    
                    # for each 29-mer from KMC3 output file
                    for line in fi:
                        # the 29-mer's sequence is in the parts[0]
                        parts = (line.strip()).split('\t')
                        # if all the characters in the 29-mer are from the AT alphabet
                        if set(parts[0]) <= at_ab:
                            # writing the 29-mer to the appropriate output file
                            at.write(f'{parts[0]}\n')
                        # if all the characters in the 29-mer are from the CtoT alphabet
                        elif set(parts[0]) <= ctot_ab:
                            # writing the 29-mer to the appropriate output file
                            ctot.write(f'{parts[0]}\n')
                        # if all the characters in the 29-mer are from the GtoA alphabet
                        elif set(parts[0]) <= gtoa_ab:
                            # writing the 29-mer to the appropriate output file
                            gtoa.write(f'{parts[0]}\n')


if __name__=="__main__":
    # getting the arguments for main function from the command line
    input_file = sys.argv[1]
    output_path = sys.argv[2]
    main(input_file, output_path)
