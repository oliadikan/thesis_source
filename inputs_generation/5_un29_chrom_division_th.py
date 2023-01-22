''' This script performs division of the unique 29-mers according to their presence in 
the transformations sequences of a certain chromosome for all chromosomes in one run.

To execute it this way, your input forward transformations sequences files have to be placed 
in the same directory, have the same name, altering only in the chromosome number, 
and end with "N_CtoT.fasta"/"N_GtoA.fasta", where N is the number of chromosome (or X/Y); 
your input reverse complements to the forward transformations sequences files have to be placed 
in the same directory, have the same name, altering only in the chromosome number, 
and end with "N_CtoT_rev.fasta"/"N_GtoA_rev.fasta", where N is the number of chromosome (or X/Y);
your files with separated unique 29-mers have to be placed in the same directory, and be named
"unique_29mers_ctot" (for CtoT ones), "unique_29mers_gtoa" (for GtoA ones) and "unique_29mers_at" (for AT ones).

Input - path to the chromosomes forward transformations fasta sequences, up to "N_CtoT.fasta"/"N_GtoA.fasta" part; 
path to the chromosomes reverse complement to the forward transformations fasta sequences, 
up to "N_CtoT_rev.fasta"/"N_GtoA_rev.fasta" part; 
path to the directory with separated unique 29-mers files;
path to the output directory;
output - files of unique 29-mers for each transformation version of each chromosome  
in the given output directory.
Usage: python3 5_un29_chrom_division_th.py /path/to/forward/transformations/before/N_CtoT(GtoA).fasta/part /path/to/reverse/transformations/before/N_CtoT(GtoA)_rev.fasta/part /path/to/the/separated/uniques/directory /path/to/the/output/directory
Can be changed to perform division for just selected chromosomes in a run. '''

import sys
import datetime

def _fasta_seqs_from_filelike(f, COMMENT=b';'[0], HEADER=b'>'[0]):
    ''' fast sequences reading from fasta file, 
    author: Prof. Dr. Sven Rahmann '''
    strip = bytes.strip
    header = seq = False
    for line in f:
        line = strip(line)
        if len(line) == 0:
            continue
        if line[0] == COMMENT:
            continue
        if line[0] == HEADER:
            if header:
                yield seq
            header = True
            seq = bytearray()
            continue
        seq.extend(line)
    yield seq


def main(dir_seq_path, rev_seq_path, uniques_path, output_path):
    '''dir_seq_path - path to the forward sequences of transformations, 
    up to "N_CtoT(GtoA).fasta" part;
    rev_seq_path - path to the sequences of the reverse complements to the transformations, 
    up to "N_CtoT(GtoA)_rev.fasta" part;
    uniques_path - path to the directory with files of separated unique 29-mers;
    output_path - path to the output directory'''

    # list with all chromosomes' numbers
    chromo = ["Y", "22", "21", "20", "19", "18", "17", "16", "15", "14", "13", "12", "11", "10", "9", "8", "X", "7", "6", "5", "4", "3", "2", "1"]

    # for each chromosome
    for chr in chromo:
        # log message
        print('\n' + chr)
        ct = datetime.datetime.now()
        print(ct)    

        # constructing filenames of files with transformations sequences
        cttd = "".join([dir_seq_path, chr, "_CtoT.fasta"]) 
        gtad = "".join([dir_seq_path, chr, "_GtoA.fasta"]) 

        cttr = "".join([rev_seq_path, chr, "_CtoT_rev.fasta"]) 
        gtar = "".join([rev_seq_path, chr, "_GtoA_rev.fasta"]) 

        # constructing filenames of separated unique 29-mers files
        uni_ctt_file = "".join([uniques_path, "/unique_29mers_ctot"])
        uni_gta_file = "".join([uniques_path, "/unique_29mers_gtoa"]) 
        uni_at_file = "".join([uniques_path, "/unique_29mers_at"])

        # constructing filenames of output files
        cttdo = "".join([output_path, "/unique_29mers_chr", chr, "_ctot"]) 
        gtado = "".join([output_path, "/unique_29mers_chr", chr, "_gtoa"]) 

        cttro = "".join([output_path, "/unique_29mers_chr", chr, "_ctot_rev"])
        gtaro = "".join([output_path, "/unique_29mers_chr", chr, "_gtoa_rev"]) 

        # lists of input and output filenames
        seqs = [cttd, gtad, cttr, gtar]
        outs = [cttdo, gtado, cttro, gtaro]

        # for each pair of input-output files
        for i in range(4):

            # reading the chromosome transformation sequence
            seq = ''

            with open(seqs[i], "rb") as f: 
                seq = next(_fasta_seqs_from_filelike(f))
                seq = seq.decode("ASCII")

            # log message
            print("sequence read")
            ct = datetime.datetime.now()
            print(ct)   
                   
            # getting the dictionary of unique 29-mers of the sequence (kmer : position)
            kmers = dict()
            for kl in range(0, len(seq)-28):
                kmer = seq[kl:kl+29]
                if kmer not in kmers.keys():
                    kmers[kmer] = kl

            # log message
            print('len kmers in seq', len(kmers))
            ct = datetime.datetime.now()
            print(ct) 

            # setting the working file of unique 29-mers, 
            # based on the current transformation version
            cur_uni_file = ''
            if i in [0, 3]:
                cur_uni_file = uni_ctt_file
            elif i in [1, 2]:
                cur_uni_file = uni_gta_file
          
            # writing the output file with all the 29-mers from the working file of unique 29-mers 
            # and from the file of AT unique 29-mers, which are among the unique 29-mers of the 
            # current sequence transformation 
            with open(outs[i], 'w') as fo: 
                with open(cur_uni_file, 'rt') as f: 
                    for km in f:
                        km = km.strip()
                        if km in kmers.keys():
                            fo.write(f'{km}\n')
                with open(uni_at_file, 'rt') as fi:
                    for km in fi:
                        km = km.strip()
                        if km in kmers.keys():
                            fo.write(f'{km}\n')

            # log message
            print("output file written")
            ct = datetime.datetime.now()
            print(ct)  
    
        
if __name__=="__main__":
    # getting the arguments for main function from the command line
    dir_seq_path = sys.argv[1]
    rev_seq_path = sys.argv[2]
    uniques_path = sys.argv[3]
    output_path = sys.argv[4]
    main(dir_seq_path, rev_seq_path, uniques_path, output_path)
