''' This script forms statistics of the amount of CpG positions of a concrete chromosome,
having a certain number of "marker" unique 29-mers in their +-50-window (in the format 
"number_of_unique_29-mers_in_the_+-50-window:number_of_CpGs_having_this_amount_of_unique_29-mers_in_their_+-50-window"),
based on the forward CtoT transformations for all chromosomes in one run.

To execute it this way, your input forward CtoT transformations sequences files have to be placed 
in the same directory, have the same name, altering only in the chromosome number, 
and end with "N_CtoT.fasta", where N is the number of chromosome (or X/Y); 
your files with chromosome-wise CpG positions have to be placed in the same directory, 
have the same name, altering only in the chromosome number, and end with "N.txt", 
where N is the number of chromosome (or X/Y);
your files with chromosome- and transformation-separated unique 29-mers have to be placed 
in the same directory, have the same name, altering only in the chromosome number, 
and end with "N_ctot", where N is the number of chromosome (or X/Y).

Input - path to the chromosomes forward transformations fasta sequences, up to "N_CtoT.fasta" part; 
path to the files with chromosome-wise CpG positions, up to "N.txt" part;
path to the chromosome- and transformation-separated unique 29-mers files, up to "N_ctot" part;
path to the output directory;
output - files in the format 
"number_of_unique_29-mers_in_the_+-50-window:number_of_CpGs_having_this_amount_of_unique_29-mers_in_their_+-50-window"
for forward CtoT transformation of each chromosome in the given output directory.
Usage: python3 unique_per_cpg_stats_th.py /path/to/forward/transformations/before/N_CtoT.fasta/part /path/to/CpGs/before/N.txt/part /path/to/the/separated/uniques/before/N_ctot/part /path/to/the/output/directory
Can be changed to be executed just for selected chromosomes in a run. '''

import sys
import datetime
import numpy as np
from numba import  njit, uint8, int64, uint64, prange
from os import cpu_count
from concurrent.futures import ProcessPoolExecutor, as_completed 
import gc

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

### Block of functions for numerical encoding of the 29-mers
# from dnaencode.py script (hackgap project, url: https://gitlab.com/rahmannlab/hackgap) ###

def _get_table_dna_to_2bits(default=4):
    b = np.full(256, default, dtype=np.uint8)
    b[ 97]=0  # a
    b[ 65]=0  # A
    b[ 99]=1  # c
    b[ 67]=1  # C
    b[103]=2  # g
    b[ 71]=2  # G
    b[116]=3  # t
    b[ 84]=3  # T
    b[117]=3  # u
    b[ 85]=3  # U
    return b
_TABLE_DNA_TO_2BITS = _get_table_dna_to_2bits()

@njit( ###__signature__ void(uint8[:], uint8[:]), 
    locals=dict(i=int64), nogil=True)
def _dna_to_2bits(x, table):
    for i in range(x.size):
        x[i] = table[x[i]]

@njit( ###__signature__ void(uint8[:]),
    nogil=True)
def quick_dna_to_2bits(x):
    for i in range(len(x)):
        x[i] = _TABLE_DNA_TO_2BITS[x[i]]

@njit( ###__signature__ void(uint8[:]),
    parallel=True, nogil=True)
def parallel_dna_to_2bits(x):
    for i in prange(x.size):
        x[i] = _TABLE_DNA_TO_2BITS[x[i]]


def dna_to_2bits(seq, table=_TABLE_DNA_TO_2BITS):
    if type(seq) == bytes:
        xx = np.frombuffer(bytearray(seq), dtype=np.uint8)
    else:
        xx = np.frombuffer(seq, dtype=np.uint8)
    _dna_to_2bits(xx, table)
    return xx


_TABLE_BITS_TO_DNASTR = ["A","C","G","T"]
def qcode_to_dnastr(qcode, q, table=_TABLE_BITS_TO_DNASTR):
    qc = int(qcode)
    return "".join([ table[((qc >> (2*(q-i-1))) & 3)] for i in range(q)])

@njit( ###__signature__ void(uint64, int64, uint8[:], int64), 
    nogil=True, locals=dict(base=uint64))
def write_qcode_to_buffer(qcode, q, buf, start):
    for i in range(q):
        base = qcode >> (2*(q-i-1)) & 3
        buf[start+i] = uint8(base)


def _get_table_2bits_to_dna(default=4):
    b = np.full(256, 35, dtype=np.uint8)  # fill with b'#'
    b[0] = 65
    b[1] = 67
    b[2] = 71
    b[3] = 84
    b[default] = 78
    return b
_TABLE_2BITS_TO_DNA = _get_table_2bits_to_dna()

@njit( ###__signature__ void(uint8[:], int64, int64),
    nogil=True)
def twobits_to_dna_inplace(buf, start=0, end=0):
    if end <= 0:
        end = len(buf) - end
    for i in range(start, end):
        buf[i] = _TABLE_2BITS_TO_DNA[buf[i]]


########### reverse complements and canonical representation ##############

@njit( ###__signature__ void(uint8[:]), 
    nogil=True,
    locals=dict(c1=uint8,c2=uint8,n=int64,drei=uint8))
def revcomp_inplace(seq):
    n = seq.size
    drei = 3
    for i in range((n+1)//2):
        j = n-1-i
        c1 = seq[i]
        c2 = seq[j]
        seq[j] = drei-c1  if c1<4  else c1
        seq[i] = drei-c2  if c2<4  else c2


@njit( ###__signature__ void(uint8[:], uint8[:]),
    nogil=True, locals=dict(
        c=uint8,n=int64,drei=uint8,rc=uint8[:]))
def revcomp_to_buffer(seq, rc):
    n = seq.size
    drei = 3
    for i in range(n):
        c = seq[n-1-i]
        rc[i] = drei - c  if c < 4  else c


@njit( ###__signature__ uint8[:](uint8[:]),
    nogil=True, locals=dict(rc=uint8[:]))
def revcomp(seq):
    rc = np.empty_like(seq, dtype=np.uint8)
    revcomp_to_buffer(seq, rc)
    return rc


@njit( ###__signature__ uint64(uint64, int64),
    nogil=True, locals=dict(
        code=uint64, drei=uint64, rc=uint64, c=uint64))
def revcomp_code(code, q):
    # only works for 0 <= q <= 31 !
    # when using uints, due to a potential bug in numpy/numba,
    # we would have to re-declare code as uint64 locally.
    drei = uint64(3)
    rc = 0
    for i in range(q):
        c = drei - (code & drei)
        rc = (rc << 2) | c
        code >>= 2
    return rc



@njit( ###__signature__ uint64[:](),
    nogil=True, locals=dict(
        code=uint64, c=uint64))
def _get_rctable():
    rctable = np.zeros(256, dtype=np.uint64)
    for c in range(256):
        rctable[c] = revcomp_code(c, 4)
    return rctable
_RCTABLE = _get_rctable()


@njit( ###__signature__ int64(uint64, int64),
    nogil=True, locals=dict(
        code=uint64, rc=uint64, c=uint64))
def revcomp_code_table(code, q):
    rc = 0
    while q >= 4:
        c = _RCTABLE[code & 255]
        rc = (rc << 8) | c
        code >>= 8
        q -= 4
    for i in range(q):
        c = 3 - (code & 3)
        rc = (rc << 2) | c
        code >>= 2
    return rc

@njit( ###__signature__ uint64(uint64, int64),
    nogil=True, locals=dict(
        code=int64, rc=int64))
def canonical_code(code, q):
    rc = revcomp_code(code, q)
    return code  if code <= rc  else rc


def generate_revcomp_and_canonical_code(q, rcmode):
    """
    return pair of functions (revcomp_code_q, canonical_code_q)
    specialized for q-gram codes for the given value of q.
    It is expected that LLVM optimization does loop unrolling.
    """
    @njit( ###__signature__ uint64(uint64),
        nogil=True, locals=dict(
            code=uint64, rc=uint64, c=uint64))
    def _rc(code):
        rc = 0
        t = q // 4
        for i in range(t):
            c = _RCTABLE[code & 255]
            rc = (rc << 8) | c
            code >>= 8
        r = q % 4
        for i in range(r):
            c = 3 - (code & 3)
            rc = (rc << 2) | c
            code >>= 2
        return rc
    
    if rcmode == "min":
        @njit( ###__signature__ uint64(uint64),
            nogil=True,locals=dict(
                code=uint64, rc=uint64))
        def _cc(code):
            rc = _rc(code)
            return code  if code <= rc  else rc
    elif rcmode == "max":
        @njit( ###__signature__ uint64(uint64),
            nogil=True, locals=dict(
                code=uint64, rc=uint64))
        def _cc(code):
            rc = _rc(code)
            return code  if code >= rc  else rc
    elif rcmode == "r":
        @njit( ###__signature__ uint64(uint64),
            nogil=True, locals=dict(
                code=uint64, rc=uint64))
        def _cc(code):
            rc = _rc(code)
            return rc
    else:  # 'f', 'both', ...
        @njit( ###__signature__ uint64(uint64), 
            nogil=True, locals=dict(
                code=uint64, rc=uint64))
        def _cc(code):
            return code

    return _rc, _cc

# translation of a DNA buffer into codes
def make_twobit_to_codes(k, rcmode, invalid=uint64(-1)):
    _, ccc = generate_revcomp_and_canonical_code(k, rcmode)

    @njit(locals=dict(code=uint64))
    def twobit_to_codes(seq, out, start=0, n=-1):
        if n == -1:
            n = len(seq) - k + 1 - start
        for i in range(start, start+n):
            code = 0
            for j in range(k):
                c = seq[i+j]
                if c >= 4:
                    out[i-start] = uint64(invalid)
                    break
                code = (code << 2) | c
            else:
                code = ccc(code)
                out[i-start] = code

    @njit(locals=dict(code=uint64))
    def twobit_to_code(seq, start):
        code = 0
        for j in range(k):
            c = seq[start+j]
            if c >= 4:
                return uint64(invalid)
            code = (code << 2) | c
        else:
            code = ccc(code)
            return code

    return twobit_to_codes, twobit_to_code

### End of block ###

def main(dir_seq_path, cpg_path, uniques_path, output_path):
    ''' dir_seq_path - path to the forward sequences of transformations, up to "N_CtoT.fasta" part;
    cpg_path - path to the files with chromosome's CpG positions, up to "N.txt";
    uniques_path - path to the files of chromosome- and transformation-separated unique 29-mers, 
    up to "N_ctot";
    output_path - path to the output directory '''

    # bits_to_qc converts two-bit representation of a 29-mer to its numerical code,
    # discriminates between one-directional 29-mer and its reverse complement
    fs, bits_to_qc = make_twobit_to_codes(29, "both")

    # list with all chromosomes' numbers 
    chromo = ["Y", "22", "21", "20", "19", "18", "17", "16", "15", "14", "13", "12", "11", "10", "9", "8", "X", "7", "6", "5", "4", "3", "2", "1"]

    # for each chromosome
    for chr in chromo:
        # calling garbage collector
        gc.collect()

        # log message
        print('\n' + chr)
        ct = datetime.datetime.now()
        print(ct)    

        # constructing filename of the CpG positions file
        cpgs_file = "".join([cpg_path, chr, ".txt"]) # (["/scratch/chroms/clean/cpgs/CpGs_", chr, ".txt"])
        # constructing filename of CtoT transformation sequence file
        ctt = "".join([dir_seq_path, chr, "_CtoT.fasta"]) # (["/scratch/chroms/clean/direct/GRCh38.p14_Chr", chr, "_CtoT.fasta"])
        # constructing filename of ctot unique 29-mers file 
        uni_ctt = "".join([uniques_path, chr, "_ctot"]) 

        # forming list of CpG position of the chromosome (forward-directional)
        cpgs = []

        with open(cpgs_file, 'rt') as f:
            for pos in f:
                if not pos.startswith("N"):
                    cpgs.append(int(pos.strip()))

        cpgs = np.array(cpgs, dtype=int)
        cpgs = np.sort(cpgs)
    
        # log message
        print("list of cpgs formed")
        ct = datetime.datetime.now()
        print(ct)   

        gc.collect()
          
        # reading the chromosome transformation sequence 
        seq = ''

        with open(ctt, "rb") as f: 
            seq = next(_fasta_seqs_from_filelike(f))
            seq = seq.decode("ASCII")

        # log message    
        print("sequence read")
        ct = datetime.datetime.now()
        print(ct)   
    
        # Python index of the last position in the sequence
        ln = len(seq) - 1
                   
        # getting the dictionary of unique 29mers of the sequence
        # (kmer : index in forward-directional transformation)
        kmers = dict()
        # for each 29-mer in the sequence
        for kl in range(0, len(seq)-28):
            # bit encoding of the 29-mer
            km = (seq[kl:kl+29]).encode('ASCII') 
            # numerical code of the 29-mer
            kmer = bits_to_qc(dna_to_2bits(km), 0)
            # if the first occurence of the 29-mer in the sequence
            if kmer not in kmers.keys():
                kmers[kmer] = kl 
        
        # log message
        print('len kmers in seq', len(kmers))
        ct = datetime.datetime.now()
        print(ct) 

        # deleting the sequence to free memory resources
        del seq
        gc.collect()

        # getting the list of unique 29-mers in chromosome positions (in_chr)
        in_chr = []

        with open(uni_ctt, 'rt') as f:
            # for each 29-mer in unique 29-mers file
            for km in f:
                # bit encoding of the unique 29-mer
                kmm = (km.strip()).encode('ASCII')
                # numerical code of the unique 29-mer
                kmmm = bits_to_qc(dna_to_2bits(kmm), 0)
                # if the unique 29-mer found in sequence 
                if kmmm in kmers.keys():
                    in_chr.append(kmers[kmmm])

        # deleting kmers to free memory resources
        del kmers
        gc.collect()
    
        # forming np.array from in_chr and sorting it
        in_chr = np.array(in_chr, dtype=int)
        in_chr = np.sort(in_chr)

        # log message
        print('len in_chr', len(in_chr))
        ct = datetime.datetime.now()
        print(ct) 
    
    
        # main dictionary 
        # (number of unique 29-mers in the +-50-window: number of CpGs, 
        # having this amount of unique 29-mers in their +-50-window)
        stat = dict()
    
        # preprocessing for the multiprocessing part

        # starting indecies of the +-50-windows for each CpG
        start = cpgs - 50
        # setting the negative position values to 0
        start[start<0] = 0
        # ending indecies of the +-50-windows for each CpG
        end = cpgs + 50
        # setting the position values outside of the sequence to the index of the last 
        # sequence position
        end[end>(ln)] = ln
        # leaving only those CpG positions, whose +-50-windows are inside 
        # the sequence area, where there are potentially "marker" unique 29-mers
        sub = np.where((start >= in_chr[0]) & (end <= in_chr[len(in_chr)-1]))
        # amount of CpG positions before selection
        in_l = len(cpgs)
        cpgsct = cpgs[sub]
        start = start[sub]
        end = end[sub]
        # amount of CpG positions after selection
        n_l = len(cpgsct)
        # amount of CpGs, already having 0 "marker" unique 29-mers in their +-50-window
        stat[0] = in_l-n_l

        # log message
        print("len cpgs", len(cpgsct))
        print("preprocessing done, start of multiprocessing now")
        ct = datetime.datetime.now()
        print(ct)   
   
        def processing(cpgs, start, end):
            '''function for processing of one chunk of cpg positions.
            cpgs - chunk of CpG positions, start - starting positions of the +-50-windows,
            end - ending positions of the +-50-windows'''
            # stat for the current chunk of CpGs
            curst = dict()
            # for each CpG in the chunk
            for j in range(len(cpgs)):
                # finding amount of all unique 29-mers, which starting positions 
                # are inside the +-50-window of the current CpG 
                # ("marker" unique 29-mers, associated with this CpG) 
                tot = len(in_chr[(in_chr >= start[j]) & (in_chr <= end[j])])
                # increasing the counter of the corresponding amount of 
                # "marker" unique 29-mer in curst
                if tot in curst.keys():
                    curst[tot] += 1
                else:
                    curst[tot] = 1
            return curst
    
        # recommended amount of processes
        cores = cpu_count() - 5 
        # partition of cpgsct, start and end into approximately equal chunks,
        # each chunk is lately processed in its own process
        cpg_parts = np.array_split(cpgsct, cores)
        st_parts = np.array_split(start, cores)
        e_parts = np.array_split(end, cores)


        if __name__ == '__main__':
            # multiprocessing using pool of processes
            with ProcessPoolExecutor(cores) as p:
                futures = {p.submit(processing, cpg_parts[i], st_parts[i], e_parts[i]) for i in range(len(cpg_parts))}
                # for each future in order of completion
                for fut in as_completed(futures):
                    # accumulating data from future result in stat
                    d = fut.result()
                    for k in d.keys():
                        if k in stat.keys():
                            stat[k] += d[k]
                        else:
                            stat[k] = d[k]
                    # log message
                    print("future processed")
                    ct = datetime.datetime.now()
                    print(ct)  
                    # deleting d to free memory resources
                    del d
                    gc.collect()
        
            gc.collect()

        # constructing the filename of the output file
        stats = "".join([output_path, '/', chr, '_ctot_stat.txt']) 
        # writing the output file
        with open(stats, 'w') as f:
            sk = {key:stat[key] for key in sorted(stat.keys())}     
            for s in sk.keys():
                line = "".join([str(s), ':', str(sk[s]),  '\n']) 
                f.write(line)
        
        # log message
        print("stats file written")
        ct = datetime.datetime.now()
        print(ct) 
    
if __name__=="__main__":
    # getting the arguments for main function from the command line
    dir_seq_path = sys.argv[1]
    cpg_path = sys.argv[2]
    uniques_path = sys.argv[3]
    output_path = sys.argv[4]
    main(dir_seq_path, cpg_path, uniques_path, output_path)
