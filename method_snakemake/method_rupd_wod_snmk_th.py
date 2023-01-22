''' The main script of the method, performs methylation calling of WGBS data files of a sample
for one chromosome in a run. All the usage specifications are given in the Snakefile '''

from datetime import datetime
from gzip import open as gzopen
from os import listdir 
from itertools import islice
from concurrent.futures import ProcessPoolExecutor, as_completed
from gc import collect
from argparse import ArgumentParser
import numpy as np 
from numba import  njit, uint8, int64, uint64, prange
import pandas as pd


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


def get_dfl(chrom, il):
    """ parameters: chrom - string, chromosome number, il - string, path to the directory with input lists; 
    return: dfl - np.ndarray(dtype={'names': ('cpg', 'code', 'shift'), 'formats': ('uint32', 'uint64', 'int8')}), contains the data from input files,
    gf - structure with index being a kmer converted to code and value being collection of dfl indices, where 'code' value equals to the gf index """
    f1 = "".join([il, '/', chrom, '_cpgs_kmers_list_upd']) # PATH TO LISTS FOR FORWARD SEQUENCE
    f2 = "".join([il, '/', chrom, '_cpgs_kmers_list_rev']) # PATH TO LISTS FOR REVERSE COMPLEMENT
    dfl1 = np.loadtxt(f1, dtype={'names': ('cpg', 'code', 'shift'), 'formats': ('uint32', 'uint64', 'int8')})
    dfl2 = np.loadtxt(f2, dtype={'names': ('cpg', 'code', 'shift'), 'formats': ('uint32', 'uint64', 'int8')})
    dfl = np.concatenate((dfl1, dfl2), dtype={'names': ('cpg', 'code', 'shift'), 'formats': ('uint32', 'uint64', 'int8')})
    df = pd.DataFrame(dfl)
    gf = df.index.to_series().groupby(df['code']).unique()
    print("nparray out of lists formed") 
    ct = datetime.now()
    print(ct)
    return dfl, gf


def get_dfc(chrom, ic):
    """ parameters: chrom - string, chromosome number, ic - string, path to the files containing cpg positions;
    return: dfc - dictionary for storing statistics of methylated/unmethylated/seen """
    dfc = dict()
    f2 = "".join([ic, chrom, ".txt"]) # PATH TO CPG POSITIONS FILE FOR CURRENT CHROMOSOME
    with open(f2, 'rt') as f:
        for pos in f:
            if not pos.startswith("N"):
                dfc[np.uint32(pos.strip())] = [0, 0, 0] 

    print("dict for cpgs formed")
    ct = datetime.now()
    print(ct)
    return dfc

def get_argument_parser():
    parser = ArgumentParser(prog="29merMethCaller", description="29mer-based methylation state caller for a collection of samples")
    parser.add_argument('-d', '--directory', dest='drt', help='name of the directory containing sample files (required)', required=True)
    parser.add_argument('-c', '--chrom', dest='chrom', choices=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y'], help='chromosome for methylation state calling (required)', required=True)
    parser.add_argument('-p', '--process', dest='proc', default=5, help='number of processes for multiprocessing part (recomendations)', type=int)
    parser.add_argument('-l', '--log', dest='log', help='path for the log file', required=True)
    parser.add_argument('-i', '--lists', dest='il', help='path to the directory with input lists', required=True)
    parser.add_argument('-k', '--cpgs', dest='ic', help='path to the cpg positions file (up to N.txt, N - number of chromosome)', required=True)
    parser.add_argument('-o', '--output', dest='out', help='path for the output file', required=True)
    parser.add_argument('-s', '--sample', dest='sd', help='path to the directory with samples subdirectories', required=True)
    return parser


parser = get_argument_parser()
args = parser.parse_args()

# bits_to_qc converts two-bit representation of a 29-mer to its numerical code,
# discriminates between one-directional 29-mer and its reverse complement
fs, bits_to_qc = make_twobit_to_codes(29, "both")

print("starting")
ct = datetime.now() # TIMESTAMP
print(ct)   

logf = ''.join([args.log, '/', str(ct),'_', args.drt, '_chr', args.chrom, '.log']) # PATH FOR THE LOG FILE

subd =  args.drt # CURRENT SAMPLE SUBDIRECTORY, EXAMPLE: "Sample_43_Hm03_BlMa_Ct2_WGBS_E_EC"
rundir = "".join([args.sd, "/", subd, "/"]) # PATH TO THE CURRENT SAMPLE SUBDIRECTORY 
runfiles = ["".join([rundir, f]) for f in listdir(rundir) if f.endswith('fastq.gz')] # ARRAY OF PATHS FOR SAMPLE FILES INSIDE THE SUBDIRECTORY
runfiles.sort()

###### LOG FILE FOR TRACING THE EXECUTION TIME #############
# BLOCKS with open(logf, 'a') as lf: ARE FOR LOGGING
with open(logf, 'a') as lf:
    line = ''.join([subd, '\n'])
    lf.write(line)


chrom = args.chrom # CURRENT CHROMOSOME
il = args.il # PATH TO THE DIRECTORY WITH INPUT LISTS
ic = args.ic # PATH TO THE DIRECTORY WITH CPG POSITIONS FILES
print(chrom)

with open(logf, 'a') as lf:
    ct = datetime.now()
    line = ' '.join([chrom, str(ct), '\n'])
    lf.write(line)
    
npfl, gf = get_dfl(chrom, il)


dfce = get_dfc(chrom, ic) # THIS ONE IS FOR THE CUMULATIVE STATISTICS OVER THE SUBDIRECTORY -> OUTPUT FILE IS FORMED OUT OF IT
dfcb = get_dfc(chrom, ic) # THIS ONE STAYS CLEAN - IT IS PASSED TO THE FUNCTION, WHICH TAKES PART IN THE MULTIPROCESSING, 
# SO EACH COPY IS UPDATED FOR ONE PARTICULAR PORTION OF READS INSIDE THE collecting_stats FUNCTION

# HERE WE PUT READING FASTQ FILES AND PROCESSING THEM ONE BY ONE

def collecting_stats(reads, dfcb):
    """ parameters: reads - array, reads from a file, dfcb - empty dictionary for storing statistics of methylated/unmethylated/seen;
    return: dfc - dictionary for storing statistics of methylated/unmethylated/seen, contains statistics for current portion of reads """
    dfc = dfcb 
    for rd in reads:
        # rd - ONE READ FROM CURRENT PORTION OF READS
        pr_in_r = [] # ARRAY FOR CPGS, PREVIOUSLY SEEN IN THE READ

        for kl in range(0, (len(rd)-28)):
            # HERE WE ARE LOOKING AT 29MER FROM THE READ, STARTING FROM INDEX kl
            kmerct = ((rd[kl:kl+29]).replace('C', 'T')).encode('ASCII') # CTOT TRANSFORMATION OF CURRENT 29MER
            kmerga = ((rd[kl:kl+29]).replace('G', 'A')).encode('ASCII') # GTOA TRANSFORMATION OF CURRENT 29MER
            ctcode = bits_to_qc(dna_to_2bits(kmerct), 0) # CODE FROM CTOT TRANSFORMATION
            gacode = bits_to_qc(dna_to_2bits(kmerga), 0) # CODE FROM GTOA TRANSFORMATION
            
            if ctcode in gf: # IF ctcode ENCODES A UNIQUE 29MER ASSOCIATED WITH A CPG FROM THE CURRENT CHROMOSOME
                for k in gf[ctcode]: # FOR EACH ENTRY OF npfl, WHERE 'code' = ctcode
                    cpg = npfl['cpg'][k] # CPG POSITION
                    difr = npfl['shift'][k] # SHIFT VALUE
                    if cpg not in pr_in_r: # IF CPG POSITION HASN'T BEEN ALREADY SEEN IN THIS READ
                        ind = int(kl + difr) # SUPPOSEDLY INDEX OF CPG POSITION INSIDE THE READ
                        if ind >= 0 and ind < len(rd): # IF THE READ ACTUALLY COVERS THIS CPG POSITION
                            if rd[ind] == 'C':
                                dfc[cpg][0] += 1 # INCREASING COUNTER FOR METHYLATED
                            elif rd[ind] == 'T':
                                dfc[cpg][1] += 1 # INCREASING COUNTER FOR UNMETHYLATED
                            dfc[cpg][2] += 1 # INCREASING COUNTER FOR SEEN
                            pr_in_r.append(cpg) # PUTTING THE CPG POSITION IN THOSE PREVIOUSLY SEEN IN THE READ

            if gacode in gf: # SAME FOR GTOA TRANSFORMATION 
                for k in gf[gacode]: 
                    cpg = npfl['cpg'][k]  
                    difr = npfl['shift'][k] 
                    if cpg not in pr_in_r:
                        ind = int(kl + difr)
                        if ind >= 0 and ind < len(rd):
                            if rd[ind] == 'G':
                                dfc[cpg][0] += 1
                            elif rd[ind] == 'A':
                                dfc[cpg][1] += 1
                            dfc[cpg][2] += 1
                            pr_in_r.append(cpg)
    return dfc


# Reading fastq files, getting all reads and their transformations 
for runfile in runfiles:

    with open(logf, 'a') as lf:
        ct = datetime.now()
        line = ' '.join([runfile, '\n', 'start', str(ct), '\n'])
        lf.write(line)


    reads = [] # ARRAY FOR THE READS FROM THE runfile
    with gzopen(runfile, 'rt') as rf:
        for line in islice(rf, 1, None, 4):
            line = line.strip()
            # line preprocessing (cut beginning till the last N, cut last symbol)
            rf = line.rfind('N')
            line = line[rf+1:-1] # got rid of the last A and the part in front till the last N
            reads.append(line) 
        
    l = " ".join([runfile, "read"])
    print(l) 
    ct = datetime.now()
    print(ct) 

    with open(logf, 'a') as lf:
        line = ' '.join(['read', str(ct), '\n'])
        lf.write(line)

    reads = np.array(reads, dtype='object') 
    cores = args.proc # AMOUNT OF THE PROCESSES FOR MULTIPROCESSING PART  
    r_parts = np.array_split(reads, cores) # SPLITTING reads INTO PORTIONS FOR MULTIPROCESSING

    with open(logf, 'a') as lf:
        line = ' '.join(['cores', str(cores), '\n'])
        lf.write(line)
    
    print("start of multiprocessing now")
    ct = datetime.now()
    print(ct)   

    with open(logf, 'a') as lf:
        line = ' '.join(['mp', str(ct), '\n'])
        lf.write(line)

    if __name__ == '__main__':
        with ProcessPoolExecutor(cores) as p:
            futures = {p.submit(collecting_stats, r_parts[i], dfcb) for i in range(len(r_parts))}
            for fut in as_completed(futures):
                d = fut.result()
                for k in d.keys():
                    # UPDATING THE COUNTERS OF MAIN STATISTICS DICTIONARY
                    dfce[k][0] += d[k][0]
                    dfce[k][1] += d[k][1]
                    dfce[k][2] += d[k][2]
                
                print("future processed")
                ct = datetime.now()
                print(ct)  
                del d
                collect()
    
    with open(logf, 'a') as lf:
        ct = datetime.now()
        line = ' '.join(['processed', str(ct), '\n'])
        lf.write(line)


print("all files processed")
ct = datetime.now()
print(ct)

with open(logf, 'a') as lf:
    line = ' '.join(['directory done', str(ct), '\n'])
    lf.write(line)

# Forming the file with results
o = args.out # PATH FOR THE SUBDIRECTORY FOR THE OUTPUT FILE

resf = "".join([o, "/", subd, "_", chrom, "_meth_results_upd_rev.tsv"]) 
with open(resf, 'w') as f:
    line = '\t'.join(["chrom", "chromStart", "chromEnd", "meth", "unmeth", "total\n"])
    f.write(line)
    sc = {key:dfce[key] for key in sorted(dfce.keys())}
    for c in sc.keys():
        line = '\t'.join([''.join(['chr', chrom]), str(c), str(c+1), str(sc[c][0]), str(sc[c][1]), str(sc[c][2])]) #  str(c), str(c+1)
        line = "".join([line, '\n'])
        f.write(line)


print("results written, finish")
ct = datetime.now()
print(ct)

with open(logf, 'a') as lf:
        line = ' '.join(['finish', str(ct), '\n'])
        lf.write(line)
