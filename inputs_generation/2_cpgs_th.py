''' This script forms a file of all CpG positions in the chromosome.
It is originally a script from my supervisor Prof. Dr. Sven Rahmann, 
with small alterations from my side.

Input - fasta file with the sequence of the chromosome of interest,
output - file with the list of all CpG positions in the given chromosome sequence.
Usage: python3 2_cpgs_th.py < /path/to/the/chromosome/fasta > /path/to/the/output/file
Can also be used with a gz compressed fasta (.gz):
python3 2_cpgs_th.py < <(gzcat /path/to/the/gz/compressed/chromosome/fasta) > /path/to/the/output/file '''

from sys import stdin
import numpy as np
from numba import njit

# generator that yields (header, sequence) from FASTA (binary)
def fasta_items_from_filelike(f, COMMENT=b';'[0], HEADER=b'>'[0]):
    # f must be an open file in READ BINARY ('rb') mode;
    # it can be sys.stdin.buffer.
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
                yield header, seq
            header = line[1:]
            seq = bytearray()
            continue
        seq.extend(line)
    yield header, seq


@njit
def find_cpg_positions(seq, pos):
    i = 0
    n = len(seq)
    for p in range(n-1):
        sp = seq[p]  # must be C (67) or c (99)
        sq = seq[p+1]  # must be G (71) or g (103)
        if (sp == 67 or sp == 99) and (sq == 71 or sq == 103):
            pos[i] = p
            i += 1
            if i >= len(pos):
                raise ValueError("pos buffer too small for all positions")
    return i


cpgpos = np.zeros(3_000_000_000, dtype=np.uint32)
for (header, seq) in fasta_items_from_filelike(stdin.buffer):
    ch = header.split()[0]
    print(ch.decode(), len(seq))  # name and length of chromosome
    ncpgs = find_cpg_positions(seq, cpgpos)
    for p in cpgpos[:ncpgs]:
        print(f"{p}")
