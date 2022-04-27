# thesis_source
 Source files created during the Master's Thesis (potentially with raw data files)
 
 - As an initial raw data file (reference human genome) is used GRCh38.p14 (from NCBI), it is preprocessed by parser_th.py to keep only the Primary Assembly sequences and remove Ns from them
 
 - On the next step a list of CpGs is created using list_of_cpgs_th.ipynb (There are cells for creating overall/chromosome-based splitted lists of CpGs with inside chromosome/throughout genome numeration and a list of start-end indices of throughout genome numeration for each chromosome)
 
 - For counting k-mers part gerbil (uni-halle tool) is requred - produces binaries (which can be converted to fasta) and counts histograms (which are later improved with a unique counts/genome's length statistic with ucgl_stats_th.py)
