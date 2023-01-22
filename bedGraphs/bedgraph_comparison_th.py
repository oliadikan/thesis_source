''' This script forms three files from the .bedGraph files with the results 
from our and another tool for further comparison and anlysis of the results. 

Input - .bedGraph file with our method's results 
and a .bedGraph file with the results from another tool; 
output - .tsv files with data on the CpG positions, present in only one of the input files,
present in both files and file combined from the former two .

Usage: python3 bedgraph_comparison_th.py /path/to/our/bedGraph /path/to/bedGraph/from/another/tool /path/to/the/output/file/with/CpGs/present/in/one /path/to/the/output/file/with/CpGs/present/in/both /path/to/the/combined/output/file '''

import sys
import pandas as pd

def main(my_i, com_i, out1, out2, out3):
    ''' my_i - path to our .bedGraph, com_i - path to a .bedGraph from another tool, 
    out1 - path to the output file with CpGs present in only one of the input files, 
    out2 - path to the output file with CpGs present in both input files, 
    out3 - path to the combined output file. '''
    # load with pandas
    my_df = pd.read_csv(my_i, sep='\t', header=None)
    my_df = my_df.drop(my_df.columns[[6]], axis=1) 

    # set column names
    header = ['chr','start','end','meth','unmeth','proc']
    my_df.columns = header[:len(my_df.columns)]

    # load with pandas and set column names
    com_df = pd.read_csv(com_i, sep='\t', header=None)
    com_df.columns = header[:len(com_df.columns)]

    print("Shape of our resulting .bedGraph")
    print(my_df.shape)
    print("Head of our resulting .bedGraph")
    print(my_df.head())
    print("Shape of resulting .bedGraph from another tool")
    print(com_df.shape)
    print("Head of resulting .bedGraph from another tool")
    print(com_df.head())

    # cutting unfound CpGs from our results
    my_df_mnfl = my_df.loc[~((my_df['meth'] == 0) & (my_df['unmeth'] == 0))]
    
    print("Shape of our resulting .bedGraph without unfound CpGs")
    print(my_df_mnfl.shape)
    print("Head of our resulting .bedGraph without unfound CpGs")
    print(my_df_mnfl.head())

    # setting file ids (where the CpG position is from)
    my_df_mnfl['file'] = 'my'
    com_df['file'] = 'conv'

    # getting the CpG positions present in our file, but not in the file from another tool
    either_mdf = my_df_mnfl.loc[~(my_df_mnfl['start'].isin(com_df['start']))]
    
    print("Shape of the dataframe with the CpG positions present in our file, but not in the file from another tool")
    print(either_mdf.shape)
    print("Head of the dataframe with the CpG positions present in our file, but not in the file from another tool")
    print(either_mdf.head())

    # getting the CpG positions present in the file from another tool, but not in our file 
    either_cdf = com_df.loc[~(com_df['start'].isin(my_df_mnfl['start']))]
    
    print("Shape of the dataframe with the CpG positions present in the file from another tool, but not in our file")
    print(either_cdf.shape)
    print("Head of the dataframe with the CpG positions present in the file from another tool, but not in our file")
    print(either_cdf.head())

    # getting the CpG positions present in only one of the input files and sorting those 
    either = pd.concat([either_mdf, either_cdf])
    either = either.sort_values('start')
    
    print("Shape of the dataframe with the CpG positions present in only one of the input files")
    print(either.shape)
    print("Head of the dataframe with the CpG positions present in only one of the input files")
    print(either.head())

    # writing the output file with CpGs present in only one of the input files
    either.to_csv(out1, header=True, index=None, sep='\t', mode='w')

    # getting the CpGs from our file that are present in the file from another tool as well
    my_df_in_com = my_df_mnfl.loc[my_df_mnfl['start'].isin(com_df['start'])]
    # getting the CpGs from the file from another tool that are present in our file as well
    com_df_in_my = com_df.loc[com_df['start'].isin(my_df_mnfl['start'])]

    print("Shape of the dataframe with the CpG positions from our file that are present in the file from another tool as well")
    print(my_df_in_com.shape)

    # getting the CpGs present in both input files
    result = pd.merge(com_df_in_my, my_df_in_com, how='inner', on='start')

    # small reformatting
    result = result.drop(result.columns[[7, 8]], axis=1) 
    result = result.rename(columns={"chr_x": "chr", "end_x": "end"})
    
    print("Shape of the dataframe with the CpG positions present in both input files")
    print(result.shape)
    print("Head of the dataframe with the CpG positions present in both input files")
    print(result.head())

    # writing the output file with CpGs present in both input files
    result.to_csv(out2, header=True, index=None, sep='\t', mode='w')

    # getting a dataframe with all the CpGs present (including small reformatting)
    either = either.rename(columns={"meth": "meth_x", "unmeth": "unmeth_x", "proc": "proc_x", "file": "file_x"})
    allt = pd.concat([result, either])
    allt['file_y'] = allt['file_y'].fillna('none')
    allt = allt.fillna(0)
    allt[['meth_y', 'unmeth_y', 'proc_y']] = allt[['meth_y', 'unmeth_y', 'proc_y']].astype(int)
    allt = allt.sort_values('start')
    
    print("Shape of the combined dataframe")
    print(allt.shape)
    print("Head of the combined dataframe")
    print(allt.head())

    # writing the combined output file
    allt.to_csv(out3, header=True, index=None, sep='\t', mode='w')

if __name__=="__main__":
    # getting the arguments for main function from the command line
    my_i = sys.argv[1]
    com_i = sys.argv[2]
    out1 = sys.argv[3]
    out2 = sys.argv[4]
    out3 = sys.argv[5]
    main(my_i, com_i, out1, out2, out3)
