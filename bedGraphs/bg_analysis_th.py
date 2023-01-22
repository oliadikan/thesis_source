''' This script produces a pie chart .png with the statistics described in '4.2 Results comparison'
for the combined .tsv output file from bedgraph_comparison_th.py.

Input - combined .tsv file from bedgraph_comparison_th.py; 
output - pie chart .png with the statistics described in '4.2 Results comparison'.

Usage: python3 bg_analysis_th.py /path/to/the/combined/file '''

import sys
import pandas as pd
import matplotlib.pyplot as plt

def main(res_file):
    ''' res_file - path to the combined file from bedgraph_comparison_th.py '''
    # load dataframe from combined .tsv
    df = pd.read_csv(res_file, sep='\t')
    # total amount of CpGs
    tot = len(df.index)
    print("Total amount of CpGs")
    print(tot)

    # CpGs present in only one of the files
    either_df = df.loc[df['file_y'] == 'none']
    either = len(either_df.index)
    print("Amount of CpGs present in only one of the files")
    print(either)

    # CpGs with fully matching results across the resulting files
    full_match_df = df.loc[(df['meth_x'] == df['meth_y']) & (df['unmeth_x'] == df['unmeth_y'])]
    full_match = len(full_match_df.index)
    print("Amount of CpGs with fully matching results across the resulting files")
    print(full_match)

    # CpGs that are present in both files, but the results do not fully match
    sense_df = df.loc[(df['file_y'] != 'none') & ~((df['meth_x'] == df['meth_y']) & (df['unmeth_x'] == df['unmeth_y']))]
    sense = len(sense_df.index)
    print("Amount of CpGs that are present in both files, but the results do not fully match")
    print(sense)

    # all agreeing CpGs
    agree_tot_df = sense_df.loc[( ( sense_df['proc_x'] > 50 ) & ( sense_df['proc_y'] > 50 ) ) | ( ( sense_df['proc_x'] <= 50 ) & ( sense_df['proc_y'] <= 50 ) )]
    agree_tot = len(agree_tot_df.index)
    print("Total amount of agreeing CpGs")
    print(agree_tot)

    # all disagreeing CpGs
    disagree_tot_df = sense_df.loc[~( ( ( sense_df['proc_x'] > 50 ) & ( sense_df['proc_y'] > 50 ) ) | ( ( sense_df['proc_x'] <= 50 ) & ( sense_df['proc_y'] <= 50 ) ) )]
    disagree_tot = len(disagree_tot_df.index)
    print("Total amount of disagreeing CpGs")
    print(disagree_tot)

    # close agreeing CpGs
    close_agree_df = agree_tot_df.loc[ (abs(agree_tot_df['meth_x'] - agree_tot_df['meth_y']) + abs(agree_tot_df['unmeth_x'] - agree_tot_df['unmeth_y'])) <= 25 ]   
    close_agree = len(close_agree_df.index)
    print("Amount of CpGs that agree and are close")
    print(close_agree)

    # other agreeing CpGs
    agree_df = agree_tot_df.loc[ (abs(agree_tot_df['meth_x'] - agree_tot_df['meth_y']) + abs(agree_tot_df['unmeth_x'] - agree_tot_df['unmeth_y'])) > 25 ]   
    agree = len(agree_df.index)
    print("Amount of CpGs that agree and are not close")
    print(agree)

    # close disagreeing CpGs
    close_disagree_df = disagree_tot_df.loc[ (abs(disagree_tot_df['meth_x'] - disagree_tot_df['meth_y']) + abs(disagree_tot_df['unmeth_x'] - disagree_tot_df['unmeth_y'])) <= 25 ]  
    close_disagree = len(close_disagree_df.index)
    print("Amount of CpGs that disagree and are close")
    print(close_disagree)

    # other disagreeing CpGs
    disagree_df = disagree_tot_df.loc[ (abs(disagree_tot_df['meth_x'] - disagree_tot_df['meth_y']) + abs(disagree_tot_df['unmeth_x'] - disagree_tot_df['unmeth_y'])) > 25 ]
    disagree = len(disagree_df.index)
    print("Amount of CpGs that disagree and are not close")
    print(disagree)

    # Pie chart, where the slices will be ordered and plotted counter-clockwise:
    labels = 'Present in one file', 'Full match', 'Agree and close', 'Agree', 'Disagree and close', 'Disagree'
    sizes = [either, full_match, close_agree, agree, close_disagree, disagree] 
    explode = (0.2, 0, 0, 0, 0, 0)  

    fig1, ax1 = plt.subplots()
    ax1.pie(sizes, labels=labels, explode=explode, autopct='%1.1f%%',
            shadow=True, startangle=90)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

    # constructing the filename for the pie chart .png
    fn = ''.join([res_file[:-4], '_pie_chart.png'])
    # saving the pie chart .png
    plt.savefig(fn)

if __name__=="__main__":
    # getting the arguments for main function from the command line
    res_file = sys.argv[1]
    main(res_file)
