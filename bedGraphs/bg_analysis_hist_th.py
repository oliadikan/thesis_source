''' This script produces a series of histogram .png with the statistics described in '4.2 Results comparison'
for the combined .tsv output file from bedgraph_comparison_th.py.

Input - combined .tsv file from bedgraph_comparison_th.py; 
output - series of histogram .png with the statistics described in '4.2 Results comparison'.

Usage: python3 bg_analysis_hist_th.py /path/to/the/combined/file '''

import sys
import numpy as np
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

    # CpGs present in both files
    both_df = df.loc[df['file_y'] != 'none']
    both = len(both_df.index)
    print("Amount of CpGs present in both")
    print(both)

    # values for the first histogram
    both_df["abs"] = abs(both_df['meth_x'] - both_df['meth_y'])  + abs(both_df['unmeth_x'] - both_df['unmeth_y'])
    both_df.loc[both_df['abs'] > 100, 'abs'] = 100

    # values for the second and the fourth histogram
    both_df["diff_meth"] = both_df['meth_y'] - both_df['meth_x']
    both_df.loc[both_df['diff_meth'] > 50, 'diff_meth'] = 50
    both_df.loc[both_df['diff_meth'] < -50, 'diff_meth'] = -50

    # values for the third and the fourth histogram
    both_df["diff_unmeth"] = both_df['unmeth_y'] - both_df['unmeth_x']
    both_df.loc[both_df['diff_unmeth'] > 50, 'diff_unmeth'] = 50
    both_df.loc[both_df['diff_unmeth'] < -50, 'diff_unmeth'] = -50 

    # getting the lists of values for further generation of histograms
    abs_val = both_df["abs"].values.tolist()
    meth_val = both_df["diff_meth"].values.tolist()
    unmeth_val = both_df["diff_unmeth"].values.tolist()

    # the first histogram: absolute values
    bins_abs = [i for i in range(101)]
    plt.hist(abs_val, bins=bins_abs, color="lightseagreen")
    plt.title("Sum of absolute differences")
    # constructing the filename for the first histogram .png
    fn = ''.join([res_file[:-4], '_abs_hist.png'])
    # saving the first histogram .png
    plt.savefig(fn)

    # the second histogram: methylated values
    plt.clf()
    bins_m = [i for i in range(-50, 51)]
    plt.hist(meth_val, bins=bins_m, color="salmon")
    plt.title("Differences of methylated counts (our results - other tool's results)")
    # constructing the filename for the second histogram .png
    fn = ''.join([res_file[:-4], '_meth_hist.png'])
    # saving the second histogram .png
    plt.savefig(fn)

    # the third histogram: unmethylated values
    plt.clf()
    plt.hist(unmeth_val, bins=bins_m, color="mediumslateblue")
    plt.title("Differences of unmethylated counts (our results - other tool's results)")
    # constructing the filename for the third histogram .png
    fn = ''.join([res_file[:-4], '_unmeth_hist.png'])
    # saving the third histogram .png
    plt.savefig(fn)

    # the forth histogram: 2d
    plt.clf()
    bins_2d_meth = [i for i in range(-20, 11)]
    bins_2d_unmeth = [i for i in range(-15, 26)]
    H, yedges, xedges = np.histogram2d(unmeth_val, meth_val, bins=(bins_2d_unmeth, bins_2d_meth))
    plt.pcolormesh(xedges, yedges, H, cmap='rainbow')
    plt.title("2D-Histogram of counts differences")
    plt.xlabel("Methylated differences")
    plt.ylabel("Unmethylated differences")
    plt.colorbar()
    # constructing the filename for the forth histogram .png
    fn = ''.join([res_file[:-4], '_2d_hist.png'])
    # saving the forth histogram .png
    plt.savefig(fn)


if __name__=="__main__":
    # getting the arguments for main function from the command line
    res_file = sys.argv[1]
    main(res_file)
