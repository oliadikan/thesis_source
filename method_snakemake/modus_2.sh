#!/usr/bin/env bash

# getting the list of samples directories
dirs=(${snakemake_params[drt]})
d0="${dirs[0]}"
d1="${dirs[1]}"
#d2="${dirs[2]}"
#d3="${dirs[3]}"
#d4="${dirs[4]}"
#d5="${dirs[5]}"
#d6="${dirs[6]}"
#d7="${dirs[7]}"
#d8="${dirs[8]}"
#d9="${dirs[9]}"
#d10="${dirs[10]}"

# creation of the modus directory (comment if not needed)
mkdir ${snakemake_params[md2_dir]}
# creation of the directory for log files (comment if not needed)
mkdir ${snakemake_params[l]}
python3 method_rupd_wod_snmk_th.py -d $d0 -c ${snakemake_params[chrom]} -p ${snakemake_params[pr]} -l ${snakemake_params[l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[md2_dir]} -s ${snakemake_input[samples]}
python3 method_rupd_wod_snmk_th.py -d $d1 -c ${snakemake_params[chrom]} -p ${snakemake_params[pr]} -l ${snakemake_params[l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[md2_dir]} -s ${snakemake_input[samples]}
#python3 method_rupd_wod_snmk_th.py -d $d2 -c ${snakemake_params[chrom]} -p ${snakemake_params[pr]} -l ${snakemake_params[l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[md2_dir]} -s ${snakemake_input[samples]}
#python3 method_rupd_wod_snmk_th.py -d $d3 -c ${snakemake_params[chrom]} -p ${snakemake_params[pr]} -l ${snakemake_params[l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[md2_dir]} -s ${snakemake_input[samples]}
#python3 method_rupd_wod_snmk_th.py -d $d4 -c ${snakemake_params[chrom]} -p ${snakemake_params[pr]} -l ${snakemake_params[l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[md2_dir]} -s ${snakemake_input[samples]}
#python3 method_rupd_wod_snmk_th.py -d $d5 -c ${snakemake_params[chrom]} -p ${snakemake_params[pr]} -l ${snakemake_params[l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[md2_dir]} -s ${snakemake_input[samples]}
#python3 method_rupd_wod_snmk_th.py -d $d6 -c ${snakemake_params[chrom]} -p ${snakemake_params[pr]} -l ${snakemake_params[l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[md2_dir]} -s ${snakemake_input[samples]}
#python3 method_rupd_wod_snmk_th.py -d $d7 -c ${snakemake_params[chrom]} -p ${snakemake_params[pr]} -l ${snakemake_params[l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[md2_dir]} -s ${snakemake_input[samples]}
#python3 method_rupd_wod_snmk_th.py -d $d8 -c ${snakemake_params[chrom]} -p ${snakemake_params[pr]} -l ${snakemake_params[l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[md2_dir]} -s ${snakemake_input[samples]}
#python3 method_rupd_wod_snmk_th.py -d $d9 -c ${snakemake_params[chrom]} -p ${snakemake_params[pr]} -l ${snakemake_params[l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[md2_dir]} -s ${snakemake_input[samples]}
#python3 method_rupd_wod_snmk_th.py -d $d10 -c ${snakemake_params[chrom]} -p ${snakemake_params[pr]} -l ${snakemake_params[l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[md2_dir]} -s ${snakemake_input[samples]}
# creation of the signal file
touch ${snakemake_output}
