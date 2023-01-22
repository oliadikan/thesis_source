#!/usr/bin/env bash

# creation of the modus directory (comment if not needed)
mkdir ${snakemake_params[md1_dir]}
# creation of the output directory (comment if not needed)
mkdir ${snakemake_params[out]}
# creation of the directory for log files (comment if not needed)
mkdir ${snakemake_params[out_l]}
python3 method_rupd_wod_snmk_th.py -d ${snakemake_params[drt]} -c Y -p 10 -l ${snakemake_params[out_l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[out]} -s ${snakemake_input[samples]}
python3 method_rupd_wod_snmk_th.py -d ${snakemake_params[drt]} -c 22 -p 10 -l ${snakemake_params[out_l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[out]} -s ${snakemake_input[samples]}
python3 method_rupd_wod_snmk_th.py -d ${snakemake_params[drt]} -c 21 -p 10 -l ${snakemake_params[out_l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[out]} -s ${snakemake_input[samples]}
python3 method_rupd_wod_snmk_th.py -d ${snakemake_params[drt]} -c 20 -p 10 -l ${snakemake_params[out_l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[out]} -s ${snakemake_input[samples]}
python3 method_rupd_wod_snmk_th.py -d ${snakemake_params[drt]} -c 19 -p 10 -l ${snakemake_params[out_l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[out]} -s ${snakemake_input[samples]}
python3 method_rupd_wod_snmk_th.py -d ${snakemake_params[drt]} -c 18 -p 10 -l ${snakemake_params[out_l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[out]} -s ${snakemake_input[samples]}
python3 method_rupd_wod_snmk_th.py -d ${snakemake_params[drt]} -c 17 -p 7 -l ${snakemake_params[out_l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[out]} -s ${snakemake_input[samples]}
python3 method_rupd_wod_snmk_th.py -d ${snakemake_params[drt]} -c 16 -p 7 -l ${snakemake_params[out_l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[out]} -s ${snakemake_input[samples]}
python3 method_rupd_wod_snmk_th.py -d ${snakemake_params[drt]} -c 15 -p 7 -l ${snakemake_params[out_l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[out]} -s ${snakemake_input[samples]}
python3 method_rupd_wod_snmk_th.py -d ${snakemake_params[drt]} -c 14 -p 7 -l ${snakemake_params[out_l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[out]} -s ${snakemake_input[samples]}
python3 method_rupd_wod_snmk_th.py -d ${snakemake_params[drt]} -c 13 -p 7 -l ${snakemake_params[out_l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[out]} -s ${snakemake_input[samples]}
python3 method_rupd_wod_snmk_th.py -d ${snakemake_params[drt]} -c X -p 7 -l ${snakemake_params[out_l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[out]} -s ${snakemake_input[samples]}
python3 method_rupd_wod_snmk_th.py -d ${snakemake_params[drt]} -c 12 -p 5 -l ${snakemake_params[out_l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[out]} -s ${snakemake_input[samples]}
python3 method_rupd_wod_snmk_th.py -d ${snakemake_params[drt]} -c 11 -p 5 -l ${snakemake_params[out_l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[out]} -s ${snakemake_input[samples]}
python3 method_rupd_wod_snmk_th.py -d ${snakemake_params[drt]} -c 10 -p 3 -l ${snakemake_params[out_l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[out]} -s ${snakemake_input[samples]}
python3 method_rupd_wod_snmk_th.py -d ${snakemake_params[drt]} -c 9 -p 5 -l ${snakemake_params[out_l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[out]} -s ${snakemake_input[samples]}
python3 method_rupd_wod_snmk_th.py -d ${snakemake_params[drt]} -c 8 -p 3 -l ${snakemake_params[out_l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[out]} -s ${snakemake_input[samples]}
python3 method_rupd_wod_snmk_th.py -d ${snakemake_params[drt]} -c 7 -p 2 -l ${snakemake_params[out_l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[out]} -s ${snakemake_input[samples]}
python3 method_rupd_wod_snmk_th.py -d ${snakemake_params[drt]} -c 6 -p 2 -l ${snakemake_params[out_l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[out]} -s ${snakemake_input[samples]}
python3 method_rupd_wod_snmk_th.py -d ${snakemake_params[drt]} -c 5 -p 2 -l ${snakemake_params[out_l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[out]} -s ${snakemake_input[samples]}
python3 method_rupd_wod_snmk_th.py -d ${snakemake_params[drt]} -c 4 -p 2 -l ${snakemake_params[out_l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[out]} -s ${snakemake_input[samples]}
python3 method_rupd_wod_snmk_th.py -d ${snakemake_params[drt]} -c 3 -p 1 -l ${snakemake_params[out_l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[out]} -s ${snakemake_input[samples]}
python3 method_rupd_wod_snmk_th.py -d ${snakemake_params[drt]} -c 2 -p 1 -l ${snakemake_params[out_l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[out]} -s ${snakemake_input[samples]}
python3 method_rupd_wod_snmk_th.py -d ${snakemake_params[drt]} -c 1 -p 1 -l ${snakemake_params[out_l]} -i ${snakemake_input[lists]} -k ${snakemake_input[cpgs]} -o ${snakemake_params[out]} -s ${snakemake_input[samples]}
# creation of the signal file
touch ${snakemake_output}
