# thesis_source
 Source files created during the Master's Thesis
 
 - As an initial raw data file (reference human genome) is used **GRCh38.p14** (from [NCBI](https://www.ncbi.nlm.nih.gov/grc/human/data)).
 
 - For obtaining the reverse complements is used the tool [seqkit](https://github.com/shenwei356/seqkit), command is given in the Appendix section.

 - For counting k-mers is used [KMC3.2.1](https://github.com/refresh-bio/KMC), commands are given in the Appendix section. 

---

1. ### inputs_generation

Here are the scripts for the generation of files with chromosome-wise CpG positions and lists of "marker" unique 29-mers, which are used in the main working script of the method. The scripts are numerated in order of their execution, additional usage of tools for obtaining the reverse complements and k-mers counting is required (in our case those are **seqkit** and **KMC3.2.1**), detailed execution example provided in the Appendix section.


---

2. ### method_snakemake

Main working directory of the method. Contains the main working python script of the method, bash scripts for modes, and the **Snakefile**, responsible for the workflow settings. Execution example provided in the Appendix section.

---

3. ### bedGraphs

Contains script for generating .bedGraph file from the resulting output file of the method, as well as the scripts for the comparison of the results of our tool and another tools. 

---

4. ### unique_per_cpg_statistics

Here is the script for the generation of the statistics regarding the amount of CpGs having a certain number of unique 29-mers in their +-50-window (chromosome-wise). 
