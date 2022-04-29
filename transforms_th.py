# Here are the transformations C -> T and G -> A

with open('/home/olyadikan/Documents/Thesis/GRCh38.p14.fasta', 'r') as fasta:
     with open('/home/olyadikan/Documents/Thesis/GRCh38.p14_CtoT.fasta', 'w') as ct:
         with open('/home/olyadikan/Documents/Thesis/GRCh38.p14_GtoA.fasta', 'w') as ga:
             line = fasta.readline()
             while line != '':
                 line_ct = line
                 line_ga = line
                 if line[0] != '>':
                     line_ct = line_ct.replace('C', 'T')
                     line_ga = line_ga.replace('G', 'A')
                 ct.write(line_ct)
                 ga.write(line_ga)
                 line = fasta.readline()
