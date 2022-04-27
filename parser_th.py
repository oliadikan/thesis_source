# Since we'd like to have an input reference genome file of a bit different format for the further processing, here is the initial parser

with open('/home/olyadikan/Documents/Thesis/GRCh38_latest_genomic.fna', 'r') as ncbi:
    with open('/home/olyadikan/Documents/Thesis/GRCh38.p14.fasta', 'w') as parsed:
        line = ncbi.readline()
        val = False
        while line != '':
            if line[0] == '>':
                if 'NC_0000' in line:
                    val = True
                    parsed.write(line)
                else:
                    val = False
            elif val:
                line = line.upper()
                line = line.replace("N", "")
                if line != "":
                    parsed.write(line)
            line = ncbi.readline()
