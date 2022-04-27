import csv

km = [21, 23, 25, 27, 29, 31]

for k in km:
    filename = str(k) + "histogram.csv"
    unique = 0
    #total = 0
    with open(filename, mode='r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=';')
        line_count = 0
        for row in csv_reader:
            if line_count != 0:
                #total += int(row[1])
                if line_count == 1:
                    unique += int(row[1])
            line_count += 1    

    ucgl = (unique/2976242597)
    with open(filename, mode='a') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=';')
        csv_writer.writerow(['count of unique/genome length', str(ucgl)])
