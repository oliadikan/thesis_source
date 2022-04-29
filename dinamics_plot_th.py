import csv
from matplotlib import markers
import matplotlib.pyplot as plt

km = [21, 23, 25, 27, 29, 31]

st = []

for k in km:
    filename = str(k) + "histogram.csv"
    with open(filename, mode='r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=';')
        for row in csv_reader:
            if row[0] == 'count of unique/genome length':
                st.append(float(row[1]))

plt.style.use('seaborn-whitegrid')

fig = plt.figure()
plt.plot(km, st, linestyle='-.', marker='H')
plt.xlim(20, 32)
plt.ylim(0.7, 0.85)
plt.title("Initial dinamics of ucgl stat")
plt.xlabel("k")
plt.ylabel("ucgl")
plt.savefig('init_ucgl_din.png')
