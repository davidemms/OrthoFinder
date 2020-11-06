import os
import sys
import csv

if len(sys.argv) != 2 or sys.argv[1] == "-h" or sys.argv[1] == "-help" or sys.argv[1] == "--help":
    print("Usage: orthogroup_gene_count.py Orthogroups.csv")
    sys.exit()

inFN = sys.argv[1]
outFN = os.path.splitext(inFN)[0] + ".GeneCount.csv"

with open(inFN, 'r') as infile, open(outFN, 'w') as outfile:
    reader = csv.reader(infile, delimiter="\t")
    writer = csv.writer(outfile, delimiter="\t")
    header = next(reader)
    n_col_skip = 3 if header[0] == "HOG" else 1
    writer.writerow(header)
    for line in reader:
        writer.writerow(line[:n_col_skip] + [0 if "" == cell else len(cell.split(", ")) for cell in line[n_col_skip:]])
print("Orthogroup gene count table has been written to %s" % outFN)
