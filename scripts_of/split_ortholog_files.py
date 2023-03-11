import csv
import os
import sys
import glob


PY2 = sys.version_info <= (3,)
csv_write_mode = 'wb' if PY2 else 'wt'
csv_read_mode = 'rb' if PY2 else 'rt'


def split_ortholog_files(d_ologs):
    if not d_ologs.endswith("/"):
        d_ologs += "/"
    filenames = list(glob.glob(d_ologs + "*.tsv"))
    species = [os.path.splitext(os.path.basename(fn))[0] for fn in filenames]
    for fn, sp0 in zip(filenames, species):
        d_out = d_ologs + "Orthologues_" + sp0 + "/"
        if not os.path.exists(d_out):
            os.mkdir(d_out)
        file_handles = []
        csv_writers = {}
        for sp1 in species:
            if sp0 == sp1:
                continue
            file_handles.append(open(d_out + '%s__v__%s.tsv' % (sp0, sp1), csv_write_mode))
            csv_writers[sp1] = csv.writer(file_handles[-1], delimiter="\t")
            csv_writers[sp1].writerow(("Orthogroup", sp0, sp1))
        with open(fn, csv_read_mode) as infile:
            reader = csv.reader(infile, delimiter="\t")
            next(reader)  # skip header
            for row in reader:
                if len(row) == 4: # OG,species,genes1,genes2
                    csv_writers[row[1]].writerow(row[:1] + row[2:])


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: split_ortholog_files.py orthologs_directory")
        sys.exit()
    d_ologs = sys.argv[1]
    if not os.path.exists(d_ologs):
        print("Directory not found: %s" % d_ologs)
        sys.exit(1)
    split_ortholog_files(d_ologs)
