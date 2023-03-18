import csv
import os
import sys
import glob
import gzip
import argparse

PY2 = sys.version_info <= (3,)
csv_write_mode = 'wb' if PY2 else 'wt'
csv_read_mode = 'rb' if PY2 else 'rt'


def split_ortholog_files(d_ologs, q_compress=False):
    if not d_ologs.endswith("/"):
        d_ologs += "/"
    filenames = list(glob.glob(d_ologs + "*.tsv") + glob.glob(d_ologs + "*.tsv.gz"))
    species = [
        os.path.splitext(os.path.splitext(os.path.basename(fn))[0])[0] if fn.endswith(".gz") else os.path.splitext(os.path.basename(fn))[0]
        for fn in filenames]
    for fn, sp0 in zip(filenames, species):
        d_out = d_ologs + "Orthologues_" + sp0 + "/"
        if not os.path.exists(d_out):
            os.mkdir(d_out)
        file_handles = []
        csv_writers = {}
        for sp1 in species:
            if sp0 == sp1:
                continue
            if q_compress:
                file_handles.append(gzip.open(d_out + '%s__v__%s.tsv.gz' % (sp0, sp1), csv_write_mode))
            else:
                file_handles.append(open(d_out + '%s__v__%s.tsv' % (sp0, sp1), csv_write_mode))
            csv_writers[sp1] = csv.writer(file_handles[-1], delimiter="\t")
            csv_writers[sp1].writerow(("Orthogroup", sp0, sp1))
        with gzip.open(fn, csv_read_mode) if fn.endswith(".gz") else open(fn, csv_read_mode) as infile:
            reader = csv.reader(infile, delimiter="\t")
            next(reader)  # skip header
            for row in reader:
                if len(row) == 4: # OG,species,genes1,genes2
                    csv_writers[row[1]].writerow(row[:1] + row[2:])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("orthologs_directory", help="Input directory")
    parser.add_argument("-c", "--compress", action="store_true", help="Compress output files")
    args = parser.parse_args()
    if not os.path.exists(args.orthologs_directory):
        print("Directory not found: %s" % args.orthologs_directory)
        sys.exit(1)
    split_ortholog_files(args.orthologs_directory, q_compress=args.compress)
