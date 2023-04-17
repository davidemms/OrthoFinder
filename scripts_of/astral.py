import os
import sys
import glob

from . import tree


def create_input_file(d, output_fn, n_skip=50):
    with open(output_fn, 'w') as outfile:
        for fn in glob.glob(d + "/*txt"):
            iog = int(os.path.basename(fn).split("_")[0][2:])
            if iog < n_skip:
                continue
            t = tree.Tree(fn)
            for n in t:
                n.name = n.name.split("_")[0]
            outfile.write(t.write(format=9) + "\n")


def get_astral_command(astral_input, species_tree, threads):
    return " ".join(["astral-pro", "-i", astral_input, "-o", species_tree, "-t", str(threads)])


if __name__ == "__main__":
    d = sys.argv[1]
    output_fn = "astral_input.nwk"
    create_input_file(d, output_fn)
