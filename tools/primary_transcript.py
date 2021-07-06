import os
import re
import sys
from collections import Counter, defaultdict

# Use the 'all' version rather than ab initio

def CheckFile(fn):
    """
    Checks for:
    - Duplicated accession lines
    """
    accs = set()
    with open(fn, 'r') as infile:
        for l in infile:
            if l.startswith(">"):
                a = l.rstrip()[1:]
                if a in accs:
                    print("\nERROR: duplicated sequence accession:\n%s" % a)
                    print("\nPlease correct this and then rerun the script.\n")
                    return False
                accs.add(a)
    return True

def ScanTags(fn):
    """
    For ensembl genomes, look for tag:id and count repeated ids
    :param fn:
    :return:
    """
    tags = set()
    tokens = []
    with open(fn, 'r') as infile:
        for line in infile:
            if not line.startswith(">"): continue
            tokens.append([t.split(":", 1) for t in line.rstrip().split() if ":" in t])
            tags.update([t[0] for t in tokens[-1]])
    for this_tag in tags:
        print(this_tag)
        # print(tokens[-1])
        c = Counter([idd for acc in tokens for t, idd in acc if t == this_tag])
        print(c.most_common(5))
        print("")

def ScanTags_NCBI(fn):
    genes = []
    with open(fn, 'r') as infile:
        for line in infile:
            if not line.startswith(">"): continue
            genes.append(line[1:].split(".", 1)[0])
    print("%d sequences, %d genes" % (len(genes), len(set(genes))))

def ScanTags_with_fn(fn, gene_name_fn):
    genes = []
    with open(fn, 'r') as infile:
        for line in infile:
            if not line.startswith(">"): continue
            genes.append(gene_name_fn(line))
    print("%d sequences, %d genes" % (len(genes), len(set(genes))))
    # print(genes[0])
    # print(sorted(genes)[:10])

def GetGeneName_Ensembl(acc_line):
    tokens = [(t.split("=") if "=" in t else t.split(":"))[1] for t in acc_line.rstrip().split() if ("gene:" in t or "gene=" in t or "locus:" in t or "locus=" in t)]
    if len(tokens) != 1: return None
    return tokens[0]

def IsNCBI(fn):
    with open(fn, 'r') as infile:
        for l in infile:
            if l.startswith(">"):
                l = l.rstrip()
                if l.startswith(">NP_") and l.endswith("]"): return True
                elif l.startswith(">XP_") and l.endswith("]"): return True
                elif l.startswith(">YP_") and l.endswith("]"): return True
                elif l.startswith(">WP_") and l.endswith("]"): return True
                return False
    return False

def GetGeneName_NCBI(acc_line):
    acc_line = acc_line[1:]
    original = acc_line
    # look for "isoform X[:d]+" or "isoform [:d]+"
    acc_line = re.sub("isoform [0-9, A-Z]+ ", "", acc_line)
    acc_line = re.sub("isoform X[0-9, A-Z]+ ", "", acc_line)
    # This last step is nasty. These are the same gene:
    # >XP_024356342.1 pyruvate decarboxylase 2-like isoform X1 [Physcomitrella patens]
    # >XP_024356343.1 pyruvate decarboxylase 2-like isoform X1 [Physcomitrella patens]
    # as the name is the same and they same 'isoform ...'
    # Whereas these are not the same gene, even though the names are identical
    # because they don't say isoform:
    # >XP_024390255.1 40S ribosomal protein S12-like [Physcomitrella patens]
    # >XP_024399722.1 40S ribosomal protein S12-like [Physcomitrella patens]
    #
    # To deal with that, we remove the ID (e.g. XP_024356342.1) if it says 'isoform'
    # so that the lines are identical, but not when it doesn't say 'isoform'
    # so that the lines are different. If I were writting the script from scratch
    # for NCBI files I'd do it a different way, but this is a way to handle it so 
    # that it works with the existing logic in the file.
    if original != acc_line:
        acc_line = acc_line.split(None, 1)[-1]
    return acc_line

def CreatePrimaryTranscriptsFile(fn, dout, gene_name_fn, q_use_original_accession_line):
    # Get genes and lengths
    max_gene_lens = defaultdict(int)
    with open(fn, 'r') as infile:
        lines = [l.rstrip() for l in infile]
    N = len(lines) - 1
    nAcc = 0
    nGeneUnidentified = 0
    acc_to_use = defaultdict(str)
    iLine = -1
    while iLine < N:
        iLine += 1
        line = lines[iLine]
        if not line.startswith(">"): continue
        nAcc += 1
        iLineAcc = iLine
        gene = gene_name_fn(line)
        if gene == None:
            nGeneUnidentified += 1
            continue
        # get length
        l = 0
        while iLine < N:
            iLine += 1
            line = lines[iLine]
            if line.startswith(">"):
                iLine -= 1
                break
            l += len(line.rstrip())
        if l > max_gene_lens[gene]:
            max_gene_lens[gene] = l
            acc_to_use[gene] = iLineAcc
    print("Found %d accessions, %d genes, %d unidentified transcripts" % (nAcc, len(max_gene_lens), nGeneUnidentified))
    # print(gene)
    # print(sorted(max_gene_lens.keys())[:10])
    # print(len(set(max_gene_lens.keys())))

    # Get longest version for each gene
    # Parse file second time and only write out sequences that are longest variant
    nGenesWriten = 0
    outfn = dout + os.path.basename(fn)
    with open(outfn, 'w') as outfile:
        iLine = -1
        while iLine < N:
            iLine += 1
            line = lines[iLine]
            if not line.startswith(">"): continue
            gene = gene_name_fn(line)
            # transcripts not identifying the gene should be written
            if gene != None and iLine != acc_to_use[gene]: continue
            if q_use_original_accession_line or gene == None:
                acc_line_out = line + "\n"
            else:
                 acc_line_out = ">%s\n" % gene
            nGenesWriten += 1
            outfile.write(acc_line_out)
            while iLine < N:
                iLine += 1
                line = lines[iLine]
                if line.startswith(">"):
                    iLine -= 1
                    break
                outfile.write(line + "\n")
    print("Wrote %d genes" % nGenesWriten)
    if nGenesWriten != len(max_gene_lens) + nGeneUnidentified:
        print("ERROR")
        raise Exception
    print(outfn)


def last_dot(text):
    return text[1:].rstrip().rsplit(".", 1)[0]


def space(text):
    return text[1:].rstrip().split(None, 1)[0]


function_dict = {"last_dot":last_dot, "space":space}

def main(args=None):
    print("")
    if args is None:
        args = sys.argv[1:]
    fn = args[0]

    if not CheckFile(fn):
        return

    dout = os.path.dirname(os.path.abspath(fn)) + "/primary_transcripts/"
    if not os.path.exists(dout):
        os.mkdir(dout)

    if len(sys.argv) == 3:
        gene_name_function_name = function_dict[sys.argv[2]]
        ScanTags_with_fn(fn, gene_name_function_name)
        CreatePrimaryTranscriptsFile(fn, dout, gene_name_function_name)
    else:
        # ScanTags(fn)
        # ScanTags_NCBI(fn)
        # ScanTags_second_dot(fn)
        if IsNCBI(fn):
            print("Identified as NCBI file")
            gene_name_function_name = GetGeneName_NCBI
            q_use_original_accession_line = True
        else:
            # This is the default, unless we can determine that it is an NCBI
            # file in which case it needs special processing
            gene_name_function_name = GetGeneName_Ensembl
            q_use_original_accession_line = False
            print('Looking for "gene=" of "gene:" to identify isoforms of same gene')
        CreatePrimaryTranscriptsFile(fn, dout, gene_name_function_name, q_use_original_accession_line)

if __name__ == "__main__":
    args = sys.argv[1:]
    main(args)
