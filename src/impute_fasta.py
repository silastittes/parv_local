import textwrap
import gzip
import argparse


parser = argparse.ArgumentParser(
    prog = "Imputed fasta",
    description = "Starting with a gbed file that contains chromosome names and lengths, and a bed file that has a extra column of single nucletodies, generates of fasta with Ns and the nucleotides and the desired positions."
)

parser.add_argument('-g', '--gbed', type=str, required = True,  
            help='An file listing chromosome names and lengths.')

parser.add_argument('-b', '--bed', type=str, required = True,
            help='A bed file with a column for nucleotides to be inserted.')

args = parser.parse_args()

def openfile(filename):
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

#construct a fasta file, based on chromosome names, lengths, and a table of new nucleotides. 
seq_dict = {}
with openfile(args.gbed) as g:
    for line in g:
        chrom, chrom_len = line.strip().split()
        seq_dict[chrom] = ['N']*int(chrom_len)

with openfile(args.bed) as f:
    for line in f:
        chrom, start, end, nuc = line.strip().split()
        pos = int(start)
        seq_dict[chrom][pos] = nuc

for seq_id, seq in seq_dict.items():
    seq_fmt = '\n'.join(textwrap.wrap(''.join(seq), 150))
    print(f">{seq_id}\n{seq_fmt}")
