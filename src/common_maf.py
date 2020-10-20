import gzip
import argparse

parser = argparse.ArgumentParser(
    
prog = "maf2treemix",    

description="Convert multiple maf files to a single treemix file.")

parser.add_argument('-f', '--mafs', type=str,
            help='Input file with two fields, the population name to appear in the final treemix file and the the name of the maf file for that pop.')

args = parser.parse_args()


#read through a list of mafs to find sites common to all of them.
#once common sites are found, generate single input file for tree mix, converted mafs to counts as it goes
#chromo  position        major   minor   ref     anc     unknownEM       nInd
#chr1    65      A       C       A       A       0.000006        1
#chr1    66      C       A       C       C       0.000006        2

def openfile(filename):
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

all_files = {}

with openfile(args.mafs) as f:
    for line in f:
        ln = line.strip().split()
        pop, file_name = ln[0:2]
        all_files[pop] = file_name

#read through files, keeping all sites of the first files, then only common sites after that
common_sites = []
for pop, file_name in all_files.items():
    with openfile(file_name) as f:
        line1 = f.readline()
        for line in f:
            ln = line.split()
            marker = f"{ln[0]}_{ln[1]}"
            common_sites.append(marker)

#bad idea?
common_sites = list(set(common_sites))
common_dict = dict.fromkeys(common_sites, [])

print(common_sites)

#make treemix input
#filter by sites with sufficient variation?
#pop_dict = dict.fromkeys(all_files.keys(), [])

for pop, file_name in all_files.items():
    with openfile(file_name) as f:
        line1 = f.readline()
        for line in f:
            ln = line.strip().split()
            marker = f"{ln[0]}_{ln[1]}"
            if marker in common_dict:
                chrom, pos, major, minor, ref, anc, freq, nInd = ln
                a1 = int(float(freq) * float(nInd))
                a2 = int((1-float(freq)) * float(nInd))
                print(pop, marker, freq, nInd)
                common_dict[marker].append(f"{a1},{a2}")
                #pop_dict[pop].append(f"{a1},{a2}")

#print(pop_dict)
print(f"{marker}\t{' '.join(list(all_files.keys()))}")
for site, markers in common_dict.items():
    print(f"{site}\t{' '.join(markers)}")

                
