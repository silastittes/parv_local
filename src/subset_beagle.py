import argparse
import gzip
import numpy as np

#from parse_read import parse_read

prog = 'Simple helper script to subset a beagle file based on the Ind[0-9][0-9]* labels in the header',

parser = argparse.ArgumentParser(description="Given a list of Ind IDs in a beagle file, subsets the beagle file for those samples. Silly but useful.")

parser.add_argument("-b", "--beagle_file", type = str, help = "Name of the beagle file to subset", required = True)

parser.add_argument("-i", "--id_file", type = str, help = "Name of the ID file with and Ind[0-9][0-9] on each line.", required = True)

args = parser.parse_args()

def openfile(filename):
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

#add ind ids to dictionary
id_list = []
with openfile(args.id_file) as f:
    for line in f:
        ind_id = line.strip().split()[0]
        id_list.append(ind_id)
id_list = np.asarray(id_list)

#filter beagle file
with openfile(args.beagle_file) as f:
    line1 = np.array(f.readline().strip().split())
    id_match = np.where(np.in1d(line1, id_list))[0]
    beagle_idx = np.concatenate((np.array([0,1,2]), id_match))

    print('\t'.join(line1[beagle_idx]))

    for line in f:
        ln = np.array(line.strip().split())
        print('\t'.join(ln[beagle_idx]))
        
    
