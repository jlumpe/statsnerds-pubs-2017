#!/usr/bin/python3

import sys
import gzip
from collections import Counter
from collections import defaultdict
import pickle
#from Bio import SeqIO
#from Bio.Seq import Seq
import random
from matplotlib import pyplot

import numpy as np
import pandas as pd
import scipy.stats

# Sources:
# https://plot.ly/python/t-test/
# https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.stats.f_oneway.html

# Scatter of RMSD and fitness_dict
# volcano plot sanity check -
# untreated fitness matplotlib
# fitness map of each group (associated with SD)

# Load in the data
if len(sys.argv) < 2:
    exit("Usage: python3 fitness_dict.py [barcode:fitness_dict_R1.pkl] [optional: barcode:fitness_dict_R2.pkl]")

# First fitness_dict
fname = sys.argv[1]
rep = fname.split(".")[0]
bc_lib = pickle.load(open("/data/buds/statsnerds/data/barcode_map_v2.pkl", "rb"))
data = pickle.load(open(fname,"rb"), encoding='latin1')

# If there's a second one to compare, load it in
#if sys.argv[2]:
#    fname2 = sys.argv[2]
#    rep2 = fname2.split(".")[0]
#    data2 = pickle.load(open(fname2,"rb"), encoding='latin1')

translate = pickle.load(open("translate.pkl", "rb"))

# # # # # # # # # # # # # FUNCTIONS # # # # # # # # # # # # #

def slope_dict(data, bc_lib):
    """
    Takes a dictionary of format {Barcode: [slope, RMSD]}
    Loads barcode dictionary of format {Barcode: [(codon, position), [quality scores]]}
    (if wild-type, codon and position = None)
    Returns dictionary of format {"pos_aa": [slopes]}
    """

    # Create empty dictionary capable of appending lists to dictionary values
    slope_dict = defaultdict(list)

    # loop through barcodes with reads in the data, matching to codon mutation & translating to aa
    for key in data:
        # If the RMSD is greater than 1.3, we throw it out
        if data[key][1] > 0.75:
            print(data[key])
            continue

        if key in bc_lib:

            slope = float(data[key][0])

            code = bc_lib[key][1]
            pos = bc_lib[key][0]

            # If the barcode refers to mutant (codon/pos != None)
            if code:
                a = [ str(pos + 1), translate[code.replace("T", "U")] ]
                mut = "_".join(a)
            else:
                mut = "WT"

            slope_dict[mut].append( slope )

    return slope_dict

def tTest(data1, data2, compName):
    twosample_results = scipy.stats.ttest_ind(data1, data2, equal_var = False)

    mean1 = sum(data1) / float(len(data1))
    mean2 = sum(data2) / float(len(data2))

    mean_dist = mean2 - mean1

    matrix_twosample = ['%r' %compName, twosample_results[0], twosample_results[1], mean_dist]

    return matrix_twosample

# # # # # # # # # # # # # COMMANDS # # # # # # # # # # # # #
#print(bc_lib)
#print(data)
slope_dict = slope_dict(data, bc_lib)
pickle.dump( slope_dict, open( "mutDict_%s.p" %rep, "wb" ) )

# If we're comparing two replicates, look at R2 & compare visually
#if sys.argv[2]:
#    slope_dict2 = slope_dict(data, bc_lib)
#    pickle.dump( slope_dict, open( "%s_dict.p" %rep2, "wb" ) )

#    plt.axis([0, 10, 0, 1])

#    for key in set(slope_dict, slope_dict2):
#        R1 = slope_dict[key]
#        R2 = slope_dict2[key]
#        plt.scatter(R1,R2)

#    savefig("FitnessCorrelation_%s_vs_%s.png" %(rep, rep2))

# Loop through each key in slope_dict (unless it's WT), save Ttest value to matrix
t_matrix = [['Mutation', 'Test Statistic', 'p-value', 'Mean Distance']]
data1 = slope_dict["WT"]

for key in slope_dict:

    if key == "WT":
        continue

    data2 = slope_dict[key]
    twosample_table = tTest(data1, data2, key)

    t_matrix.append(twosample_table)

# save out pandas DataFrame
headers = t_matrix.pop(0)
df = pd.DataFrame(t_matrix, columns=headers)
df.set_index('Mutation')
df.to_csv("%s_Ttest.csv" %rep)
