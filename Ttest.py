#!/usr/bin/python

import cPickle as pickle
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
import random
from matplotlib import pyplot

import numpy as np
import pandas as pd
import scipy.stats

import random

def tTest(data1, data2, compName):
    twosample_results = scipy.stats.ttest_ind(data1, data2, equal_var = False)

    mean1 = sum(data1) / float(len(data1))
    mean2 = sum(data2) / float(len(data2))

    mean_dist = mean1 - mean2

    matrix_twosample = ['%r' %compName, twosample_results[0], twosample_results[1], mean_dist]

    return matrix_twosample

# Load in the data - slope_dict
slope_dict = {
"1_M" : random.sample(range(0, 100), 20),
"2_S" : random.sample(range(0, 100), 19),
"3_K" : random.sample(range(0, 100), 21),
"4_S" : random.sample(range(0, 100), 10),
"5_L" : random.sample(range(0, 100), 19),
"6_T" : random.sample(range(0, 100), 14),
"7_H" : random.sample(range(0, 100), 20),
"8_S" : random.sample(range(0, 100), 13),
"9_M" : random.sample(range(0, 100), 4),
"WT" : random.sample(range(0, 100), 90)
}

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
df.to_csv("dummy_Ttest.csv")
