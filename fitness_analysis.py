from __future__ import division
import cPickle as pic
from Bio import Seq
from Bio import SeqIO
import numpy as np
from scipy.stats import linregress

AA_TO_NUM = {'A': 1, 'C': 2, 'E': 3, 'D':4, 'G': 5, 'F': 6, 'I': 7, 'H': 8, 'K': 9, 'M': 10, 'L': 11, 'N': 12, 
	      'Q': 13, 'P': 14, 'S': 15, 'R': 16, 'T': 17, 'W': 18, 'V': 19, 'Y': 20}

def parseFASTQ2(f, revcomp):

	reads = []
	qscores = {}
	
	for record in SeqIO.parse(f, "fastq"):
		if revcomp:
			seq = reverse_complement(record.seq)
		else:
			seq = str(record.seq)
		score = record.letter_annotations['phred_quality']

		reads.append(seq[:18])
		qscores[record.id] = score

	return reads, qscores

def calc_freqs(reads, scale):
	"""
	Calculates frequncies to then be put into the fitness score calculations
	Frequenices are normalized by the total number of barcodes in pools
	Scaled by SCALE (which is related to the OD)
	Log scaled with log_2(1+x) 
	"""

	# iterate through reads and keep count of the total number of times
	# you see a particular barcode
	freqs = {}
	for r in reads:
		if freqs.has_key(r):
			freqs[r] += 1
		else:
			freqs[r] = 1	

	# normalize by the size of the pool
	for f in freqs:
		freqs[f] = np.log2( (freqs[f]+1) / len(reads) * (scale))

	return freqs

def calc_fitness(freqs1, freqs2, freqs3):
	
	fitness = {}
	
	for bc in freqs1:
			
		# For now, if a barcode doesn't exist at timepoint 2, skip the barcode
		if bc not in freqs2:
			continue
	
		f1 = freqs1[bc]
		f2 = freqs2[bc]

		# if barcode doesn't exist at time point 3, we can assume it dropped out of the population
		if bc not in freqs3:
			f3 = 0
		else:
			f3 = freqs3[bc]

		f = np.array([f1, f2, f3])
		x = np.array([0,12,24])
			
		# Fit a line to fitness points	
		slope, intercept, r_value, p_value, std_err = linregress(x,f)

		# Calculate RMSD
		fn = slope*x + intercept
		rmsd = np.sqrt(np.mean((f - fn)**2))	
		
		fitness[bc] = [slope, rmsd]

	return fitness
		

		

if __name__ == "__main__":


	#mut_dict = pic.load(open("final_mutation_dictionary.pkl", "rb"))
	
	rep1_od1 = .29
	rep1_od2_1 = 2.18
	rep1_od2_2 = .42
	rep1_od3 = 1.86

	rep2_od1 = .27
	rep2_od2_1 = 1.64
	rep2_od2_2 = .43
	rep2_od3 = 1.8

	print("Read in Yeast Mode reads")
	#ym_1_0, q_1_0 = parseFASTQ2("../BUDS-1-0_S7_L008_R1_001.fastq", False)
	#ym_1_12, q_1_12 = parseFASTQ2("../YM-1-12_S24_L008_R1_001.fastq", False) 
	#ym_1_24, q_1_24 = parseFASTQ2("../YM-1-24_S25_L008_R1_001.fastq", False)
	
	#ym_2_0, q_2_0 = parseFASTQ2("../YM-2-0_S26_L008_R1_001.fastq", False)
	#ym_2_12, q_2_12 = parseFASTQ2("../YM-2-12_S27_L008_R1_001.fastq", False)
	#ym_2_24, q_2_24 = parseFASTQ2("../YM-2-24_S28_L008_R1_001.fastq", False)

	#print("Read in Robert reads")
	r_1_0, rq_1_0 = parseFASTQ2("../RWN-1-0_S1_L008_R1_001.fastq", False)
	r_1_12, rq_1_12 = parseFASTQ2("../RWN-1-12_S2_L008_R1_001.fastq", False) 
	r_1_24, rq_1_24 = parseFASTQ2("../RWN-1-24_S3_L008_R1_001.fastq", False)
	
	r_2_0, rq_2_0 = parseFASTQ2("../RWN-2-0_S4_L008_R1_001.fastq", False)
	r_2_12, rq_2_12 = parseFASTQ2("../RWN-2-12_S5_L008_R1_001.fastq", False)
	r_2_24, rq_2_24 = parseFASTQ2("../RWN-2-24_S6_L008_R1_001.fastq", False)

	print("Calc mut freqs")
	freqs_ym_1_0 = calc_freqs(ym_1_0, 1)
	freqs_ym_1_12 = calc_freqs(ym_1_12, rep1_od2_1/rep1_od1)
	freqs_ym_1_24 = calc_freqs(ym_1_24, (rep1_od2_1 / rep1_od1) * (rep1_od3 / rep1_od2_2))
	freqs_ym_2_0 = calc_freqs(ym_2_0, 1)
	freqs_ym_2_12 = calc_freqs(ym_2_12, rep2_od2_1/rep2_od1)
	freqs_ym_2_24 = calc_freqs(ym_2_24, (rep2_od2_1 / rep2_od1) * (rep2_od3 / rep2_od2_2))
	
	#freqs_r_1_0 = calc_mut_freqs(mut_dict, r_1_0)
	#freqs_r_1_12 = calc_mut_freqs(mut_dict, r_1_12)
	#freqs_r_1_24 = calc_mut_freqs(mut_dict, r_1_24)
	#freqs_r_2_0 = calc_mut_freqs(mut_dict, r_2_0)
	#freqs_r_2_12 = calc_mut_freqs(mut_dict, r_2_12)
	#freqs_r_2_24 = calc_mut_freqs(mut_dict, r_2_24)

	print("Calc fitness")
	fitness_ym_1 = calc_fitness(freqs_ym_1_0, freqs_ym_1_12, freqs_ym_1_24)
	fitness_ym_2 = calc_fitness(freqs_ym_2_0, freqs_ym_2_12, freqs_ym_2_24)

	print("Write fitness arrays")
	pic.dump(fitness_ym_1, open("fitness_ym_1.2.pkl", "wb"))
	pic.dump(fitness_ym_2, open("fitness_ym_2.2.pkl", "wb"))

	print("done!")
	

