import numpy
import cPickle as pic
from Bio import SeqIO

BP_COMP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

def reverse_complement(x):
	comp = ""
	rev = x[::-1]
	for i in rev:
		comp += BP_COMP[i]
	return comp

def parseFASTQ2(f, revcomp):

	reads = {}
	qscores = {}
	
	for record in SeqIO.parse(f, "fastq"):
		if revcomp:
			seq = reverse_complement(record.seq)
		else:
			seq = str(record.seq)
		score = record.letter_annotations['phred_quality']

		reads[record.id] = seq
		qscores[record.id] = score

	return reads, qscores

def createRefDict(r1s, r2s, r3s):
	
	output = []

	for k in r1s.keys():
		output.append([r1s[k], r2s[k], r3s[k]])
	return output

def createBarcodeToMutationDict(reads):

	mutdict ={}
	for read in reads:
		bc = read[0]
		seq = read[1]
		qs = read[2]
		if "N" not in bc and "N" not in seq:
			mutdict[bc] = (seq, qs)
	
	return mutdict

def getCodonMutations(seq1, seq2):
	
	muts = []
	for i in range(0, len(seq1), 3):
		c1 = seq1[i:i+3]
		c2 = seq2[i:i+3]
		if c1 != c2:
			muts.append((c2, i/3))
	return muts

def filterDict(WT, mutdict):

	fin_mutdict = {}
	for k in mutdict.keys():
	
		seq = mutdict[k][0]
		qs = mutdict[k][1]
		muts = getCodonMutations(WT, seq)
		
		
		if len(muts) == 1:
			print "THERE IS ONE THAT IS ONE LENGTH"
			fin_mutdict[k] = (muts[0], qs) 
	return fin_mutdict

def constructSeqs(wt, reads, scores):
	
	wt_length = len(wt)

	good_reads = []
	ambig_reads = []
	for i in range(len(reads)):
		r1 = reads[i][0]
		bc = reads[i][1]
		r3 = reads[i][2]
		q1 = scores[i][0]
		q3 = scores[i][2]

		for k in range(len(r1), wt_length):
			r1 += "x"
			q1.append(-1)
		r3_diff = wt_length - len(r3)
		r3 = ("x")*r3_diff + r3
		q3 = [-1]*r3_diff + q3

		seq = ""
		score = []
		for j in range(len(r1)):
			bp1 = r1[j]
			bp2 = r3[j]
			s1 = q1[j]
			s2 = q3[j]
			if s1 > s2:
				seq += bp1
			elif s2 > s1:
				seq += bp2	
			else:
				if bp1 == wt[j]:
					seq += bp1
				else:
					seq += bp2 	
			score.append((int(s1) + int(s2)) / 2)	
		if ("N" in seq) or ("N" in bc):
			ambig_reads.append((bc, seq, score))
		else:
			good_reads.append((bc, seq, score))
	return good_reads, ambig_reads

if __name__ == "__main__":
	
	WT = "ATGGATGTATTCATGAAAGGACTTTCAAAGGCCAAGGAGGGAGTTGTGGCTGCTGCTGAGAAAACCAAACAGGGTGTGGCAGAAGCAGCAGGAAAGACAAAAGAGGGTGTTCTCTATGTAGGCTCCAAAACCAAGGAGGGAGTGGTGCATGGTGTGGCAACAGTGGCTGAGAAGACCAAAGAGCAAGTGACAAATGTTGGAGGAGCAGTGGTGACGGGTGTGACAGCAGTAGCCCAGAAGACAGTGGAGGGAGCAGGGAGCATTGCAGCAGCCACTGGCTTTGTCAAAAAGGACCAGTTGGGCAAGAATGAAGAAGGAGCCCCACAGGAAGGAATTCTGGAAGATATGCCTGTGGATCCTGACAATGAGGCTTATGAAATGCCTTCTGAGGAAGGGTATCAAGACTACGAACCTGAAGCC"

	print("Read 1")
	r1s, q1s = parseFASTQ2("../Undetermined_S0_L001_R1_001.fastq", False)
	print("Read 2")
	r2s, q2s = parseFASTQ2("../Undetermined_S0_L001_R2_001.fastq", False)
	print("Read 3")
	r3s, q3s = parseFASTQ2("../Undetermined_S0_L001_R3_001.fastq", True)

	#print("Read 1")
	#r1s, q1s = parseFASTQ2("../FASTQ_1/r1_test.fastq", False)
	#print("Read 2")
	#r2s, q2s = parseFASTQ2("../FASTQ_1/r2_test.fastq", False)
	#print("Read 3")
	#r3s, q3s = parseFASTQ2("../FASTQ_1/r3_test.fastq", True)

	read_list = createRefDict(r1s, r2s, r3s)
	print("Created read reference")

	qscore_list = createRefDict(q1s, q2s, q3s)
	print("Created qscore reference")

	greads, areads = constructSeqs(WT, read_list, qscore_list)
	print("Constructed seqs")
	
	#greads = pic.load(open("good_reads.pkl", "rb"))
	#print greads
	mutdict = createBarcodeToMutationDict(greads)
	print("Created mutation dictionary")
	
	fin_mutdict = filterDict(WT, mutdict)
	print("Created final mutation dictionary")

	pic.dump(greads, open("good_reads.pkl", 'w'))
	pic.dump(areads, open("ambig_reads.pkl", "w"))
	pic.dump(mutdict, open("mutation_dictionary.pkl","wb"))
	#print mutdict
	pic.dump(fin_mutdict, open("final_mutation_dictionary.pkl", "wb"))
	#print fin_mutdict
	print "done"
