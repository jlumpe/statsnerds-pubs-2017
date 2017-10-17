import gzip
from Bio.Seq import Seq
from Bio import SeqIO
import pickle
import sys

r1_fastq, r3_fastq = sys.argv[1], sys.argv[2]
wt = "ATGGATGTATTCATGAAAGGACTTTCAAAGGCCAAGGAGGGAGTTGTGGCTGCTGCTGAGAAAACCAAACAGGGTGTGGCAGAAGCAGCAGGAAAGACAAAAGAGGGTGTTCTCTATGTAGGCTCCAAAACCAAGGAGGGAGTGGTGCATGGTGTGGCAACAGTGGCTGAGAAGACCAAAGAGCAAGTGACAAATGTTGGAGGAGCAGTGGTGACGGGTGTGACAGCAGTAGCCCAGAAGACAGTGGAGGGAGCAGGGAGCATTGCAGCAGCCACTGGCTTTGTCAAAAAGGACCAGTTGGGCAAGAATGAAGAAGGAGCCCCACAGGAAGGAATTCTGGAAGATATGCCTGTGGATCCTGACAATGAGGCTTATGAAATGCCTTCTGAGGAAGGGTATCAAGACTACGAACCTGAAGCC"

# Create iterators for each read
read1 = SeqIO.parse(gzip.open(r1_fastq), "fastq")
read3 = SeqIO.parse(gzip.open(r3_fastq), "fastq")

def index_mutations(r1, q1, r3, q3, head, ls):
    sub_ls = [head]
    for i in range(0, len(r1)-80):
        if r1[i] != wt[i]:
            sub_ls.append((i, r1[i], q1[i]))
    
    ovlp_r1, ovlp_r3 = r1[-80:], r3[:80]
    ovlp_q1, ovlp_q3 = q1[-80:], q3[:80]
    ovlp_wt = wt[170:250]
    for i in range(0, len(ovlp_r1)):
        if ovlp_r1[i] == ovlp_r3[i]:
            if ovlp_r1[i] != wt[i]:
                sub_ls.append((i+170, ovlp_r1[i], max(ovlp_q1[i], ovlp_q3[i])))
        else:
            if ovlp_q1[i] > ovlp_q3[i] and ovlp_r1[i] != ovlp_wt[i]:
                sub_ls.append((i+170, ovlp_r1[i], ovlp_q1[i]))
            elif ovlp_q3[i] > ovlp_q1[i] and ovlp_r3[i] != ovlp_wt[i]:
                sub_ls.append((i+170, ovlp_r3[i], ovlp_q3[i]))

    for i in range(80, len(r3)):
        if r3[i] != wt[i+170]:
            sub_ls.append((i+170, r3[i], q3[i]))
    ls.append(sub_ls)
    return ls
                
positions = []
for record in read1:
    head = str(record.id)
    r1, q1 = str(record.seq), record.letter_annotations["phred_quality"]
    third = read3.next()
    r3,q3 = str(third.seq.reverse_complement()), record.letter_annotations["phred_quality"][::-1]

    positions = index_mutations(r1, q1, r3, q3, head, positions)
print(positions[0])
with open("mut_positions.pkl", "w") as out:
    pickle.dump(positions, out)
