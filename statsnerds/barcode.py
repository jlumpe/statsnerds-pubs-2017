"""Process barcode and mutation reads into dictionary and perform QC.

.. data::

"""

from warnings import warn

from . import ASYN_WT_SEQ
from .io import MultiSequenceReader


# Default parameter values for process_read_files. You could load alternative
# values from a JSON file, for example.
DEFAULT_PARAMS = {
	'wt_seq': ASYN_WT_SEQ,  # Wild-type DNA sequence as string
	'barcode_len': 18,  # Number of nucelotides in barcode
	'point_mutation_high_qual': 30,  # Min PHRED score for high quality point mutation
	'point_mutation_low_qual': 20,  # Min PHRED score for low quality point mutation
	'max_point_mutations': 10,  # Mutation finding will abort if more than this may point mutations are found
	'max_mutated_codons': 5,  # Abort after seeing this many mutated codons for one read pair
	'codon_mutation_min_qual': 15,  # Won't consider codons where no bases are mutated with a quality >= this value
	'wt_bias': 15,  # Quality bias of each base towards wild-type over mutation
}


class BarcodeError(RuntimeError):
	"""
	Special error class raised when a function can't complete successfully
	because an unexpected condition was found with the reads.
	"""


def process_read_files(files, params=None, progress=False):
	"""Parse and process reads from a triplet of files.

	:param files: Sequence of 3 sequence files, as file path strings or readable
		file-like objects.
	:param dict params: Parameter values. If None will use default.
	:param bool progress: If True display a progress bar with the :mod:`tqdm`,
		module, if installed.
	:param \\**kwargs: Additional parameter values as keyword arguments.

	:returns: TODO
	"""

	# Check 3 files given
	if len(files) != 3:
		raise ValueError('Must give three files')

	# Merge parameters with default
	# If there are duplicate keys, later ones override earlier
	params = {**DEFAULT_PARAMS, **params}

	# Check parameter keys
	for key in params:
		if key not in DEFAULT_PARAMS:
			raise KeyError('Unknown parameter name: {!r}'.format(key))

	# Iterator over read triplets
	reader = MultiSequenceReader(files, 'fastq')
	reads_iter = reader

	if progress:
		# Wrap read iterator in tqdm progress bar if requested
		try:
			from tqdm import tqdm

		except ImportError:
			warn('Unable to import tqdm, cannot display progress')

		else:
			reads_iter = tqdm(reader)

	# TODO - may want to do parsing in a separate thread?
	for i, triplet in enumerate(reads_iter):

		# TODO - distribute this to different processes
		try:
			result = process_triplet(triplet, params)

		except Exception as exc:
			# TODO - log or save this exception somehow
			continue

		# TODO - save result somehow

	# TODO - return collected results


def process_triplet(triplet, params=DEFAULT_PARAMS):
	"""Process a triplet of reads and get mutation, barcode, and stats.

	:param tuple triplet: 3-tuple of :class:`Bio.SeqIO.SeqRecord` objects for
		reads 1, 2, and 3.
	:param dict params: Dictionary of parameter values.
	:returns: TODO
	"""

	read1, read2, read3 = triplet

	# Find mutation
	mutation = find_mutation(read1, read3, params)

	# Find barcode sequence
	barcode = get_barcode(read2, params)

	# TODO - Collect more statistics?

	# TODO - What is the return type?


def find_mutation(read1, read3, params=DEFAULT_PARAMS):
	"""Find the mutation in the WT sequence from reads 1 and 3.

	:param read1: Read 1 of triplet (protein forward).
	:type read1: Bio.SeqIO.SeqRecord
	:param read3: Read 3 of triplet (protein reverse).
	:type read3: Bio.SeqIO.SeqRecord
	:param dict params: Dictionary of parameter values.

	:returns: 2-tuple of ``(residue_index, new_codon)``. If no mutation was
		found and the sequence is believed to be wild type, both tuple elements
		will be None.
	:rtype: tuple[int, str]
	"""
	point_mutations = find_point_mutations(read1, read3, params)

	max_point_muts = params['max_point_mutations']
	if max_point_muts is not None and len(point_mutations) > max_point_muts:
		raise BarcodeError('Maximum number of point mutations exceeded')

	return codon_mutation_from_point(point_mutations, params)


def find_point_mutations(read1, read3, params=DEFAULT_PARAMS):
	"""Find point mutations in reads 1 and 3 relative to wild-type.

	:param read1: Forward read (1st of triplet).
	:type read1: Bio.SeqIO.SeqRecord
	:param read3: Reverse read (3rd of triplet).
	:type read3: Bio.SeqIO.SeqRecord
	:param dict params: Parameter dictionary.

	:returns: List of ``(index, to_nucleotide, quality)`` tuples for each
		point mutation.
	:rtype: list[tuple[int, str, int]]
	"""
	wt = params['wt_seq']

	# Get sequence and quality from each read (RC of read3)
	r1 = read1.seq
	q1 = read1.letter_annotations['phred_quality']

	read3rc = read3.reverse_complement()
	r3 = read3rc.seq
	q3 = read3rc.letter_annotations['phred_quality']

	# Figure out where the overlap occurs
	overlap_start = len(wt) - len(read3)
	overlap_end = len(read1)
	overlap_len = overlap_end - overlap_start
	assert overlap_len > 0

	sub_ls = []

	# Check read 1 before overlap
	for i in range(0, len(r1) - overlap_len):
		if r1[i] != wt[i]:
			sub_ls.append((i, r1[i], q1[i]))

	# Check overlap region
	ovlp_r1, ovlp_r3 = r1[-overlap_len:], r3[:overlap_len]
	ovlp_q1, ovlp_q3 = q1[-overlap_len:], q3[:overlap_len]
	ovlp_wt = wt[overlap_start:overlap_end]

	for i in range(0, len(ovlp_r1)):
		if ovlp_r1[i] == ovlp_r3[i]:
			# Reads agree
			if ovlp_r1[i] != wt[i + overlap_start]:
				sub_ls.append((i + overlap_start, ovlp_r1[i], max(ovlp_q1[i], ovlp_q3[i])))

		else:
			# Reads disagree
			if ovlp_q1[i] > ovlp_q3[i] and ovlp_r1[i] != ovlp_wt[i]:
				sub_ls.append((i + overlap_start, ovlp_r1[i], ovlp_q1[i]))
			elif ovlp_q3[i] > ovlp_q1[i] and ovlp_r3[i] != ovlp_wt[i]:
				sub_ls.append((i + overlap_start, ovlp_r3[i], ovlp_q3[i]))

	# Check read 3 after overlap
	for i in range(overlap_len, len(r3)):
		if r3[i] != wt[i + overlap_start]:
			sub_ls.append((i + overlap_start, r3[i], q3[i]))

	return sub_ls


def codon_mutation_from_point(point_mutations, params=DEFAULT_PARAMS):
	"""Given point mutations found in the reads, determine the codon mutation.

	:param list point_mutations: List of point mutations in format returned by
		:func:`.find_point_mutations`.
	:param dict params: Parameter dictionary

	:returns: 2-tuple of ``(residue_index, new_codon)``. If no mutation was
		found and the sequence is believed to be wild type, both tuple elements
		will be None.
	:rtype: tuple[int, str]

	:raises BacodeError: If the mutation is ambiguous.
	"""
	wt = params['wt_seq']
	q1 = params['point_mutation_high_qual']
	q2 = params['point_mutation_low_qual']

	codon_index = None
	mutations = [None, None, None]

	# First pass to look for high-quality mutations
	for i, b, q in point_mutations:
		if q < q1:
			continue

		# Codon index and base within codon
		c = i // 3
		n = i % 3

		if codon_index is None:
			# First seen, use this as the codon
			codon_index = c
			mutations[n] = b

		elif c == codon_index:
			# In previously seen codon
			mutations[n] = b

		else:
			# There are high-quality mutations in mutliple codons. Give up.
			raise BarcodeError('HQ point mutations in different codons')

	# Second pass to look for lower-quality mutations consistent with this
	if codon_index is None:
		# There were no high-quality mutations, codon position unknown
		# Only accept if all low-quality mutations are in the same codon
		for i, b, q in point_mutations:
			if not (q1 <= q < q2):
				continue

			c = i // 3
			n = i % 3

			if codon_index is None:
				codon_index = c
				mutations[n] = b

			elif c == codon_index:
				mutations[n] = b

			else:
				# There are mutations in different codons, give up.
				raise BarcodeError('LQ point mutations in different codons')

		if codon_index is None:
			# No mutations found, return wild-type
			return (None, None)

	else:
		# Codon position is known from high-quality mutations
		for i, b, q in point_mutations:
			if not (q1 <= q < q2):
				continue

			c = i // 3
			n = i % 3

			# Only keep it if it is in the known codon
			if c == codon_index:
				mutations[n] = b

	# Get codon bases after mutation
	old_codon = wt[codon_index * 3:codon_index * 3 + 3]
	new_codon = ''.join(m or c for m, c in zip(mutations, old_codon))

	return codon_index, new_codon


def get_barcode(read2, params=DEFAULT_PARAMS):
	"""Get the barcode sequence from read 2.

	:param read2: 2nd read of triplet.
	:type read2: Bio.SeqIO.SeqRecord
	:param dict params: Dictionary of parameter values.
	:returns: Barcode sequence. If no valid barcode found, returns None.
	:rtype: str
	"""
	barcode = str(read2[:18])
	return None if 'N' in barcode else barcode


def find_possible_codon_mutations(wt, forward, reverse, params=DEFAULT_PARAMS):
	"""Locate possible mutated codons and get base reads for them.

	:param str wt: Wild-type DNA sequence.
	:param forward: Forward read.
	:type forward: Bio.SeqIO.SeqRecord
	:param reverse: Reverse read.
	:type reverse: Bio.SeqIO.SeqRecord
	:param dict params: Parameter dictionary.

	:returns: Generator yielding tuples of ``(residue_index, codon_reads)``
		where ``codon_reads`` is a length-3 list where each element is a list of
		``(base, quality)`` pairs of base reads for the corresponding position
		in the codon.
	"""

	threshold = params['codon_mutation_min_qual']

	reverse_rc = reverse.reverse_complement()

	forward_qual = forward.letter_annotations['phred_quality']
	reverse_qual = reverse_rc.letter_annotations['phred_quality']

	overlap_start = len(wt) - len(reverse)
	overlap_end = len(forward)

	# For each codon in WT sequence
	for res in range(len(wt) // 3):

		ibegin = res * 3
		iend = ibegin + 3

		wt_codon = wt[ibegin:iend]

		mutations_found = False

		# Check forward read
		if ibegin < overlap_end:
			# Check forward codon (could be partial) against WT
			f_codon = str(forward.seq[ibegin:iend])
			if f_codon != wt_codon[:len(f_codon)]:
				mutations_found = True

		# Check reverse read
		if iend > overlap_start:

			# Begin/end index from start of RC of reverse read
			rbegin = max(ibegin - overlap_start, 0)
			rend = iend - overlap_start

			# Sequence and start index of reverse codon (may be partial)
			r_codon = str(reverse_rc.seq[rbegin:rend])
			r_offset = 3 + rbegin - rend

			# Check reverse codon vs WT
			if r_codon != wt_codon[r_offset:]:
				mutations_found = True

		if mutations_found:
			# Compile reads for each base in codon into list

			codon_reads = [[] for i in range(3)]
			best_mutation_score = 0

			# Add forward base reads
			if ibegin < overlap_end:
				f_codon_qual = forward_qual[ibegin:iend]
				for i, (f_b, q) in enumerate(zip(f_codon, f_codon_qual)):
					codon_reads[i].append((f_b, q))
					if f_b != wt_codon[i] and q > best_mutation_score:
						best_mutation_score = q

			# Add reverse base reads
			if iend > overlap_start:
				r_codon_qual = reverse_qual[rbegin:rend]
				for i, r_b, q in zip(range(r_offset, 3), r_codon, r_codon_qual):
					codon_reads[i].append((r_b, q))
					if r_b != wt_codon[i] and q > best_mutation_score:
						best_mutation_score = q

			if best_mutation_score >= threshold:
				yield res, codon_reads


def score_codon(codon, codon_reads):
	"""Score a codon sequence against reads of its bases.

	:param str codon: Codon sequence.
	:param list codon_reads: Reads for each position of codon. See 2nd element
		of return value of :fund:`.find_possible_codon_mutations`.
	:returns: Score of codon.
	:rtype: int
	"""
	score = 0

	for base, codon_reads in zip(codon, codon_reads):
		for read_base, read_qual in codon_reads:
			if read_base == base:
				score += read_qual

	return score


def get_best_codon_match(wt_codon, allowed, codon_reads):
	"""Get the best codon match for a set of codon reads.

	:param str wt_codon: Sequence of WT codon.
	:param allowed: Set of (non-WT) codons to choose from.
	:type allowed: set[str]
	:param list codon_reads: Reads for each position of codon. See 2nd element
		of return value of :fund:`.find_possible_codon_mutations`.

	:returns: Tuple of ``(best_codon, codon_score, wt_score)``.
	:rtype: tuple[str, int, int]
	"""

	wt_score = score_codon(wt_codon, codon_reads)

	best_mut = None
	best_mut_score = None

	for mut_codon in allowed:
		score = score_codon(mut_codon, codon_reads)
		if best_mut_score is None or score > best_mut_score:
			best_mut = mut_codon
			best_mut_score = score

	return best_mut, best_mut_score, wt_score


def find_mutation_from_allowed(allowed, forward, reverse, params=DEFAULT_PARAMS):
	"""Find a sequence mutation given the allowed mutations per residue/codon.

	:param dict allowed: Mapping from residue/codon index to sets of allowed
		codon sequences.
	:param forward: Forward read.
	:type forward: Bio.SeqIO.SeqRecord
	:param reverse: Reverse read.
	:type reverse: Bio.SeqIO.SeqRecord
	:param dict params: Parameter dictionary.

	:returns: ``(residue_index, codon_seq)`` tuple. If no mutation is detected
		(wild-type) both values will be None.
	:rtype: tuple[int, str]

	:raises BarcodeError: If there seems to be something wrong with the reads
		and the mutation cannot be called with confidence.
	"""

	wt = params['wt_seq']
	max_mutated_codons = params['max_mutated_codons']
	wt_bias = params['wt_bias']

	# Find possible mutated codons
	mutated_codons = []

	for res, reads in find_possible_codon_mutations(wt, forward, reverse):
		if len(mutated_codons) >= max_mutated_codons:
			raise BarcodeError('Maximum number of codon mutations exceeded')
		mutated_codons.append((res, reads))

	best_score = None
	best_codon = None
	best_res = None

	for res, reads in mutated_codons:
		wt_codon = wt[res * 3:res * 3 + 3]
		best_mut, best_mut_score, wt_score = get_best_codon_match(wt_codon, allowed[res], reads)

		score = best_mut_score - wt_score
		if best_score is None or score > best_score:
			best_score = score
			best_res = res
			best_codon = best_mut

	if best_res is not None and best_score >= 3 * wt_bias:
		return (best_res, best_codon)
	else:
		return (None, None)
