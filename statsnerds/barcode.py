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
