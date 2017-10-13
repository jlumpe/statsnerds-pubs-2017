"""Process barcode and mutation reads into dictionary and perform QC.

.. data::

"""

from warnings import warn

from .io import MultiSequenceReader


# Default parameter values for process_read_files. You could load alternative
# values from a JSON file, for example.
DEFAULT_PARAMS = {
	# 'some_parameter': 42,
}


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

	if mutation is None:
		# TODO - couldn't get a valid mutation, log error somehow
		return None

	# Find barcode sequence
	barcode = get_barcode(read2, params)

	if barcode is None:
		# TODO - couldn't get a valid barcode, log error somehow
		return None

	# TODO - Collect more statistics?

	# TODO - What is the return type?


def find_mutation(read1, read3, params=DEFAULT_PARAMS):
	"""Find the mutation in the WT sequence from reads 1 and 3.

	:param read1: Read 1 of triplet (protein forward).
	:type read1: Bio.SeqIO.SeqRecord
	:param read3: Read 3 of triplet (protein reverse).
	:type read3: Bio.SeqIO.SeqRecord
	:param dict params: Dictionary of parameter values.
	:returns: 2-tuple of ``(residue_index, new_codon)``. If no valid mutation
		was found, both tuple elements will be None.
	"""
	raise NotImplementedError()  # TODO


def get_barcode(read2, params=DEFAULT_PARAMS):
	"""Get the barcode sequence from read 2.

	:param read2: 2nd read of triplet.
	:type read2: Bio.SeqIO.SeqRecord
	:param dict params: Dictionary of parameter values.
	:returns: Barcode sequence. If no valid barcode found, returns None.
	:rtype: str
	"""
	raise NotImplementedError()  # TODO
