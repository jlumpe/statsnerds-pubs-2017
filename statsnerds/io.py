"""Tools for reading/writing/parsing sequence data, etc."""

import re
from collections import namedtuple
from warnings import warn

from Bio import SeqIO


# Regular expression for Illumina FASTQ header
ILLUMINA_HEADER_RE = re.compile(r'([^:]+):(\d+):([^:]+):(\d+):(\d+):(\d+):(\d+) (\d+):([YN]):(\d+):(\d+)')


# Named tuple for illumina header data
IlluminaHeader = namedtuple(
	'IlluminaHeader',
	'instrument run flowcell_id lane tile xpos ypos read is_filtered controlnum, samplenum',
	module=__name__,
)


def parse_illumina_header(header):
	"""Parse the header line from a read in an Illumia FASTQ file.

	Source: http://support.illumina.com/content/dam/illumina-support/help/BaseSpaceHelp_v2/Content/Vault/Informatics/Sequencing_Analysis/BS/swSEQ_mBS_FASTQFiles.htm

	:param str header: Header line text, first "@" character optional.

	:rtype: .IlluminaHeader
	"""
	# Strip whitespace and remove leading @
	header = header.strip()
	if header[0] == '@':
		header = header[1:]

	# Parse w/ regex
	match = ILLUMINA_HEADER_RE.match(header)

	if match is None:
		raise ValueError('header not in expected format')

	# Unpack
	inst, run, fcell, lane, tile, x, y, read, filtered, control, sample = match.groups()

	try:
		return IlluminaHeader(
			instrument=inst,
			run=int(run),
			flowcell_id=fcell,
			lane=int(lane),
			tile=int(tile),
			xpos=int(x),
			ypos=int(y),
			read=int(read),
			is_filtered=filtered == 'Y',
			controlnum=int(control),
			samplenum=int(sample),
		)

	except ValueError as exc:
		# Can't parse int
		raise ValueError('Invalid header line: {}'.format(exc)) from None


class MultiSequenceReader:
	"""Reads sequences one at a time from several files in parallel.

	Acts as an iterator over tuples of reads, as :class:`Bio.SeqIO.SeqRecord`
	objets.

	This class is not thread-safe.

	Example usage:

	>>> files = ['reads1.fastq', 'reads2.fastq', 'reads3.fastq']
	>>> MultiSequenceReader(files, 'fastq') as reader:
	>>>     for read1, read2, read3 in reader:
	>>>         # Do stuff with reads...

	TODO: enable easy subsampling, tracking of indices, seeking, etc.

	:param files: Sequence of three sequence files, as file path strings or
		readable file-like objects
	:param str format_: Sequence file format as understood by
		:func:`Bio.SeqIO.parse`.
	"""

	def __init__(self, files, format_):

		if len(files) == 0:
			raise ValueError('files cannot be empty')

		self.format_ = format_
		self.nfiles = len(files)

		# Convert files argument to file objects
		file_objs = []

		for file in files:
			if hasattr(file, '__fspath__'):
				# Implements path-like protocol - new in 3.6
				# Convert to string
				file = str(file.__fspath__())

			if isinstance(file, str):
				# File path as string
				file_objs.append(open(file))

			else:
				# Assume a file-like object
				file_objs.append(file)

		self.files = tuple(file_objs)

		# Biopython parser iterators
		self._parsers = tuple(SeqIO.parse(fobj, format_) for fobj in self.files)

	def __iter__(self):
		# It's an iterator already, so return self
		return self

	def __next__(self):
		# Get the next set tuple of records

		records = []
		finished = False

		# Parse one record from each parser
		for parser in self._parsers:

			try:
				rec = next(parser)

			except StopIteration:
				# Reached end of this file
				finished = True
				continue

			records.append(rec)

		if finished:
			# At least one file ran out

			if records:
				# Some files still had records left
				warn('Not all sequence files contain the same number of records')

			# This is a special exception raised by iterators when there is no
			# next item
			raise StopIteration()

		else:
			# Got a record from each file
			return tuple(records)

	def skip(self, n=1):
		"""Skip over one or more records in each of the files.

		For certain file formats, this is quicker than parsing a record and then
		discarding the result.

		:param int n: Number of reads to skip over in each file.
		:returns: Number of reads actually skipped. Will be less than n if there
			are fewer than n records remaining in the files.
		:rtype: int
		"""
		nskipped = 0

		for i in range(n):
			try:
				# TODO - can implement this faster for FASTQ by reading over
				# next four lines of each file without parsing
				next(self)
			except StopIteration:
				break

			nskipped += 1

		return nskipped

	def tell(self):
		"""Get the index of the next set of sequence records to be read.

		:rtype: int
		"""
		raise NotImplementedError()  # TODO

	def seek(self, index):
		"""Seek all files to the record with the given index.

		:param int index: Index of next set of records to read.
		"""
		raise NotImplementedError()  # TODO

	def reset(self):
		"""Seek back to the beginning of all files."""
		raise NotImplementedError()  # TODO
