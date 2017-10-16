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


def to_file_obj(file, mode='rt'):
	"""Converts a file argument from several formats to an open file object.

	:param file: File path string, path-like object, or file-like object.
	:returns: File-like object.
	"""
	if hasattr(file, '__fspath__'):
		# Implements path-like protocol - new in 3.6
		# Convert to string
		file = str(file.__fspath__())

	if isinstance(file, str):
		# File path as string
		return open(file, mode)

	else:
		# Assume a file-like object
		return file


class MultiSequenceReader:
	"""Reads sequences one at a time from several files in parallel.

	Acts as an iterator over tuples of sequences as :class:`Bio.SeqIO.SeqRecord`
	objets.

	This class is not thread-safe.

	Example usage:

	>>> files = ['reads1.fastq', 'reads2.fastq', 'reads3.fastq']
	>>> MultiSequenceReader(files, 'fastq') as reader:
	>>>     for read1, read2, read3 in reader:
	>>>         # Do stuff with reads...

	TODO: enable easy subsampling, tracking of indices, seeking, etc.

	:param files: Sequence of readable files, as file path strings or file-like
		objects.
	:param str format_: Sequence file format as understood by
		:func:`Bio.SeqIO.parse`.
	"""

	def __init__(self, files, format_):

		if len(files) == 0:
			raise ValueError('files cannot be empty')

		self.format_ = format_
		self.nfiles = len(files)

		# Convert files argument to file objects
		self.files = tuple(map(to_file_obj, files))

		# Create BioPython parser iterators for files.
		self._parsers = tuple(
			SeqIO.parse(fobj, self.format_)
			for fobj in self.files
		)

	def __iter__(self):
		# It's an iterator already, so return self
		return self

	def __enter__(self):
		return self

	def __exit__(self, *args):
		"""Closes all files."""
		self.close()

	def close(self):
		"""Closes all filles the instance has open for reading."""
		for file in self.files:
			file.close()

	def __next__(self):
		"""Get the next set tuple of records.

		:returns: Tuple containing sequence record for each file.
		:rtype: tuple[Bio.SeqIO.SeqRecord]
		"""
		# Check we're not already at the end

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
				warn('Not all sequence files contain the same number of records')

			raise StopIteration()

		else:
			# Got a record from each file
			return tuple(records)

	def _skip_next(self):
		"""Skip over the next set of sequences.

		:returns: True if skipped, False if already at end.
		"""
		try:
			# TODO - can implement this faster for FASTQ by reading over
			# next four lines of each file without parsing
			next(self)
			return True

		except StopIteration:
			return False

	def skip(self, n=1):
		"""Skip over one or more records in each of the files.

		For certain file formats, this is quicker than parsing a record and then
		discarding the result.

		:param int n: Number of sequences to skip over in each file.
		:returns: Number of sequences actually skipped. Will be less than n if
			there are fewer than n records remaining in the files.
		:rtype: int
		"""
		if n < 0:
			raise ValueError('n must be positive')

		nskipped = 0

		for i in range(n):
			if not self._skip_next():
				break

			nskipped += 1

		return nskipped
