"""Tests for statsnerds.io module."""

import random
from pathlib import Path
import io

import pytest
from Bio.Seq import Seq
from Bio import SeqIO

from statsnerds import io as snio


def test_parse_header():
	"""Test parse_illumina_header function."""

	# Example from Illumina help page
	header = '@SIM:1:FCX:1:15:6329:1045 1:N:0:2'
	expected = ('SIM', 1, 'FCX', 1, 15, 6329, 1045, 1, False, 0, 2)

	# Basic parsing
	parsed = snio.parse_illumina_header(header)
	assert isinstance(parsed, snio.IlluminaHeader)
	assert parsed == expected

	# Without @
	assert snio.parse_illumina_header(header[1:]) == expected

	# Invalid format - doesn't match regex
	with pytest.raises(ValueError):
		snio.parse_illumina_header('@this:is:not:a:proper:header')

	# Invalid format - not integers
	with pytest.raises(ValueError):
		snio.parse_illumina_header('@SIM:1:FCX:1:15:notanint:1045 1:N:0:2')

	# Invalid format - not Y/N
	with pytest.raises(ValueError):
		snio.parse_illumina_header('@SIM:1:FCX:1:15:6329:1045 1:Q:0:2')


class TestMultiSequenceReader:
	"""Test the MultiSequenceReader class."""

	N_FILES = 3

	def random_seq(self, minlength=100, maxlength=200):
		"""Create a random sequence.

		:rtype: Bio.Seq.Seq
		"""
		length = random.randint(minlength, maxlength)
		seqstr = ''.join(random.choices('ATGC', k=length))
		return Seq(seqstr)

	def record_eq(self, rec1, rec2):
		"""Check that two sequence records are equal (sequence and ID)."""
		return str(rec1.seq) == str(rec2.seq) and rec1.id == rec2.id

	def tuple_eq(self, tup1, tup2):
		"""Check that two sequence record tuples are equal."""
		return len(tup1) == len(tup2) and \
			all(self.record_eq(t1, t2) for t1, t2 in zip(tup1, tup2))

	@pytest.fixture(scope='class')
	def seq_tuples(self):
		"""List of tuples of sequence records."""
		random.seed(0)

		tuples = []

		for i in range(100):
			records = []

			for j in range(self.N_FILES):
				seq = self.random_seq()
				id_ = 'test-{}-{}'.format(i, j)

				records.append(SeqIO.SeqRecord(seq, id=id_))

			tuples.append(tuple(records))

		return tuples

	@pytest.fixture()
	def seq_files(self, seq_tuples):
		"""Readable file-like objects with sequences written to them."""
		files = []

		for records in zip(*seq_tuples):
			file = io.StringIO()
			SeqIO.write(records, file, 'fasta')
			file.seek(0)
			files.append(file)

		return files

	@pytest.fixture()
	def reader(self, seq_files):
		"""MultiSequenceReader instance on sequence files."""
		return snio.MultiSequenceReader(seq_files, 'fasta')

	def test_iter(self, seq_tuples, reader):
		"""Test basic iteration."""

		seqs_seen = 0

		for tup1, tup2 in zip(seq_tuples, reader):
			assert self.tuple_eq(tup1, tup2)
			seqs_seen += 1

		assert seqs_seen == len(seq_tuples)

		with pytest.raises(StopIteration):
			next(reader)

	@pytest.mark.parametrize('type_', ['str', 'file', 'Path'])
	def test_files_arg_types(self, tmpdir, seq_tuples, type_):
		"""Test different types passed as files argument to constructor."""

		# Write sequences to temporary files
		paths = []

		seq_lists = tuple(map(list, zip(*seq_tuples)))

		for i, seqlist in enumerate(seq_lists):
			path = tmpdir.join('test{}.fasta'.format(i + 1))
			paths.append(path)

			with path.open('wt') as fobj:
				SeqIO.write(seqlist, fobj, 'fasta')

		# Create reader
		if type_ == 'str':
			files_arg = [path.strpath for path in paths]

		elif type_ == 'file':
			files_arg = [path.open() for path in paths]

		elif type_ == 'Path':
			files_arg = [Path(path.strpath) for path in paths]

		reader = snio.MultiSequenceReader(files_arg, 'fasta')

		# Check contents equal
		ntups = 0

		for tup1, tup2 in zip(seq_tuples, reader):
			assert self.tuple_eq(tup1, tup2)
			ntups += 1

		assert ntups == len(seq_tuples)

	def test_close(self, reader):
		"""Test the close() method."""

		for file in reader.files:
			assert not file.closed

		reader.close()

		for file in reader.files:
			assert file.closed

	def test_context(self, reader):
		"""Test using reader as context manager."""

		with reader as rval:
			assert rval is reader

			for file in reader.files:
				assert not file.closed

		for file in reader.files:
			assert file.closed

	def test_skip(self, seq_tuples, reader):
		"""Test the skip() method."""

		n_vals = [5, 2, 13, 1]
		assert sum(n_vals) < len(seq_tuples)

		pos = 0

		# Check skipping before hitting end
		for n in n_vals:

			# Skip and check return value
			assert reader.skip(n) == n

			pos += n

			# Check getting the next tuple works
			expected_next = seq_tuples[pos]
			next_val = next(reader)
			assert self.tuple_eq(next_val, expected_next)

			pos += 1

		# Try skipping past the end
		rval = reader.skip(2 * len(seq_tuples))
		assert rval == len(seq_tuples) - pos

		# Check it's done
		with pytest.raises(StopIteration):
			next(reader)

		# Negative skip count should be an error
		with pytest.raises(ValueError):
			reader.skip(-1)
