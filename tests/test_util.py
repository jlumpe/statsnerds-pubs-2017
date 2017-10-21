"""Test statsnerds.util module."""

import pytest
import numpy as np

from statsnerds import util


@pytest.mark.parametrize('partial', [False, True])
def test_zip_by_index(partial):
	"""Test the zip_by_index() function."""

	random = np.random.RandomState(0)

	l = 100
	p = .5

	# A bunch of random tests
	for _ in range(1000):

		# Two random enumerated sequences
		seq1 = [(i, (0, i)) for i in range(l) if random.rand() < p]
		seq2 = [(i, (1, i)) for i in range(l) if random.rand() < p]

		seq1_indices = {i for i, a in seq1}
		seq2_indices = {i for i, a in seq2}

		# Check we have the union or intersection of indices
		if partial:
			expected_indices = sorted(seq1_indices | seq2_indices)
		else:
			expected_indices = sorted(seq1_indices & seq2_indices)

		# Go through each
		zipped = util.zip_by_index(seq1, seq2, partial)

		for expected_i, (i, val1, val2) in zip(expected_indices, zipped):

			# Check index is expected value
			assert i == expected_i

			# Check first value
			if val1 is not None:
				a, b = val1
				assert a == 0
				assert b == i
				assert i in seq1_indices
			else:
				assert partial
				assert i not in seq1_indices

			# Check second value
			if val2 is not None:
				a, b = val2
				assert a == 1
				assert b == i
				assert i in seq2_indices
			else:
				assert partial
				assert i not in seq2_indices

		# Check we went through all expected indices
		assert i == expected_indices[-1]

		# Check the zipped iterator is done
		with pytest.raises(StopIteration):
			next(zipped)
