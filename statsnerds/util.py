"""Misc utility code."""


def zip_by_index(iterable1, iterable2, partial=False):

	# Convert to iterators
	iters = [iter(iterable1), iter(iterable2)]

	last_i = [-1, -1]  # Last index of each
	last_vals = [None, None]  # Last value of each
	last_used = [True, True]  # Whether the last value was yielded yet

	# Get partial yield value (index not seen in other)
	def _partial_yield_val(which_iter):
		if which_iter == 0:
			return (last_i[0], last_vals[0], None)
		else:
			return (last_i[1], None, last_vals[1])

	while True:

		# Choose the iterator with the lowest index to go next
		next_iter = 0 if last_i[0] <= last_i[1] else 1
		const_iter = 0 if next_iter == 1 else 1

		# Should have already used last before getting a new one
		assert last_used[next_iter]

		# Try getting the next item
		try:
			i, val = next(iters[next_iter])

		except StopIteration:
			# This one ran out

			if partial and not last_used[const_iter]:
				# Yield the last value from the other one if we haven't yet
				yield _partial_yield_val(const_iter)

			break

		# Check we did move ahead
		assert i > last_i[next_iter]

		# Update values for iterator which was advanced
		last_i[next_iter] = i
		last_vals[next_iter] = val
		last_used[next_iter] = False

		# Check the new index
		if last_i[next_iter] == last_i[const_iter]:
			# Same index, yield both
			assert not last_used[const_iter]
			yield (i, *last_vals)
			last_used = [True, True]

		elif last_i[next_iter] < last_i[const_iter]:
			# Still behind the other one

			if partial:
				yield _partial_yield_val(next_iter)
			last_used[next_iter] = True

		else:
			# Skipped over the other one
			if partial and not last_used[const_iter]:
				yield _partial_yield_val(const_iter)
			last_used[const_iter] = True

	# Yield remaining in other iterator
	if partial:
		for i, val in iters[const_iter]:
			if const_iter == 0:
				yield (i, val, None)
			else:
				yield (i, None, val)
