"""Team Stats Nerds pipeline for PUBS 2017."""

from pkg_resources import resource_string

# Load wild-type sequence
try:
	ASYN_WT_SEQ = resource_string(__name__, 'alphasynuclien_wt_seq_dna.txt').decode()

except FileNotFoundError:
	ASYN_WT_SEQ = None
