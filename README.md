# statsnerds-pubs-2017
Stats nerds supergroup in PUBS 2017.

See wiki for info on installation, package structure, and more.


## Requirements

* Python 3.5 or above
* Numpy
* BioPython


## Example usage

### Barcode analysis

```python

files = ['reads1.fastq', 'reads2.fastq', 'reads3.fastq']

from statsnerds.barcode import process_read_files
results = process_read_files(files)

# Do stuf with your results...
```
