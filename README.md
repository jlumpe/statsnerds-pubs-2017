# statsnerds-pubs-2017
Stats nerds supergroup in PUBS 2017.


## Requirements

* Python 3.5 or above
* Numpy
* BioPython


## Installation

The package is installable through ``setuptools``, which should come with just about any Python installation. You can install by ``cd``-ing to the directory the repo was cloned into and running

    python setup.py install
    
After this, you can just do

```python
from statsnerds import myfunction
````
    
from a Jupyter notebook or anywhere else.
    
Recommended you do a development install, which means you can modify the files in-place and the changes will be available without having to run the setup script again. Just run:

    python setup.py develop
    
Note that to use your changes you will have to restart your interpreter or IPython/Jupyter kernel. You can also do

```python
import statsnerds
# ...Edit a file...
from importlib import reload
reload(statsnerds)
```
    
but this can be finnicky.
