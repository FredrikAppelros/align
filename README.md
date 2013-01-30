align
=====

Performs global and local sequence alignment.

Introduction
------------

align is a Python module that provides a general (read: not restricted to bioinformatics)
implementation of both global and local sequence alignment. It is implemented in C for
speed and wrapped with Cython to provide the ease of use of Python. The algorithms used
are Needleman-Wunsch for global alignment and Smith-Waterman for local alignment.

align uses a quadratic ```numpy.ndarray``` as scoring matrix which needs to be of dtype ```numpy.int32```.

Installation
------------

Install align with ```python setup.py install```

Usage
-----
```
python
>>> import numpy
>>> import align
>>> S = -numpy.ones((256, 256)) + 2 * numpy.identity(256)
>>> S = S.astype(numpy.int32)
>>> (s, a1, a2) = align.align('AATGT', 'ATGAC', -2, S, local=True)
>>> print('%s\n%s\nScore: %d' % (align.alignment_to_string(a1), align.alignment_to_string(a2), s))
ATG
ATG
Score: 3
```

License
-------

Distributed under the MIT license. See the ```LICENSE``` file.
