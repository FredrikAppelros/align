align
=====

Performs sequence alignment. The alignment can either be global or local
in combination with either mutual or non-mutual.

Introduction
------------

align is a Python module that provides a general (read: not restricted to bioinformatics)
implementation of both global and local sequence alignment. The alignments can further be
both mutual and non-mutual. Mutual alignment is sequence alignment where both sequences are
candidates for gap insertion whilst non-mutual alignment only allows for gap insertion into
the second sequence.

It is implemented in C for speed and wrapped with Cython to provide the ease of use of Python.
The algorithms used are Needleman-Wunsch for global alignment and Smith-Waterman for local alignment.

align uses a symmetric ```numpy.ndarray``` as scoring matrix which needs to be of dtype ```numpy.int16```.

It also supports iterative alignments (aligning a sequence with another sequence that already
contains gaps). To use this the scoring matrix needs to be of the shape (257, 257) with
the last row and column being used to score matchings of the fixed gaps and the symbols
of the alphabet.

Installation
------------

Install align with ```python setup.py install```

Usage
-----
```python
>>> import numpy
>>> import align
>>> S = -numpy.ones((256, 256)) + 2 * numpy.identity(256)
>>> S = S.astype(numpy.int16)
>>> s1 = align.string_to_alignment('AATGT')
>>> s2 = align.string_to_alignment('ATGAC')
>>> (s, a1, a2) = align.align(s1, s2, -2, -2, S, local=True)
>>> print('%s\n%s\nScore: %d' % (align.alignment_to_string(a1), align.alignment_to_string(a2), s))
ATG
ATG
Score: 3
```

License
-------

Distributed under the MIT license. See the ```LICENSE``` file.
