align
=====

Performs global and local sequence alignment.

Introduction
------------
align is a Python module that provides a general (read: not restricted to bioinformatics)
implementation of both global and local sequence alignment. It is implemented in C for
speed and wrapped with Cython to provide the ease of use of Python. The algorithms used
are Needleman-Wunsch for global alignment and Smith-Waterman for local alignment.

align uses a quadratic ```numpy.ndarray``` as scoring matrix.

Installation
------------

Install align with ```python setup.py install```

Usage
-----
```
```

License
-------
Distributed under the MIT license. See the ```LICENSE``` file.
