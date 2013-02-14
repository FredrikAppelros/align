from libc.stdlib cimport malloc, free
from calign cimport align_t
from calign cimport align as calign_
from numpy cimport ndarray

def align(a, b, d, ndarray[short, ndim=2, mode='c'] S not None, local=False):
    """
    Performs local or global sequence alignment.

    Parameters
    ----------
     * a        - first sequence
     * b        - second sequence
     * d        - gap penalty
     * S        - scoring matrix
     * local    - true for local alignment, false otherwise

    Returns
    -------
     * s        - alignment score
     * a1       - alignment of the first sequence
     * a2       - alignment of the second sequence

    """
    cdef size_t len_a = len(a)
    cdef size_t len_b = len(b)
    cdef size_t len_S = len(S)
    cdef size_t i
    cdef align_t al
    cdef short* ca = <short*> malloc(len_a * sizeof(short))
    cdef short* cb = <short*> malloc(len_b * sizeof(short))
    if not ca or not cb:
        raise MemoryError()

    try:
        for i in range(len_a):
            ca[i] = a[i]
        for i in range(len_b):
            cb[i] = b[i]
        al = calign_(len_a, ca, len_b, cb, d, len_S, &S[0,0], local)
    finally:
        free(ca)
        free(cb)

    a1 = [al.a1[i] for i in range(al.len)]
    a2 = [al.a2[i] for i in range(al.len)]

    free(al.a1)
    free(al.a2)

    return (al.s, a1, a2)

def string_to_alignment(s):
    return map(ord, s)

def alignment_to_string(al, hex_=False):
    def conv(c):
        if c != 256:
            if hex_:
                return '%02x' % c
            else:
                return chr(c)
        else:
            if hex_:
                return '--'
            else:
                return '-'
    return ''.join(map(conv, al))

