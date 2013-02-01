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
            ca[i] = ord(a[i])
        for i in range(len_b):
            cb[i] = ord(b[i])
        al = calign_(len_a, ca, len_b, cb, d, len_S, &S[0,0], local)
    finally:
        free(ca)
        free(cb)

    a1 = []
    a2 = []
    for i in range(al.len):
        a1i = al.a1[i]
        a2i = al.a2[i]
        if a1i != -1:
            a1.append(a1i)
        else:
            a1.append(None)
        if a2i != -1:
            a2.append(a2i)
        else:
            a2.append(None)

    free(al.a1)
    free(al.a2)

    return (al.s, a1, a2)

def alignment_to_string(al):
    def conv(c):
        if c is not None:
            return chr(c)
        return '-'
    return ''.join(map(conv, al))

