from libc.stdlib cimport malloc, free
from calign cimport align_t
from calign cimport align as calign_
from numpy cimport ndarray

def align(a, b, d_a, d_b, ndarray[short, ndim=2, mode='c'] S not None, local=False, mutual=True):
    """
    Performs sequence alignment. The alignment can either be global or
    local in combination with either mutual or non-mutual.

    Mutual alignment is sequence alignment where both sequences are
    candidates for gap insertion whilst non-mutual alignment only allows
    for gap insertion into the second sequence.

    NOTE: In non-mutual alignment the longest sequence needs to be
    supplied in sequence a.

    Parameters
    ----------
     * a        - first sequence
     * b        - second sequence
     * d_a      - gap penalty for the first sequence
     * d_b      - gap penalty for the second sequence
     * S        - scoring matrix
     * local    - true for local alignment, false otherwise
     * mutual   - true for mutual alignment, false otherwise

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
    if not mutual:
        assert len_b <= len_a, ("The longest sequence needs to be "
        "supplied in parameter a for non-mutual alignment.")

    try:
        for i in range(len_a):
            ca[i] = a[i]
        for i in range(len_b):
            cb[i] = b[i]
        al = calign_(len_a, ca, len_b, cb, d_a, d_b, len_S, &S[0,0], local, mutual)
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

