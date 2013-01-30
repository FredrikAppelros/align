from libcpp cimport bool

cdef extern from 'calign.h':
    cdef struct align_t:
        int s
        size_t len
        short* a1
        short* a2

    align_t align(size_t len_a, short* a, size_t len_b, short* b, int d, size_t len_S, int* S, bool local)

