from libcpp cimport bool

cdef extern from 'calign.h':
    cdef struct align_t:
        short s
        size_t len
        short* a1
        short* a2

    align_t align(size_t len_a, short* a, size_t len_b, short* b, short d_a, short d_b, size_t len_S, short* S, bool local, bool mutual)

