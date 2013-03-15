#ifndef CALIGN_H
#define CALIGN_H

#include <stdlib.h>
#include <stdbool.h>

typedef struct align_t {
    short s;
    size_t len;
    short* a1;
    short* a2;
} align_t;

align_t align(size_t len_a, const short* a, size_t len_b, const short* b, short d_a, short d_b, size_t len_S, const short* S, bool local, bool mutual);

#endif

