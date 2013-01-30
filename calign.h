#include <stdlib.h>
#include <stdbool.h>

typedef struct align_t {
    int s;
    size_t len;
    short* a1;
    short* a2;
} align_t;

align_t align(size_t len_a, const short* a, size_t len_b, const short* b, int d, size_t len_S, const int* S, bool local);
