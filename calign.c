#include "calign.h"

#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif

/**

Performs local or global sequence alignment.

Parameters
----------
 * len_a    - length of the first sequence
 * a        - first sequence
 * len_b    - length of the second sequence
 * b        - second sequence
 * d        - gap penalty
 * len_S    - length of alphabet
 * S        - scoring matrix
 * local    - true for local alignment, false otherwise

Returns
-------
 * s        - alignment score
 * len      - length of the alignment
 * a1       - alignment of the first sequence
 * a2       - alignment of the second sequence

*/
align_t align(size_t len_a, const short* a, size_t len_b, const short* b, int d, size_t len_S, const int* S, bool local)
{
    size_t len_al;
    size_t i;
    size_t j;
    size_t max_i = 0;
    size_t max_j = 0;
    size_t cnt1 = 0;
    size_t cnt2 = 0;
    int match;
    int delete;
    int insert;
    int max_val = 0;
    int F[len_b + 1][len_a + 1];
    short tmp;
    short* a1;
    short* a2;
    align_t ret;

    // Initialize the matrix
    if (!local) {
        for (i = 0; i < len_a + 1; i++) {
            F[0][i] = d * i;
        }
        for (i = 0; i < len_b + 1; i++) {
            F[i][0] = d * i;
        }
    } else {
        for (i = 0; i < len_a + 1; i++) {
            F[0][i] = 0;
        }
        for (i = 0; i < len_b + 1; i++) {
            F[i][0] = 0;
        }
    }

    // Fill the matrix
    for (i = 1; i < len_b + 1; i++) {
        for (j = 1; j < len_a + 1; j++) {
            match = F[i-1][j-1] + S[len_S * a[j-1] + b[i-1]];
            delete = F[i-1][j] + d;
            insert = F[i][j-1] + d;
            F[i][j] = max(match, delete);
            F[i][j] = max(F[i][j], insert);
            if (local) {
                F[i][j] = max(F[i][j], 0);
            }
            if (F[i][j] > max_val) {
                max_val = F[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }

    // Backtrace through the matrix
    a1 = malloc((len_a + len_b) * sizeof(short));
    a2 = malloc((len_a + len_b) * sizeof(short));
    if (local) {
        i = max_i;
        j = max_j;
    } else {
        i = len_b;
        j = len_a;
    }
    while (i > 0 || j > 0) {
        if (local && !F[i][j]) {
            break;
        }
        if (i > 0 && j > 0 && F[i][j] == F[i-1][j-1] + S[len_S * a[j-1] + b[i-1]]) {
            a1[cnt1++] = a[j-1];
            a2[cnt2++] = b[i-1];
            i--;
            j--;
        } else if (i > 0 && F[i][j] == F[i-1][j] + d) {
            a1[cnt1++] = -1;
            a2[cnt2++] = b[i-1];
            i--;
        } else {
            a1[cnt1++] = a[j-1];
            a2[cnt2++] = -1;
            j--;
        }
    }

    // Reverse the alignments
    len_al = cnt1;
    for (i = 0; i <= len_al / 2; i++) {
        tmp = a1[i];
        a1[i] = a1[len_al - (i + 1)];
        a1[len_al - (i + 1)] = tmp;
        tmp = a2[i];
        a2[i] = a2[len_al - (i + 1)];
        a2[len_al - (i + 1)] = tmp;
    }

    ret.s = max_val;
    ret.len = len_al;
    ret.a1 = a1;
    ret.a2 = a2;
    return ret;
}

