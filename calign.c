#include "calign.h"

#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif

/**

Performs sequence alignment. The alignment can either be global or
local in combination with either mutual or non-mutual.

Mutual alignment is sequence alignment where both sequences are
candidates for gap insertion whilst non-mutual alignment only allows
for gap insertion into the second sequence.

NOTE: In non-mutual alignment the longest sequence needs to be supplied
in sequence a.

Parameters
----------
 * len_a    - length of the first sequence
 * a        - first sequence
 * len_b    - length of the second sequence
 * b        - second sequence
 * d_a      - gap penalty for the first sequence
 * d_b      - gap penalty for the second sequence
 * len_S    - length of alphabet
 * S        - scoring matrix
 * local    - true for local alignment, false otherwise
 * mutual   - true for mutual alignment, false otherwise

Returns
-------
 * s        - alignment score
 * len      - length of the alignment
 * a1       - alignment of the first sequence
 * a2       - alignment of the second sequence

*/
align_t align(size_t len_a, const short* a, size_t len_b, const short* b, short d_a, short d_b, size_t len_S, const short* S, bool local, bool mutual)
{
    size_t len_al;
    size_t i;
    size_t j;
    size_t max_i = 0;
    size_t max_j = 0;
    size_t cnt1 = 0;
    size_t cnt2 = 0;
    short match;
    short delete;
    short insert;
    short max_val = 0;
    short F[len_b + 1][len_a + 1];
    short tmp;
    short* a1;
    short* a2;
    align_t ret;

    // Initialize the matrix
    if (!local) {
        for (i = 0; i < len_a + 1; i++) {
            F[0][i] = d_b * i;
        }
        for (i = 0; i < len_b + 1; i++) {
            F[i][0] = d_a * i;
        }
    } else {
        for (i = 0; i < len_a + 1; i++) {
            F[0][i] = 0;
        }
        for (i = 0; i < len_b + 1; i++) {
            F[i][0] = 0;
        }
    }
    if (!mutual) {
        for (i = 1; i < len_b + 1; i++) {
            F[i][i] = F[i-1][i-1] + S[len_S * a[i-1] + b[i-1]];
        }
    }

    // Fill the matrix
    for (i = 1; i < len_b + 1; i++) {
        for (j = 1 + (mutual ? 0 : i); j < len_a + 1; j++) {
            match = F[i-1][j-1] + S[len_S * a[j-1] + b[i-1]];
            delete = F[i-1][j] + d_a;
            insert = F[i][j-1] + d_b;
            F[i][j] = max(match, insert);
            if (mutual) {
                F[i][j] = max(F[i][j], delete);
            }
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
    if (a1 == NULL || a2 == NULL) {
        exit(1);
    }

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
        } else if (!mutual || j > 0 && F[i][j] == F[i][j-1] + d_b) {
            a1[cnt1++] = a[j-1];
            a2[cnt2++] = 256;
            j--;
        } else {
            a1[cnt1++] = 256;
            a2[cnt2++] = b[i-1];
            i--;
        }
    }

    // Reverse the alignments
    len_al = cnt1;
    for (i = 0; i < len_al / 2; i++) {
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

