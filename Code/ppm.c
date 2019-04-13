/* Courtesy: Nominal Animal <question@nominal-animal.net>
*/

#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include "prng.h"
#include "matrix.h"

static_inline int is_spanning(matrix *const m, const cell  v)
{
    if (m && m->span) {
        const cell *const ends = m->span + m->spans[0] + m->spans[1];
        const cell       *curr = m->span;
        while (curr < ends)
            if (*curr == v)
                return 1;
            else
                curr++;
    }
    return 0;
}

static_inline unsigned int  rgb(const double r, const double g, const double b)
{
    unsigned int  result = 0x000000;

    if (r >= 255.0/256.0)
        result |= 0xFF0000;
    else
    if (r > 0.0)
        result |= (unsigned int)(256.0 * r) << 16;

    if (g >= 255.0/256.0)
        result |= 0x00FF00;
    else
    if (g > 0.0)
        result |= (unsigned int)(256.0 * g) << 8;

    if (b >= 255.0/256.0)
        result |= 0x0000FF;
    else
    if (b > 0.0)
        result |= (unsigned int)(256.0 * b);

    return result;
}

int main(void)
{
    unsigned int      *color;
    prng_seed_string   seed;

    int     n = 100;
    matrix  m = MATRIX_INITIALIZER;
    int     result, r, c;

    result = matrix_init(&m, n, STATS_ALL);
    if (result) {
        fprintf(stderr, "Cannot initialize matrix: %s.\n", matrix_strerror(result));
        return EXIT_FAILURE;
    }

    color = malloc((size_t)n * (size_t)n * 2 * sizeof (unsigned int));
    if (!color) {
        fprintf(stderr, "Not enough memory for color map.\n");
        return EXIT_FAILURE;
    }

    m.nonzero = 0.5;
    m.diagonal = 1.0;
    m.diagonal_nonzero = 0.5;

    fprintf(stderr, "Seed: %s\n", prng_get_seed(&(m.rng), seed));
    fprintf(stderr, "Size: %u x %u cells\n", n, n);

    matrix_generate(&m);

    fprintf(stderr, "White cells: %" PRI_CELL " (%.3f%%)\n", m.fill[0], 100.0 * (double)(m.fill[0]) / (double)(n * n));
    fprintf(stderr, "Black cells: %" PRI_CELL " (%.3f%%)\n", m.fill[1], 100.0 * (double)(m.fill[1]) / (double)(n * n));
    fprintf(stderr, "White clusters: %" PRI_CELL " (%.3f%%)\n", m.unique[0], 100.0 * (double)(m.unique[0]) / (double)(m.unique[0] + m.unique[1]));
    fprintf(stderr, "Black clusters: %" PRI_CELL " (%.3f%%)\n", m.unique[1], 100.0 * (double)(m.unique[1]) / (double)(m.unique[0] + m.unique[1]));

    for (c = 0; c < 2*n*n; c += 2) {
        const double  p = prng_unit(&(m.rng));

        /* White cluster: */
        if (is_spanning(&m, c+0))
            color[c+0] = rgb(0.6 + 0.3*p, 0.6 + 0.3*p, 1.0 - 0.2*p);
        else
            color[c+0] = rgb(0.6 + 0.4*p, 0.6 + 0.4*p, 0.6 + 0.4*p);

        /* Black cluster: */
        if (is_spanning(&m, c+1))
            color[c+1] = rgb(1.0 - 0.4*p, 0.3*p, 0.3*p);
        else
            color[c+1] = rgb(0.4*p, 0.4*p, 0.4*p);
    }
    
    if (m.diagonal > 0.0) {
        fprintf(stderr, "Diagonal cluster joins: %" PRI_CELL " out of %" PRI_CELL " (%.3f%%)\n",
                        m.djoins[0] + m.djoins[1], m.djoins[0] + m.djoins[1] + m.djoins[2],
                        100.0*(double)(m.djoins[0] + m.djoins[1]) / (double)(m.djoins[0] + m.djoins[1] + m.djoins[2]));
        fprintf(stderr, "White clusters joined diagonally: %" PRI_CELL " (%.3f%%)\n",
                        m.djoins[0], 100.0*(double)(m.djoins[0]) / (double)(m.djoins[0] + m.djoins[1]));
        fprintf(stderr, "Black clusters joined diagonally: %" PRI_CELL " (%.3f%%)\n",
                        m.djoins[1], 100.0*(double)(m.djoins[1]) / (double)(m.djoins[0] + m.djoins[1]));
    } else
        fprintf(stderr, "No diagonally joined clusters\n");

    if (m.spans[0] + m.spans[1] > 0)
        fprintf(stderr, "Spanning clusters: %" PRI_CELL " (%" PRI_CELL " or %.3f%% black, %" PRI_CELL " or %.3f%% white)\n",
                        m.spans[0] + m.spans[1],
                        m.spans[1], 100.0 * (double)(m.spans[1]) / (double)(m.spans[0] + m.spans[1]),
                        m.spans[0], 100.0 * (double)(m.spans[0]) / (double)(m.spans[0] + m.spans[1]));
    else
        fprintf(stderr, "No spanning clusters.\n");

    printf("P6\n%d %d\n255\n", n, n);
    for (r = 0; r < n; r++) {
        for (c = 0; c < n; c++) {
            const unsigned int  col = color[m.map[r*m.size + c]];
            fputc((col >> 16) & 255, stdout);
            fputc((col >>  8) & 255, stdout);
            fputc( col        & 255, stdout);
        }
    }

    return EXIT_SUCCESS;
}
    
    
