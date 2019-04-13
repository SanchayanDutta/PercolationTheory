#ifndef   MATRIX_H
#define   MATRIX_H
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include "prng.h"

/* Written by Nominal Animal <question@nominal-animal.net>.
   This is in public domain. No guarantees, no warranties. */

#ifndef  static_inline
#define  static_inline  static inline
#endif

typedef  uint64_t   count;
#define  PRI_COUNT  PRIu64
#define  SCN_COUNT  SCNu64

typedef  uint16_t   cell;
#define  PRI_CELL   PRIu16
#define  SCN_CELL   SCNu16

#define  CELL_INDEX(v)      ((cell)(v) >> 1)
#define  CELL_COLOR(v)      ((cell)(v) & 1)
#define  CELL_VALUE(i, c)   (((cell)(i) << 1) | ((cell)(c) & 1))
#define  SAME_COLOR(v1, v2) (!(((cell)(v1) ^ (cell)(v2)) & 1))

static_inline int cell_cmp(const void *ptr1, const void *ptr2)
{
    const cell  val1 = *(const cell *)ptr1;
    const cell  val2 = *(const cell *)ptr2;
    return (val1 < val2) ? -1 :
           (val1 > val2) ? +1 : 0;
}

static_inline size_t cell_common(cell *const dest, cell *const array1, cell *const array2, const size_t count)
{
    const cell       *src1 = array1;
    const cell *const end1 = array1 + count;
    const cell       *src2 = array2;
    const cell *const end2 = array2 + count;
    size_t            have = 0;

    qsort(array1, count, sizeof (cell), cell_cmp);
    qsort(array2, count, sizeof (cell), cell_cmp);

    while (src1 < end1 && src2 < end2) {
        while (src1 < end1 && *src1 < *src2)
            src1++;
        while (src2 < end2 && *src2 < *src1)
            src2++;

        while (src1 < end1 && src2 < end2 && *src1 == *src2) {
            const cell  same = *src1;

            do {
                src1++;
            } while (src1 < end1 && *src1 == same);

            do {
                src2++;
            } while (src2 < end2 && *src2 == same);

            dest[have++] = same;
        }
    }

    return have;
}

typedef struct {
    prng        rng;
    size_t      size;
    cell       *map;                 /* Cell map, size*size cells */
    double      nonzero;             /* Probability of a cell to be nonzero */
    double      diagonal;            /* Probability of connecting clusters diagonally */
    double      diagonal_nonzero;    /* Probability of diagonal connection being between nonzero clusters */
    cell        unique[2];           /* Number of unique clusters. Only if counts is available. */
    cell        fill[2];             /* Number of cells of each color */
    cell        spans[2];            /* Number of zero/nonzero spanning clusters */
    cell       *span;                /* 3*size array for spanning testing */
    cell       *counts;              /* Cluster value occurrences, (size*size)*2 */
    cell        djoins[3];           /* Number of diagonal joins. [2] is omitted joins. Only updated if diagonal > 0. */
} matrix;
#define  MATRIX_INITIALIZER  { {0}, 0, }

static_inline void matrix_free(matrix *m)
{
    if (m) {
        free(m->counts);
        free(m->span);
        free(m->map);
        m->size   = 0;
        m->map    = NULL;
        m->span   = NULL;
        m->counts = NULL;
    }
}

static_inline const char *matrix_strerror(const int retval)
{
    switch (retval) {
    case 0: return "OK";
    case 1: return "No matrix specified";
    case 2: return "Invalid matrix size";
    case 3: return "Matrix size is too large";
    case 4: return "Out of memory";
    default: return "(Unknown error)";
    }
}

#define  STATS_NONE        0
#define  STATS_ALL        (~0u)
#define  STATS_SPANNING   (1u << 0)
#define  STATS_CLUSTERS   (1u << 1)

int matrix_init(matrix *m, const size_t size, const unsigned int statistics)
{
    const cell    test = CELL_VALUE(size * size, 0);
    const size_t  cells = CELL_INDEX(test);
    const size_t  twice = size * size * 2;

    if (!m)
        return 1; /* No matrix specified */

    /* Initialize fields, so we can safely call matrix_free() for cleanup. */
    m->size   = 0;
    m->map    = NULL;
    m->span   = NULL;
    m->counts = NULL;

    if (size < 2)
        return 2; /* Invalid size */

    if (cells / size != size ||
        (size_t)((twice / size) / 2) != size)
        return 3; /* Size is too large */

    m->map = calloc(cells, sizeof (cell));
    if (!m->map)
        return 4; /* Not enough memory */

    prng_init(&(m->rng));
    memset(m->spans, 0, sizeof m->spans);

    m->size             = size;
    m->nonzero          = 0.5;
    m->diagonal         = 0.0;
    m->diagonal_nonzero = 0.0;

    if (statistics & STATS_SPANNING) {
        m->span = calloc(3 * size, sizeof (cell));
        if (!m->span) {
            matrix_free(m);
            return 4; /* Not enough memory */
        }
    }

    if (statistics & STATS_CLUSTERS) {
        m->counts = calloc(twice, sizeof (cell));
        if (!m->counts) {
            matrix_free(m);
            return 4; /* Not enough memory */
        }
    }

    return 0;
}

/* Disjoint set operations on cells. Define VERIFY to add color checks. */

static_inline cell  djs_root(cell *const djs, size_t index)
{
    cell  prev, curr = djs[index];

    do {
        index = CELL_INDEX(curr);
        prev = curr;
        curr = djs[index];
#ifdef VERIFY
        if (CELL_COLOR(curr) != CELL_COLOR(prev)) {
            fprintf(stderr, "djs_root(): Cluster contains cells with different colors.\n");
            exit(EXIT_FAILURE);
        }
#endif
    } while (prev != curr);

    return curr;
}

static_inline void  djs_path(cell *const djs, size_t index, const cell root)
{
    cell  curr = djs[index];

    while (curr != root) {
#ifdef VERIFY
        if (CELL_COLOR(curr) != CELL_COLOR(root)) {
            fprintf(stderr, "djs_path(): Cluster contains cells with different colors.\n");
            exit(EXIT_FAILURE);
        }
#endif
        djs[index] = root;
        index = CELL_INDEX(curr);
        curr = djs[index];
    }
}

static_inline cell  djs_flatten(cell *const djs, size_t index)
{
    cell  root;
    root = djs_root(djs, index);
    djs_path(djs, index, root);
    return root;
}

static_inline cell  djs_join2(cell *const djs, size_t index1, size_t index2)
{
    cell  root, temp;

    root = djs_root(djs, index1);

    temp = djs_root(djs, index2);
#ifdef VERIFY
    if (CELL_COLOR(temp) != CELL_COLOR(root)) {
        fprintf(stderr, "djs_join2(): Cluster contains cells with different colors: index %" PRI_CELL " has color %" PRI_CELL ", but index %" PRI_CELL " color %" PRI_CELL ".\n",
                        CELL_INDEX(temp), CELL_COLOR(temp), CELL_INDEX(root), CELL_COLOR(root));
        exit(EXIT_FAILURE);
    }
#endif
    if (root > temp)
        root = temp;

    djs_path(djs, index1, root);
    djs_path(djs, index2, root);

    return root;
}    

static_inline cell  djs_join3(cell *const djs, size_t index1, size_t index2, size_t index3)
{
    cell  root, temp;

    root = djs_root(djs, index1);

    temp = djs_root(djs, index2);
#ifdef VERIFY
    if (CELL_COLOR(temp) != CELL_COLOR(root)) {
        fprintf(stderr, "djs_join3(): Cluster contains cells with different colors.\n");
        exit(EXIT_FAILURE);
    }
#endif
    if (root > temp)
        root = temp;

    temp = djs_root(djs, index3);
#ifdef VERIFY
    if (CELL_COLOR(temp) != CELL_COLOR(root)) {
        fprintf(stderr, "djs_join3(): Cluster contains cells with different colors.\n");
        exit(EXIT_FAILURE);
    }
#endif
    if (root > temp)
        root = temp;

    djs_path(djs, index1, root);
    djs_path(djs, index2, root);
    djs_path(djs, index3, root);

    return root;
}    

static_inline cell  djs_join4(cell *const djs, size_t index1, size_t index2, size_t index3, size_t index4)
{
    cell  root, temp;

    root = djs_root(djs, index1);

    temp = djs_root(djs, index2);
#ifdef VERIFY
    if (CELL_COLOR(temp) != CELL_COLOR(root)) {
        fprintf(stderr, "djs_join4(): Cluster contains cells with different colors.\n");
        exit(EXIT_FAILURE);
    }
#endif
    if (root > temp)
        root = temp;

    temp = djs_root(djs, index3);
#ifdef VERIFY
    if (CELL_COLOR(temp) != CELL_COLOR(root)) {
        fprintf(stderr, "djs_join4(): Cluster contains cells with different colors.\n");
        exit(EXIT_FAILURE);
    }
#endif
    if (root > temp)
        root = temp;

    temp = djs_root(djs, index4);
#ifdef VERIFY
    if (CELL_COLOR(temp) != CELL_COLOR(root)) {
        fprintf(stderr, "djs_join4(): Cluster contains cells with different colors.\n");
        exit(EXIT_FAILURE);
    }
#endif
    if (root > temp)
        root = temp;

    djs_path(djs, index1, root);
    djs_path(djs, index2, root);
    djs_path(djs, index3, root);
    djs_path(djs, index4, root);

    return root;
}


void matrix_generate(matrix *const m)
{
    prng *const       rng = &(m->rng);
    cell *const       map = m->map;
    const size_t      size = m->size;
    const prng_limit  p_1 = prng_set_probability(m->nonzero);

    /* First row. */
    {
        cell    prevvalue, prevcolor, currvalue, currcolor;
        size_t  c;
        
        currcolor = prng_probability(rng, p_1);
        map[0] = currvalue = CELL_VALUE(0, currcolor);
        for (c = 1; c < size; c++) {
            prevcolor = currcolor;
            prevvalue = currvalue;
            currcolor = prng_probability(rng, p_1);
            map[c] = currvalue = ((prevcolor == currcolor) ? prevvalue : CELL_VALUE(c, currcolor));
        }
    }

    /* Other rows. */
    {
        size_t  r, index;

        for (r = 1; r < size; r++) {
            const size_t  endindex = r * size + size;
            map[r*size] = CELL_VALUE(r*size, prng_probability(rng, p_1));

            for (index = r * size + 1; index < endindex; index++) {
                const cell  color = prng_probability(rng, p_1);

                switch ( ((CELL_COLOR(map[index-1]) == color) ? 1 : 0)
                       + ((CELL_COLOR(map[index-size]) == color) ? 2 : 0) ) {
                case 0: /* Different color than left or up. */
                    map[index] = CELL_VALUE(index, color); break;
                case 1: /* Join left. */
                    map[index] = djs_flatten(map, index-1); break;
                case 2: /* Join up. */
                    map[index] = djs_flatten(map, index-size); break;
                case 3: /* Join up and left. */
                    map[index] = djs_join2(map, index-1, index-size); break;
                } 
            }
        }
    }

    /* Diagonal connection pass? */
    if (m->diagonal > 0.0) {
        const prng_limit  p_d = prng_set_probability(m->diagonal);
        const prng_limit  p_d_1 = prng_set_probability(m->diagonal_nonzero);
        const size_t      last = size - 1;
        size_t            joins[3] = { 0, 0, 0 };
        size_t            r, index;
        cell              value;

        for (r = 0; r < last; r++) {
            const size_t  endindex = r * size + last;
            for (index = r*size; index < endindex; index++) {
                const cell  target    = djs_flatten(map, index);
                const cell  right     = djs_flatten(map, index + 1);
                const cell  down      = djs_flatten(map, index + size);
                const cell  downright = djs_flatten(map, index + size + 1);

                if (target != downright &&
                    right != down &&
                    SAME_COLOR(target, downright) &&
                    SAME_COLOR(right, down) &&
                    !SAME_COLOR(target, right)) {
                    /* Possible diagonal connection case. */
                    if (prng_probability(rng, p_d)) {
                        /* Connect diagonally. */
                        if (prng_probability(rng, p_d_1) == CELL_COLOR(target))
                            value = djs_join2(map, index, index + size + 1);
                        else
                            value = djs_join2(map, index + 1, index + size);
                        /* Update diagonal count based on color. */
                        joins[CELL_COLOR(value)]++;
                    } else {
                        /* Diagonal connection was possible, but was omitted. */
                        joins[2]++;
                    }
                }
            }
        }

        /* Save counters. */
        m->djoins[0] = joins[0];
        m->djoins[1] = joins[1];
        m->djoins[2] = joins[2];
    } else {
        /* No diagonal joining; clear counters. */
        m->djoins[0] = 0;
        m->djoins[1] = 0;
        m->djoins[2] = 0;
    }

    /* Flatten clusters. */
    if (m->counts) {
        cell *const  counts = m->counts;
        const size_t total = 2*size*size;
        size_t       unique[2] = { 0, 0 };
        size_t       fill[2] = { 0, 0 };
        size_t       index;

        index = total;
        while (index-->0)
            counts[index] = 0;

        index = size*size;
        while (index-->0) {
            const cell  c = djs_flatten(map, index);
#ifdef VERIFY
            if ((size_t)c >= total) {
                fprintf(stderr, "matrix_generate(): Invalid cell value (%" PRI_CELL ", maximum %lu).\n", c, (unsigned long)(total-1));
                exit(EXIT_FAILURE);
            }
#endif
            unique[CELL_COLOR(c)] += !(counts[c]++);
            fill[CELL_COLOR(c)]++;
        }

        m->unique[0] = unique[0];
        m->unique[1] = unique[1];
        m->fill[0] = fill[0];
        m->fill[1] = fill[1];

    } else {
#ifdef VERIFY
        const size_t  total = 2 * size * size;
#endif
        size_t        index = size * size;
        size_t        fill[2] = { 0, 0 };

        while (index-->0) {
            const cell  c = djs_flatten(map, index);
#ifdef VERIFY
            if ((size_t)c >= total) {
                fprintf(stderr, "matrix_generate(): Invalid cell value (%" PRI_CELL ", maximum %lu).\n", c, (unsigned long)(total-1));
                exit(EXIT_FAILURE);
            }
#endif
            fill[CELL_COLOR(c)]++;
        }

        m->fill[0] = fill[0];
        m->fill[1] = fill[1];
        m->unique[0] = 0;
        m->unique[1] = 0;
    }

    /* Spanning test? */
    if (m->span) {
        cell *const span = m->span;
        size_t spans[2] = { 0, 0 };
        size_t n;

        /* Because the matrix is square, we can grab a minor speedup by
           checking for vertical spanning. */
        memcpy(span + 1*size, map,                 size * sizeof (cell));
        memcpy(span + 2*size, map + (size-1)*size, size * sizeof (cell));
        n = cell_common(span, span + size, span + 2*size, size);

        /* In case of further user analysis, we append ~(cell)0 to the list. */
        span[n] = ~(cell)0;

        /* Count spanning clusters by color. */
        while (n-->0)
            spans[CELL_COLOR(span[n])]++;

        /* Update counts. */
        m->spans[0] = spans[0];
        m->spans[1] = spans[1];
    } else {
        /* No spanning tests; clear counters. */
        m->spans[0] = 0;
        m->spans[1] = 0;
    }

    /* Done. */
}

#endif /* MATRIX_H */
