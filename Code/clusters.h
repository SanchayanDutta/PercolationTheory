#ifndef   CLUSTERS_H
#define   CLUSTERS_H
#include <stdlib.h>
#include <inttypes.h>
#include <limits.h>
#include <time.h>

/* This file is in public domain. No guarantees, no warranties.
   Written by Nominal Animal <question@nominal-animal.net>.
*/

/* For pure C89 compilers, use '-DSTATIC_INLINE=static' at compile time. */
#ifndef  STATIC_INLINE
#define  STATIC_INLINE  static inline
#endif

#define  ERR_INVALID   -1   /* Invalid function parameter */
#define  ERR_TOOLARGE  -2   /* Matrix size is too large */
#define  ERR_NOMEM     -3   /* Out of memory */

typedef  unsigned char  cluster_color;

typedef  uint32_t  cluster_label;
typedef  uint64_t  cluster_count;

#define  CLUSTER_WHITE  0
#define  CLUSTER_BLACK  1
#define  CLUSTER_NONE   UCHAR_MAX   /* Reserved */

#define  FMT_COLOR  "u"
#define  FMT_LABEL  PRIu32
#define  FMT_COUNT  PRIu64


typedef struct {
    uint64_t        state;
} prng;

typedef struct {
    /* Pseudo-random number generator used */
    prng            rng;

    /* Actual size of the matrix */
    cluster_label   rows;
    cluster_label   cols;

    /* Number of matrices the histograms have been collected from */
    cluster_count   iterations;

    /* Probability of each cell being black */
    uint64_t        p_black;

    /* Probability of a diagonal connection */
    uint64_t        p_diag;
    uint64_t        p_diag_black;

    /* Cluster colormap contains (rows+2) rows and (cols+1) columns */
    cluster_color  *map;

    /* Disjoint set of (rows) rows and (cols) columns */
    cluster_label  *djs;

    /* Number of occurrences per disjoint set root */
    cluster_label  *white_roots;
    cluster_label  *black_roots;

    /* Histograms of white and black clusters */
    cluster_count  *white_histogram;
    cluster_count  *black_histogram;
} cluster;
#define  CLUSTER_INITIALIZER  { {0}, 0, 0, 0, 0.0, 0.0, 0.0, NULL, NULL, NULL, NULL, NULL, NULL }

/* Calculate uint64_t limit corresponding to probability p. */
STATIC_INLINE uint64_t  probability_limit(const double p)
{
    if (p <= 0.0)
        return UINT64_C(0);
    else
    if (p <= 0.5)
        return (uint64_t)(p * 18446744073709551615.0);
    else
    if (p >= 1.0)
        return UINT64_C(18446744073709551615);
    else
        return UINT64_C(18446744073709551615) - (uint64_t)((double)(1.0 - p) * 18446744073709551615.0);
}

/* Return true at probability corresponding to limit 'limit'.
   This implements a Xorshift64* pseudo-random number generator. */
STATIC_INLINE int  probability(prng *const rng, const uint64_t limit)
{
    uint64_t  state = rng->state;
    uint64_t  value;

    /* To correctly cover the entire range, we ensure we never generate a zero. */
    do {
        state ^= state >> 12;
        state ^= state << 25;
        state ^= state >> 27;
        value  = state * UINT64_C(2685821657736338717);
    } while (!value);

    rng->state = state;

    return (value <= limit) ? CLUSTER_BLACK : CLUSTER_WHITE;
}

/* Generate a random seed for the Xorshift64* pseudo-random number generator. */
static uint64_t  randomize(prng *const rng)
{
    unsigned int  rounds = 127;
    uint64_t      state = UINT64_C(3069887672279) * (uint64_t)time(NULL)
                        ^ UINT64_C(60498839) * (uint64_t)clock();
    if (!state)
        state = 1;

    /* Churn the state a bit. */
    while (rounds-->0) {
        state ^= state >> 12;
        state ^= state << 25;
        state ^= state >> 27;
    }

    if (rng)
        rng->state = state;

    return state;
}

/* Free all resources related to a cluster. */
STATIC_INLINE void free_cluster(cluster *c)
{
    if (c) {
        /* Free dynamically allocated pointers. Note: free(NULL) is safe. */
        free(c->map);
        free(c->djs);
        free(c->white_roots);
        free(c->black_roots);
        free(c->white_histogram);
        free(c->black_histogram);
        c->rng.state       = 0;
        c->rows            = 0;
        c->cols            = 0;
        c->iterations      = 0;
        c->p_black         = 0;
        c->p_diag          = 0;
        c->p_diag_black    = 0;
        c->map             = NULL;
        c->djs             = NULL;
        c->white_roots     = 0;
        c->black_roots     = 0;
        c->white_histogram = NULL;
        c->black_histogram = NULL;
    }
}

/* Initialize cluster structure, for a matrix of specified size. */
static int init_cluster(cluster *c, const int rows, const int cols,
                        const double p_black,
                        const double p_diag, const double p_diag_black)
{
    const cluster_label  label_cols = cols;
    const cluster_label  label_rows = rows;
    const cluster_label  color_rows = rows + 2;
    const cluster_label  color_cols = cols + 1;
    const cluster_label  color_cells = color_rows * color_cols;
    const cluster_label  label_cells = label_rows * label_cols;
    const cluster_label  labels = label_cells + 2; /* One extra! */

    if (!c)
        return ERR_INVALID;

    c->rng.state       = 0; /* Invalid seed for Xorshift64*. */
    c->rows            = 0;
    c->cols            = 0;
    c->iterations      = 0;
    c->p_black         = 0;
    c->p_diag          = 0;
    c->p_diag_black    = 0;
    c->map             = NULL;
    c->djs             = NULL;
    c->white_roots     = NULL;
    c->black_roots     = NULL;
    c->white_histogram = NULL;
    c->black_histogram = NULL;

    if (rows < 1 || cols < 1)
        return ERR_INVALID;

    if ((unsigned int)color_rows <= (unsigned int)rows ||
        (unsigned int)color_cols <= (unsigned int)cols ||
        (cluster_label)(color_cells / color_rows) != color_cols ||
        (cluster_label)(color_cells / color_cols) != color_rows ||
        (cluster_label)(label_cells / label_rows) != label_cols ||
        (cluster_label)(label_cells / label_cols) != label_rows)
        return ERR_TOOLARGE;

    c->black_histogram = calloc(labels, sizeof (cluster_count));
    c->white_histogram = calloc(labels, sizeof (cluster_count));
    c->black_roots = calloc(labels, sizeof (cluster_label));
    c->white_roots = calloc(labels, sizeof (cluster_label));
    c->djs = calloc(label_cells, sizeof (cluster_label));
    c->map = calloc(color_cells, sizeof (cluster_color));
    if (!c->map || !c->djs ||
        !c->white_roots || !c->black_roots ||
        !c->white_histogram || !c->black_histogram) {
        free(c->map);
        free(c->djs);
        free(c->white_roots);
        free(c->black_roots);
        free(c->white_histogram);
        free(c->black_histogram);
        return ERR_NOMEM;
    }

    c->rows = rows;
    c->cols = cols;

    c->p_black      = probability_limit(p_black);
    c->p_diag       = probability_limit(p_diag);
    c->p_diag_black = probability_limit(p_diag_black);

    /* Initialize the color map to NONE. */
    {
        cluster_color        *ptr = c->map;
        cluster_color *const  end = c->map + color_cells;
        while (ptr < end)
            *(ptr++) = CLUSTER_NONE;
    }

    /* calloc() initialized the other arrays to zeros already. */
    return 0;
}

/* Disjoint set: find root. */
STATIC_INLINE cluster_label  djs_root(const cluster_label *const  djs, cluster_label  from)
{
    while (from != djs[from])
        from = djs[from];
    return from;
}

/* Disjoint set: path compression. */
STATIC_INLINE void  djs_path(cluster_label *const  djs, cluster_label  from, const cluster_label  root)
{
    while (from != root) {
        const cluster_label  temp = djs[from];
        djs[from] = root;
        from = temp;
    }
}

/* Disjoint set: Flatten. Returns the root, and flattens the path to it. */
STATIC_INLINE cluster_label  djs_flatten(cluster_label *const  djs, cluster_label  from)
{
    const cluster_label  root = djs_root(djs, from);
    djs_path(djs, from, root);
    return root;
}

/* Disjoint set: Join two subsets. */
STATIC_INLINE void  djs_join2(cluster_label *const  djs,
                              cluster_label  from1,  cluster_label  from2)
{
    cluster_label  root, temp;

    root = djs_root(djs, from1);

    temp = djs_root(djs, from2);
    if (root > temp)
        temp = root;

    djs_path(djs, from1, root);
    djs_path(djs, from2, root);
}

/* Disjoint set: Join three subsets. */
STATIC_INLINE void  djs_join3(cluster_label *const  djs,
                              cluster_label  from1, cluster_label  from2,
                              cluster_label  from3)
{
    cluster_label  root, temp;

    root = djs_root(djs, from1);

    temp = djs_root(djs, from2);
    if (root > temp)
        root = temp;

    temp = djs_root(djs, from3);
    if (root > temp)
        root = temp;

    djs_path(djs, from1, root);
    djs_path(djs, from2, root);
    djs_path(djs, from3, root);
}

/* Disjoint set: Join four subsets. */
STATIC_INLINE void  djs_join4(cluster_label *const  djs,
                              cluster_label  from1, cluster_label  from2,
                              cluster_label  from3, cluster_label  from4)
{
    cluster_label  root, temp;

    root = djs_root(djs, from1);

    temp = djs_root(djs, from2);
    if (root > temp)
        root = temp;

    temp = djs_root(djs, from3);
    if (root > temp)
        root = temp;

    temp = djs_root(djs, from4);
    if (root > temp)
        root = temp;

    djs_path(djs, from1, root);
    djs_path(djs, from2, root);
    djs_path(djs, from3, root);
    djs_path(djs, from4, root);
}

static void iterate(cluster *const cl)
{
    prng          *const  rng = &(cl->rng);
    uint64_t       const  p_black = cl->p_black;
    uint64_t       const  p_diag = cl->p_diag;
    uint64_t       const  p_diag_black = cl->p_diag_black;

    cluster_color *const  map = cl->map + cl->cols + 2;
    cluster_label  const  map_stride = cl->cols + 1;

    cluster_label *const  djs = cl->djs;

    cluster_label        *roots[2];

    cluster_label  const  rows = cl->rows;
    cluster_label  const  cols = cl->cols;

    int                   r, c;

    roots[CLUSTER_WHITE] = cl->white_roots;
    roots[CLUSTER_BLACK] = cl->black_roots;

    for (r = 0; r < rows; r++) {
        cluster_label  const  curr_i = r * cols;
        cluster_color *const  curr_row = map + r * map_stride;
        cluster_color *const  prev_row = curr_row - map_stride;

        for (c = 0; c < cols; c++) {
            cluster_color  color = probability(rng, p_black);
            cluster_label  label = curr_i + c;
            unsigned int   joins = 0;

            /* Assign the label and color of the current cell, */
            djs[label] = label;
            curr_row[c] = color;

            /* We use pattern  C B  where X is current cell, and
                               A X
               A, B, C are 1 if they have the same color as X;
               this is then encoded in binary in joins as 00000CBA. */

            joins |= (curr_row[c-1] == color) << 0; /* A */
            joins |= (prev_row[c  ] == color) << 1; /* B */
            joins |= (prev_row[c-1] == color) << 2; /* C */

            /* Do the corresponding joins. */
            switch (joins) {
            case 1: /* X==A */
                djs_join2(djs, label, label - 1);
                break;
            case 2: /* X==B */
                djs_join2(djs, label, label - cols);
                break;
            case 3: /* X==A and X==B */
                djs_join3(djs, label, label - 1, label - cols);
                break;
            case 4: /* X==C and A==B. Diagonal connection case. Handled later. */
                break;
            case 5: /* X==A and X==C. A==C was already connected on previous step. */
                djs_join3(djs, label, label - 1, label - cols - 1);
                break;
            case 6: /* X==B and X==C. B==C was already connected on previous row. */
                djs_join3(djs, label, label - cols, label - cols - 1);
                break;
            case 7: /* X==A, X==B, and X==C. A==C was connected on previous step
                       and B==C on previous row. */
                djs_join4(djs, label, label - 1, label - cols, label - cols - 1);
                break;
            }
        }
    }

    /* Diagonal connections? */
    if (p_diag > 0) {
        /* First, we need to flatten the disjoint sets. */
        {
            size_t  i = rows * cols;
            while (i-->0)
                djs_flatten(djs, i);
        }

        /* Next, we examine the cluster labels (up and left). */
        for (r = 1; r < rows; r++) {
            const cluster_color *const prev_color = map + (r-1)*map_stride;
            const cluster_color *const curr_color = map + r*map_stride;
            const cluster_label        prev_label = (r-1)*cols;
            const cluster_label        curr_label = r*cols;
            for (c = 1; c < cols; c++) {
                /* D C   Color of A is curr_color[c],
                   B>A<  label of A is djs[curr_label+c].
                */
                if (curr_color[c] == prev_color[c-1] &&
                    curr_color[c] != prev_color[c] &&
                    curr_color[c-1] == prev_color[c]) {
                    /* Possible case: A==D, B==C */
                    if (djs[curr_label+c] != djs[prev_label+c-1] &&
                        djs[curr_label+c-1] != djs[prev_label+c]) {
                        /* We can join A==D and B==C. We do join one. */
                        if ((curr_color[c] & 1) == probability(rng, p_diag_black))
                            djs_join2(djs, curr_label + c, prev_label + c - 1);
                        else    
                            djs_join2(djs, curr_label + c - 1, prev_label + c);
                    } else
                    if (djs[curr_label+c] != djs[prev_label+c-1]) {
                        /* We can join A==D. B==C is already joined. */
                        if ((curr_color[c] & 1) == probability(rng, p_diag_black))
                            djs_join2(djs, curr_label + c, prev_label + c - 1);
                    } else
                    if (djs[curr_label+c-1] != djs[prev_label+c]) {
                        /* We can join B==C. A==D is already joined. */
                        if ((curr_color[c-1] & 1) == probability(rng, p_diag_black))
                            djs_join2(djs, curr_label + c - 1, prev_label + c);                
                    } /* Otherwise, both A==D and B==C are already joined. */
                }
            }
        }
    }

    /* Count the occurrences of each disjoint-set root label. */
    if (roots[0] && roots[1]) {
        const size_t  labels = rows * cols + 2;

        /* Clear the counts. */
        memset(roots[0], 0, labels * sizeof (cluster_label));
        memset(roots[1], 0, labels * sizeof (cluster_label));

        for (r = 0; r < rows; r++) {
            const cluster_color *const  curr_row = map + r * map_stride;
            const cluster_label         curr_i   = r * cols;
            for (c = 0; c < cols; c++)
                roots[curr_row[c]][djs_flatten(djs, curr_i + c)]++;
        }
    } else {
        size_t  i = rows * cols;
        while (i-->0)
            djs_flatten(djs, i);
    }

    /* Collect the statistics. */
    if (cl->white_histogram && roots[0]) {
        const cluster_label *const  root_count = roots[0];
        cluster_count *const        histogram = cl->white_histogram;
        size_t                      i = rows * cols + 1;
        while (i-->0)
            histogram[root_count[i]]++;
    }
    if (cl->black_histogram && roots[1]) {
        const cluster_label *const  root_count = roots[1];
        cluster_count *const        histogram = cl->black_histogram;
        size_t                      i = rows * cols + 1;
        while (i-->0)
            histogram[root_count[i]]++;
    }

    /* Note: index zero and (rows*cols+1) are zero in the histogram, for ease of scanning. */
    if (cl->white_histogram || cl->black_histogram) {
        const size_t  n = rows * cols + 1;
        if (cl->white_histogram) {
            cl->white_histogram[0] = 0;
            cl->white_histogram[n] = 0;
        }
        if (cl->black_histogram) {
            cl->black_histogram[0] = 0;
            cl->black_histogram[n] = 0;
        }
    }

    /* One more iteration performed. */
    cl->iterations++;
}

#endif /* CLUSTERS_H */
