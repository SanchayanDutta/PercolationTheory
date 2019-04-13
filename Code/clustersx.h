#ifndef   CLUSTERS_H
#define   CLUSTERS_H
#include <stdlib.h>
#include <inttypes.h>
#include <limits.h>
#include <time.h>
#include <string.h>
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

	/* Number of times when at least one cluster spanned the matrix */
	cluster_count   white_spans;
	cluster_count   black_spans;

	/* Probability of each cell being black */
	uint64_t        p_black;

	/* Probability of diagonal connections */
	uint64_t        d_black;
	uint64_t        d_white;

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
#define  CLUSTER_INITIALIZER  { {0}, 0, 0, 0, 0, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL }

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
		value = state * UINT64_C(2685821657736338717);
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
		c->rng.state = 0;
		c->rows = 0;
		c->cols = 0;
		c->iterations = 0;
		c->white_spans = 0;
		c->black_spans = 0;
		c->p_black = 0;
		c->d_white = 0;
		c->d_black = 0;
		c->map = NULL;
		c->djs = NULL;
		c->white_roots = 0;
		c->black_roots = 0;
		c->white_histogram = NULL;
		c->black_histogram = NULL;
	}
}

/* Initialize cluster structure, for a matrix of specified size. */
static int init_cluster(cluster *c, const int rows, const int cols,
	const double p_black,
	const double d_white, const double d_black)
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

	c->rng.state = 0; /* Invalid seed for Xorshift64*. */
	c->rows = 0;
	c->cols = 0;
	c->iterations = 0;
	c->white_spans = 0;
	c->black_spans = 0;
	c->p_black = 0;
	c->d_white = 0;
	c->d_black = 0;
	c->map = NULL;
	c->djs = NULL;
	c->white_roots = NULL;
	c->black_roots = NULL;
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

	c->black_histogram = (cluster_count*)calloc(labels, sizeof(cluster_count));
	c->white_histogram = (cluster_count*)calloc(labels, sizeof(cluster_count));
	c->black_roots = (cluster_label*)calloc(labels, sizeof(cluster_label));
	c->white_roots = (cluster_label*)calloc(labels, sizeof(cluster_label));
	c->djs = (cluster_label*)calloc(label_cells, sizeof(cluster_label));
	c->map = (cluster_color*)calloc(color_cells, sizeof(cluster_color));
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

	c->p_black = probability_limit(p_black);
	c->d_white = probability_limit(d_white);
	c->d_black = probability_limit(d_black);

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
	cluster_label  from1, cluster_label  from2)
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

/* Disjoint set: Join five subsets. */
STATIC_INLINE void  djs_join5(cluster_label *const  djs,
	cluster_label  from1, cluster_label  from2,
	cluster_label  from3, cluster_label  from4,
	cluster_label  from5)
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

	temp = djs_root(djs, from5);
	if (root > temp)
		root = temp;

	djs_path(djs, from1, root);
	djs_path(djs, from2, root);
	djs_path(djs, from3, root);
	djs_path(djs, from4, root);
	djs_path(djs, from5, root);
}

static int compare_cluster_label(const void *ptr1, const void *ptr2)
{
	const cluster_label  val1 = *(const cluster_label *)ptr1;
	const cluster_label  val2 = *(const cluster_label *)ptr2;
	return (val1 < val2) ? -1 :
		(val1 > val2) ? +1 : 0;
}

/* Return nonzero if the two arrays contain any same elements. */
static int have_same_labels(cluster_label *set1, cluster_label *const end1,
	cluster_label *set2, cluster_label *const end2)
{
	if (end1 <= set1 || end2 <= set2)
		return 0;

	/* Sort the arrays. */
	if (set1 + 1 < end1)
		qsort(set1, (size_t)(end1 - set1), sizeof set1[0], compare_cluster_label);
	if (set2 + 1 < end2)
		qsort(set2, (size_t)(end2 - set2), sizeof set2[0], compare_cluster_label);

	while (set1 < end1 && set2 < end2)
		if (*set1 < *set2) {
			while (set1 < end1 && *set1 < *set2)
				set1++;
		}
		else
			if (*set2 < *set1) {
				while (set2 < end2 && *set2 < *set1)
					set2++;
			}
			else
				return 1;

	return 0;
}


static void iterate(cluster *const cl)
{
	prng          *const  rng = &(cl->rng);
	uint64_t       const  p_black = cl->p_black;
	uint64_t              d_color[2];

	cluster_color *const  map = cl->map + cl->cols + 2;
	cluster_label  const  map_stride = cl->cols + 1;

	cluster_label *const  djs = cl->djs;

	cluster_label        *roots[2];

	cluster_label  const  rows = cl->rows;
	cluster_label  const  cols = cl->cols;

	int                   r, c;

	d_color[CLUSTER_WHITE] = cl->d_white;
	d_color[CLUSTER_BLACK] = cl->d_black;
	roots[CLUSTER_WHITE] = cl->white_roots;
	roots[CLUSTER_BLACK] = cl->black_roots;

	for (r = 0; r < rows; r++) {
		cluster_label  const  curr_i = r * cols;
		cluster_color *const  curr_row = map + r * map_stride;
		cluster_color *const  prev_row = curr_row - map_stride;

		for (c = 0; c < cols; c++) {
			cluster_color  color = probability(rng, p_black);
			cluster_label  label = curr_i + c;
			uint64_t       diag = d_color[color];
			unsigned int   joins = 0;

			/* Assign the label and color of the current cell, */
			djs[label] = label;
			curr_row[c] = color;

			/* Because we join left, up-left, up, and up-right, and
			all those have been assigned to already, we can do
			the necessary joins right now. */

			/* Join left? */
			joins |= (curr_row[c - 1] == color) << 0;

			/* Join up? */
			joins |= (prev_row[c] == color) << 1;

			/* Join up left? */
			joins |= (prev_row[c - 1] == color && probability(rng, diag)) << 2;

			/* Join up right? */
			joins |= (prev_row[c + 1] == color && probability(rng, diag)) << 3;

			/* Do the corresponding joins. */
			switch (joins) {
			case 1: /* Left */
				djs_join2(djs, label, label - 1);
				break;
			case 2: /* Up */
				djs_join2(djs, label, label - cols);
				break;
			case 3: /* Left and up */
				djs_join3(djs, label, label - 1, label - cols);
				break;
			case 4: /* Up-left */
				djs_join2(djs, label, label - cols - 1);
				break;
			case 5: /* Left and up-left */
				djs_join3(djs, label, label - 1, label - cols - 1);
				break;
			case 6: /* Up and up-left */
				djs_join3(djs, label, label - cols, label - cols - 1);
				break;
			case 7: /* Left, up, and up-left */
				djs_join4(djs, label, label - 1, label - cols, label - cols - 1);
				break;
			case 8: /* Up-right */
				djs_join2(djs, label, label - cols + 1);
				break;
			case 9: /* Left and up-right */
				djs_join3(djs, label, label - 1, label - cols + 1);
				break;
			case 10: /* Up and up-right */
				djs_join3(djs, label, label - cols, label - cols + 1);
				break;
			case 11: /* Left, up, and up-right */
				djs_join4(djs, label, label - 1, label - cols, label - cols + 1);
				break;
			case 12: /* Up-left and up-right */
				djs_join3(djs, label, label - cols - 1, label - cols + 1);
				break;
			case 13: /* Left, up-left, and up-right */
				djs_join4(djs, label, label - 1, label - cols - 1, label - cols + 1);
				break;
			case 14: /* Up, up-left, and up-right */
				djs_join4(djs, label, label - cols, label - cols - 1, label - cols + 1);
				break;
			case 15: /* Left, up, up-left, and up-right */
				djs_join5(djs, label, label - 1, label - cols, label - cols - 1, label - cols + 1);
				break;
			}
		}
	}

	/* Count the occurrences of each disjoint-set root label. */
	if (roots[0] && roots[1]) {
		const size_t  labels = rows * cols + 2;

		/* Clear the counts. */
		memset(roots[0], 0, labels * sizeof(cluster_label));
		memset(roots[1], 0, labels * sizeof(cluster_label));

		for (r = 0; r < rows; r++) {
			const cluster_color *const  curr_row = map + r * map_stride;
			const cluster_label         curr_i = r * cols;
			for (c = 0; c < cols; c++)
				roots[curr_row[c]][djs_flatten(djs, curr_i + c)]++;
		}
	}
	else {
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

	/* Check if a white cluster spans the matrix. */
	{
		cluster_label *const  set1 = roots[0];
		cluster_label        *end1 = set1;
		cluster_label *const  set2 = roots[1];
		cluster_label        *end2 = set2;
		cluster_color *const  map_lastrow = map + (rows - 1)*map_stride;
		cluster_label *const  djs_lastrow = djs + (rows - 1)*cols;

		for (r = 0; r < rows; r++) {
			if (map[r*map_stride] == CLUSTER_WHITE)
				*(end1++) = djs[r*cols];
			if (map[r*map_stride + cols - 1] == CLUSTER_WHITE)
				*(end2++) = djs[r*cols + cols - 1];
		}

		if (have_same_labels(set1, end1, set2, end2))
			cl->white_spans++;
		else {
			end1 = set1;
			end2 = set2;

			for (c = 0; c < cols; c++) {
				if (map[c] == CLUSTER_WHITE)
					*(end1++) = djs[c];
				if (map_lastrow[c] == CLUSTER_WHITE)
					*(end2++) = djs_lastrow[c];
			}

			if (have_same_labels(set1, end1, set2, end2))
				cl->white_spans++;
		}
	}

	/* Check if a black cluster spans the matrix. */
	{
		cluster_label *const  set1 = roots[0];
		cluster_label        *end1 = set1;
		cluster_label *const  set2 = roots[1];
		cluster_label        *end2 = set2;
		cluster_color *const  map_lastrow = map + (rows - 1)*map_stride;
		cluster_label *const  djs_lastrow = djs + (rows - 1)*cols;

		for (r = 0; r < rows; r++) {
			if (map[r*map_stride] == CLUSTER_BLACK)
				*(end1++) = djs[r*cols];
			if (map[r*map_stride + cols - 1] == CLUSTER_BLACK)
				*(end2++) = djs[r*cols + cols - 1];
		}

		if (have_same_labels(set1, end1, set2, end2))
			cl->black_spans++;
		else {
			end1 = set1;
			end2 = set2;

			for (c = 0; c < cols; c++) {
				if (map[c] == CLUSTER_BLACK)
					*(end1++) = djs[c];
				if (map_lastrow[c] == CLUSTER_BLACK)
					*(end2++) = djs_lastrow[c];
			}

			if (have_same_labels(set1, end1, set2, end2))
				cl->black_spans++;
		}
	}

	/* One more iteration performed. */
	cl->iterations++;
}

#endif /* CLUSTERS_H */
