Courtesy: Nominal Animal <question@nominal-animal.net>

#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <stdio.h>
#include "clusters.h"

#define  DEFAULT_ROWS          100
#define  DEFAULT_COLS          100
#define  DEFAULT_P_BLACK       0.0
#define  DEFAULT_P_DIAG        0.0
#define  DEFAULT_P_DIAG_BLACK  0.0
#define  DEFAULT_ITERS         1

int usage(const char *argv0)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: %s [ -h | --help ]\n", argv0);
    fprintf(stderr, "       %s OPTIONS [ > output.txt ]\n", argv0);
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "       rows=SIZE    Set number of rows. Default is %d.\n", DEFAULT_ROWS);
    fprintf(stderr, "       cols=SIZE    Set number of columns. Default is %d.\n", DEFAULT_ROWS);
    fprintf(stderr, "       L=SIZE       Set rows=SIZE and cols=SIZE.\n");
    fprintf(stderr, "       black=P      Set the probability of a cell to be black. Default is %g.\n", DEFAULT_P_BLACK);
    fprintf(stderr, "                    All non-black cells are white.\n");
    fprintf(stderr, "       diag=P       Set the probability of cells connecting diagonally.\n");
    fprintf(stderr, "                    Default is %g.\n", DEFAULT_P_DIAG);
    fprintf(stderr, "       diagblack=P  Set the probability of a diagonal connection to be between\n");
    fprintf(stderr, "                    black cells. Default is %g.\n", DEFAULT_P_DIAG_BLACK);
    fprintf(stderr, "       N=COUNT      Number of iterations for gathering statistics. Default is %d.\n", DEFAULT_ITERS);
    fprintf(stderr, "       seed=U64     Set the Xorshift64* pseudorandom number generator seed; nonzero.\n");
    fprintf(stderr, "                    Default is to pick one randomly (based on time).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "The output consists of comment lines and data lines.\n");
    fprintf(stderr, "Comment lines begin with a #:\n");
    fprintf(stderr, "   # This is a comment line.\n");
    fprintf(stderr, "Each data line contains a cluster size, the number of white clusters of that size\n");
    fprintf(stderr, "observed during iterations, the number of black clusters of that size observed\n");
    fprintf(stderr, "during iterations, and the number of any clusters of that size observed:\n");
    fprintf(stderr, "   SIZE  WHITE_CLUSTERS  BLACK_CLUSTERS  TOTAL_CLUSTERS\n");
    fprintf(stderr, "\n");
    return EXIT_SUCCESS;
}

int main(int argc, char *argv[])
{
    int      rows = DEFAULT_ROWS;
    int      cols = DEFAULT_COLS;
    double   p_black      = DEFAULT_P_BLACK;
    double   p_diag       = DEFAULT_P_DIAG;
    double   p_diag_black = DEFAULT_P_DIAG_BLACK;
    long     iters = DEFAULT_ITERS;
    uint64_t seed = 0;
    cluster  c = CLUSTER_INITIALIZER;

    int      arg, itemp;
    uint64_t u64temp;
    double   dtemp;
    long     ltemp;
    char     dummy;

    size_t   n;
    size_t   i;

    if (argc < 2)
        return usage(argv[0]);

    for (arg = 1; arg < argc; arg++)
        if (!strcmp(argv[arg], "-h") || !strcmp(argv[arg], "/?") || !strcmp(argv[arg], "--help"))
            return usage(argv[0]);
        else
        if (sscanf(argv[arg], "L=%d %c", &itemp, &dummy) == 1 ||
            sscanf(argv[arg], "l=%d %c", &itemp, &dummy) == 1 ||
            sscanf(argv[arg], "size=%d %c", &itemp, &dummy) == 1) {
            rows = itemp;
            cols = itemp;
        } else
        if (sscanf(argv[arg], "seed=%" SCNu64 " %c", &u64temp, &dummy) == 1 ||
            sscanf(argv[arg], "seed=%" SCNx64 " %c", &u64temp, &dummy) == 1 ||
            sscanf(argv[arg], "s=%" SCNu64 " %c", &u64temp, &dummy) == 1 ||
            sscanf(argv[arg], "s=%" SCNx64 " %c", &u64temp, &dummy) == 1) {
            seed = u64temp;
        } else
        if (sscanf(argv[arg], "N=%ld %c", &ltemp, &dummy) == 1 ||
            sscanf(argv[arg], "n=%ld %c", &ltemp, &dummy) == 1 ||
            sscanf(argv[arg], "count=%ld %c", &ltemp, &dummy) == 1) {
            iters = ltemp;
        } else
        if (sscanf(argv[arg], "rows=%d %c", &itemp, &dummy) == 1 ||
            sscanf(argv[arg], "r=%d %c", &itemp, &dummy) == 1 ||
            sscanf(argv[arg], "height=%d %c", &itemp, &dummy) == 1 ||
            sscanf(argv[arg], "h=%d %c", &itemp, &dummy) == 1) {
            rows = itemp;
        } else
        if (sscanf(argv[arg], "columns=%d %c", &itemp, &dummy) == 1 ||
            sscanf(argv[arg], "cols=%d %c", &itemp, &dummy) == 1 ||
            sscanf(argv[arg], "c=%d %c", &itemp, &dummy) == 1 ||
            sscanf(argv[arg], "width=%d %c", &itemp, &dummy) == 1 ||
            sscanf(argv[arg], "w=%d %c", &itemp, &dummy) == 1) {
            cols = itemp;
        } else
        if (sscanf(argv[arg], "black=%lf %c", &dtemp, &dummy) == 1 ||
            sscanf(argv[arg], "p0=%lf %c", &dtemp, &dummy) == 1 ||
            sscanf(argv[arg], "b=%lf %c", &dtemp, &dummy) == 1 ||
            sscanf(argv[arg], "P=%lf %c", &dtemp, &dummy) == 1 ||
            sscanf(argv[arg], "p0=%lf %c", &dtemp, &dummy) == 1 ||
            sscanf(argv[arg], "p=%lf %c", &dtemp, &dummy) == 1) {
            p_black = dtemp;
        } else
        if (sscanf(argv[arg], "white=%lf %c", &dtemp, &dummy) == 1 ||
            sscanf(argv[arg], "p1=%lf %c", &dtemp, &dummy) == 1) {
            p_black = 1.0 - dtemp;
        } else
        if (sscanf(argv[arg], "d=%lf %c", &dtemp, &dummy) == 1 ||
            sscanf(argv[arg], "diag=%lf %c", &dtemp, &dummy) == 1 ||
            sscanf(argv[arg], "pd=%lf %c", &dtemp, &dummy) == 1 ||
            sscanf(argv[arg], "pdiag=%lf %c", &dtemp, &dummy) == 1) {
            p_diag = dtemp;
        } else
        if (sscanf(argv[arg], "d0=%lf %c", &dtemp, &dummy) == 1 ||
            sscanf(argv[arg], "db=%lf %c", &dtemp, &dummy) == 1 ||
            sscanf(argv[arg], "diag0=%lf %c", &dtemp, &dummy) == 1 ||
            sscanf(argv[arg], "diagblack=%lf %c", &dtemp, &dummy) == 1 ||
            sscanf(argv[arg], "diag_black=%lf %c", &dtemp, &dummy) == 1) {
            p_diag_black = dtemp;
        } else {
            fprintf(stderr, "%s: Unknown option.\n", argv[arg]);
            return EXIT_FAILURE;
        }

    switch (init_cluster(&c, rows, cols, p_black, p_diag, p_diag_black)) {
    case 0: break; /* OK */
    case ERR_INVALID:
        fprintf(stderr, "Invalid size.\n");
        return EXIT_FAILURE;
    case ERR_TOOLARGE:
        fprintf(stderr, "Size is too large.\n");
        return EXIT_FAILURE;
    case ERR_NOMEM:
        fprintf(stderr, "Not enough memory.\n");
        return EXIT_FAILURE;
    }

    if (!seed)
        seed = randomize(NULL);

    c.rng.state = seed;

    /* The largest possible cluster has n cells. */
    n = (size_t)rows * (size_t)cols;

    /* Print the comments describing the initial parameters. */
    printf("# seed: %" PRIu64 " (Xorshift 64*)\n", seed);
    printf("# size: %d rows, %d columns\n", rows, cols);
    printf("# P(black): %.6f (%" PRIu64 "/18446744073709551615)\n", p_black, c.p_black);
    printf("# P(connecting diagonally): %.6f (%" PRIu64 "/18446744073709551615)\n", p_diag, c.p_diag);
    printf("# P(black connecting diagonally): %.6f (%" PRIu64 "/18446744073709551615)\n", p_diag_black, c.p_diag_black);
    fflush(stdout);

    while (iters-->0)
        iterate(&c);

    printf("# Iterations: %" PRIu64 "\n", c.iterations);
    printf("#\n");
    printf("# size  white_clusters(size) black_clusters(size) clusters(size)\n");

    /* Note: c._histogram[0] == c._histogram[n] == 0, for ease of scanning. */
    for (i = 1; i <= n; i++)
        if (c.white_histogram[i-1] || c.white_histogram[i] || c.white_histogram[i+1] ||
            c.black_histogram[i-1] || c.black_histogram[i] || c.black_histogram[i+1])
            printf("%lu %" FMT_COUNT " %" FMT_COUNT " %" FMT_COUNT "\n",
                   (unsigned long)i,
                   c.white_histogram[i],
                   c.black_histogram[i],
                   c.white_histogram[i]+c.black_histogram[i]);

    /* Since we are exiting anyway, this is not really necessary. */
    free_cluster(&c);

    /* All done. */
    return EXIT_SUCCESS;
}
