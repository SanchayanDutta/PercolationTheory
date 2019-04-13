/* Courtesy: Nominal Animal <question@nominal-animal.net>
*/

#ifndef   PRNG_H
#define   PRNG_H
/*
Xorshift64* pseudo-random number generator.
*/
#include <inttypes.h>
#include <time.h>

#ifndef  static_inline
#define  static_inline  static inline
#endif

/* Generator state. */
typedef struct {
    uint64_t    state;
} prng;

/* Initialize generator.  Randomizes the state based on time. */
static_inline void  prng_init(prng *const  rng)
{
    int       churn = 128;
    uint64_t  state = UINT64_C(3069887672279) * (uint64_t)time(NULL)
                    ^ UINT64_C(60498839) * (uint64_t)clock();

    /* Zero state is invalid for Xorshift64*. */
    if (!state)
        state = 1;

    /* Churn the initial state, to reduce the time dependency. */
    while (churn-->0) {
        state ^= state >> 12;
        state ^= state << 25;
        state ^= state >> 27;
    }

    rng->state = state;
}

/* Char array large enough to hold the generator state in human-readable form. */
typedef  char  prng_seed_string[24];

/* Return the current generator state as a human-readable string,
   stored somewhere in the prng_seed_string. */
static_inline const char *prng_get_seed(prng *const  rng,
                                        char *const  to)
{
    uint64_t  state = rng->state;
    char     *s = to + sizeof (prng_seed_string);

    /* Save the state as a positive decimal integer. */
    *(--s) = '\0';
    do {
        *(--s) = '0' + (state % 10u);
        state /= 10u;
    } while (state);

    /* Done. */
    return (const char *)s;
}

/* Set the current generator based on a human-readable string.
   Returns a pointer to the first unparsed character, or NULL if an error occurs. */
static_inline const char *prng_set_seed(prng *const  rng,
                                        const char  *from)
{
    uint64_t  state;

    if (!from || !*from)
        return NULL;

    /* We'll ignore leading whitespace, if any. */
    while (*from == '\t' || *from == '\n' || *from == '\v' ||
           *from == '\f' || *from == '\r' || *from == ' ')
        from++;

    /* We do need at least one decimal digit. */
    if (!(*from >= '0' && *from <= '9'))
        return NULL;
    else
        state = *(from++) - '0';

    /* Parse the rest of the digits. */
    while (*from >= '0' && *from <= '9') {
        const uint64_t  oldstate = state;

        state = (state * 10) + (*(from++) - '0');

        /* Overflow? */
        if ((uint64_t)(state / 10) != oldstate)
            return NULL;
    }

    /* Zero state is not valid. */
    if (!state)
        return NULL;

    rng->state = state;

    /* Done. */
    return from;
}    

/* Probability limit type. */
typedef  uint64_t  prng_limit;

/* Define a probability limit used by prng_probability(). */
static_inline prng_limit  prng_set_probability(const double  p)
{
    if (p <= 0.0)
        return UINT64_C(0);
    else
    if (p <= 0.5)
        return 1 + (uint64_t)(18446744073709551615.0 * p);
    else
    if (p < 1.0)
        return UINT64_C(18446744073709551615)
             - (uint64_t)(18446744073709551615.0 * (double)(1 - p));
    else
        return UINT64_C(18446744073709551615);
}

/* Evaluate a probability. Returns 1 at the specified probability limit. */
static_inline int  prng_probability(prng *const       rng,
                                    const prng_limit  limit)
{
    uint64_t  state = rng->state;
    uint64_t  value;

    do {
        state ^= state >> 12;
        state ^= state << 25;
        state ^= state >> 27;
        value = state * UINT64_C(2685821657736338717);
    } while (!value);

    rng->state = state;

    return (value <= limit);
}

/* Return [0, 1] at uniform probablility. */
static_inline double prng_unit(prng *const rng)
{
    uint64_t  state = rng->state;
    uint64_t  value;

    do {
        state ^= state >> 12;
        state ^= state << 25;
        state ^= state >> 27;
        value = state * UINT64_C(2685821657736338717);
    } while (!value || value > UINT64_C(9223372036854775809));

    rng->state = state;

    return (double)(value - 1) / 9223372036854775808.0;
}

/* Return a value between min and max, inclusive, at uniform probability. */
static_inline double prng_drange(prng *const rng, const double min, const double max)
{
    uint64_t  state = rng->state;
    uint64_t  value;

    do {
        state ^= state >> 12;
        state ^= state << 25;
        state ^= state >> 27;
        value = state * UINT64_C(2685821657736338717);
    } while (!value || value > UINT64_C(9223372036854775809));

    rng->state = state;

    {
        const double  phase = (double)(value - 1) / (double)9223372036854775808.0;
        return (double)(phase * max) + (double)((1.0 - phase) * min);
    }
}

#endif /* PRNG_H */
