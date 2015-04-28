#ifndef JSW_RAND_H
#define JSW_RAND_H

#define JSW_STATE_SIZE 624
typedef struct jsw_rand_t jsw_rand_t;

struct jsw_rand_t
{
  unsigned long x[JSW_STATE_SIZE];
  int next;
} ;

/* Seed the RNG. Must be called first */
void          jsw_seed(jsw_rand_t *R, unsigned long s);

/* Return a 32-bit random number */
unsigned long jsw_rand(jsw_rand_t *R);

/* Seed with current system time */
unsigned      jsw_time_seed(jsw_rand_t *R);


unsigned long jsw_random_int(jsw_rand_t *R, unsigned long a, unsigned long b);
double jsw_random_double(jsw_rand_t *R, double a, double b);
float jsw_random_float(jsw_rand_t *R, float a, float b);


#endif
