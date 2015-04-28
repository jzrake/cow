#include <limits.h>
#include <time.h>
#include "jsw_rand.h"

#define N JSW_STATE_SIZE /* 624 */
#define M 397
#define A 0x9908b0dfUL
#define U 0x80000000UL
#define L 0x7fffffffUL

/* Internal state */
//static unsigned long x[N];
//static int next; /* initialized to zero (ISO/IEC 9899) */

#define x R->x
#define next R->next

/* Initialize internal state */
void jsw_seed (jsw_rand_t *R, unsigned long s)
{
  int i;

  next = 0;
  x[0] = s & 0xffffffffUL;

  for ( i = 1; i < N; i++ ) {
    x[i] = ( 1812433253UL 
      * ( x[i - 1] ^ ( x[i - 1] >> 30 ) ) + i );
    x[i] &= 0xffffffffUL;
  }
}

/* Mersenne Twister */
unsigned long jsw_rand (jsw_rand_t *R)
{
  unsigned long y, a;
  int i;

  /* Refill x if exhausted */
  if ( next == N ) {
    next = 0;

    for ( i = 0; i < N - 1; i++ ) {
      y = ( x[i] & U ) | (x[i + 1] & L);
      a = ( y & 0x1UL ) ? A : 0x0UL;
      x[i] = x[( i + M ) % N] ^ ( y >> 1 ) ^ a;
    }

    y = ( x[N - 1] & U ) | (x[0] & L);
    a = ( y & 0x1UL ) ? A : 0x0UL;
    x[N - 1] = x[M - 1] ^ ( y >> 1 ) ^ a;
  }

  y = x[next++];

  /* Improve distribution */
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);

  return y;
}

/* Portable time seed */
unsigned jsw_time_seed(jsw_rand_t *R)
{
  time_t now = time ( 0 );
  unsigned char *p = (unsigned char *)&now;
  unsigned seed = 0;
  size_t i;

  for ( i = 0; i < sizeof now; i++ )
    seed = seed * ( UCHAR_MAX + 2U ) + p[i];

  return seed;
}




unsigned long jsw_random_int(jsw_rand_t *R, unsigned long a, unsigned long b)
{
  return a + jsw_rand(R) % (b - a);
}

double jsw_random_double(jsw_rand_t *R, double a, double b)
{
  return a + (double) jsw_rand(R) / UINT_MAX * (b - a);
}

float jsw_random_float(jsw_rand_t *R, float a, float b)
{
  return a + (float) jsw_rand(R) / UINT_MAX * (b - a);
}
