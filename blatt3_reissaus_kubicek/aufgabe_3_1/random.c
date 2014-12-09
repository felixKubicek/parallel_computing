#include <pthread.h>
#include <stdio.h>
#include "random.h"

#define RNG_MOD 714025


int pr_random(void)
{
  static int state = 0;
  return (state = (1366 * state + 150889) % RNG_MOD);

}

double pr_random_f(double range)
{  
	int pr_random_safe = 0;
  #pragma omp critical(random)
	{
  	pr_random_safe = pr_random();
  }
  return ((double) pr_random_safe / (double) RNG_MOD) * range;
}
