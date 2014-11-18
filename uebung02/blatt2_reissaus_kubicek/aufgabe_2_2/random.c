#include <pthread.h>
#include <stdio.h>
#include "random.h"

#define RNG_MOD 714025


static pthread_mutex_t get_random_mutex = PTHREAD_MUTEX_INITIALIZER;


int pr_random(void)
{
  static int state = 0;
    
  return (state = (1366 * state + 150889) % RNG_MOD);
}

double pr_random_f(double range)
{  
  pthread_mutex_lock(&get_random_mutex);
  //printf("begin random_mutex: %d\n",(int) pthread_self ());
  int pr_random_safe = pr_random();
  //printf("end   random_mutex: %d\n",(int) pthread_self ());
  pthread_mutex_unlock(&get_random_mutex);
  
  return ((double) pr_random_safe / (double) RNG_MOD) * range;
}
