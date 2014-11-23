#include <pthread.h>

#define RNG_MOD 714025

pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;

int pr_random(void)
{
    pthread_mutex_lock(&mut);
    static int state = 0;
    int rand = (state = (1366 * state + 150889) % RNG_MOD);
    pthread_mutex_unlock(&mut);
    return rand;
}

double pr_random_f(double range)
{
    return ((double) pr_random() / (double) RNG_MOD) * range;
}
