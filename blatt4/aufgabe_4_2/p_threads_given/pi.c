#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include <pthread.h>

#define RNG_MOD 714025


static int pr_random(int* state)
{
    return (*state = (1366 * *state + 150889) % RNG_MOD);
}

static double pr_random_f(double range, int *state)
{
    return ((double) pr_random(state) / ((double) RNG_MOD)) * range;
}

static void* monto_carlo_thread_fn(void* arg)
{
    int samples = *((int*) arg);
    int hits, i, state = 0;
    double x, y;
    
    hits = 0;
    for (i = 0; i < samples; i++) {
        x = pr_random_f(1.0, &state);
        y = pr_random_f(1.0, &state);
        if (x * x + y * y < 1.0) {
            hits++;
        }
    }

    *((int*) arg) = hits;

    return NULL;
}

void 
get_difference(struct timespec *start, struct timespec *end, struct timespec *diff){

  if(start->tv_nsec > end->tv_nsec){
    diff->tv_nsec = (1000000000 + end->tv_nsec) - start->tv_nsec;
    diff->tv_sec = end->tv_sec - start->tv_sec - 1;
  }
  else {
    diff->tv_nsec = end->tv_nsec - start->tv_nsec;
    diff->tv_sec = end->tv_sec - start->tv_sec;
  }
}

void 
print_difference(struct timespec *diff){

  double seconds = (double) diff->tv_sec;
  double nanoseconds = (double) diff->tv_nsec / 1000000000.0;

  printf("time in seconds: %.2f\n", seconds + nanoseconds);
}

int main(int argc, char* argv[])
{

    struct timespec start, end, diff;

    //get start time
    if(clock_gettime(CLOCK_MONOTONIC, &start) != 0){
        fprintf(stderr, "Error: could not get clock time (clock_gettime).");
        exit(1);
    }
    int i, num_threads, num_samples, status, hits;
    int* thread_args;
    pthread_t* threads;
    double error, pi_calc;

    if (argc != 3) {
        fprintf(stderr, "usage: %s threads samples\n", argv[0]);
        return EXIT_FAILURE;
    }

    if ((num_threads = atoi(argv[1])) < 1) {
        fprintf(stderr, "invalid number of threads '%s'\n", argv[1]);
        return EXIT_FAILURE;
    }

    if ((num_samples = atoi(argv[2])) < 1) {
        fprintf(stderr, "invalid number of samples '%s'\n", argv[2]);
        return EXIT_FAILURE;
    }

    threads = calloc(num_threads, sizeof(*threads));
    if (!threads) {
        perror("calloc");
        return EXIT_FAILURE;
    }

    thread_args = calloc(num_threads, sizeof(*thread_args));
    if (!thread_args) {
        perror("calloc");
        return EXIT_FAILURE;
    }

    /* create threads */
    for (i = 0; i < num_threads; i++) {
        thread_args[i] = (num_samples / num_threads) + (num_samples % num_threads) * ((i+1) / num_threads);

        status = pthread_create(&threads[i], NULL, monto_carlo_thread_fn, &thread_args[i]);
        if (status) {
            fprintf(stderr, "failed to create pthread (code %d)", status);
            return EXIT_FAILURE;
        }
    }

    /* wait for threads */
    hits = 0;
    for (i = 0; i < num_threads; i++) {
        status = pthread_join(threads[i], NULL);
        if (status) {
            fprintf(stderr, "failed to join pthread (code %d)", status);
            return EXIT_FAILURE;
        }
        hits += thread_args[i];
    }

    pi_calc = 4.0 * hits / num_samples;
    error = (pi_calc - M_PI) / M_PI;

    printf("%d threads, x = %lf, pi = %lf, r_error = %le\n", num_threads, pi_calc, M_PI, error);
    
    free(threads);
    free(thread_args);

      //get end time
    if(clock_gettime(CLOCK_MONOTONIC, &end) != 0){
        fprintf(stderr, "Error: could not get clock time (clock_gettime).");
        exit(1);
    }

    get_difference(&start, &end, &diff);
    print_difference(&diff);

    return EXIT_SUCCESS;
}

