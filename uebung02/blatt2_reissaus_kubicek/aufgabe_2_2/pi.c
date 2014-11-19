#include "random.h"
#include "math.h"
#include <stdio.h>
#include <pthread.h>

#define THREAD_ARG 0
#define CIRCLE_RADIUS 1

static pthread_mutex_t circle_misses_mutex = PTHREAD_MUTEX_INITIALIZER;
static int circle_misses = 0;


int
circle_miss(double x, double y)
{
  return ((x * x) + (y * y) > 1);
}

static void *
thread_routine (void * arg)
{
  int *samples_to_compute = arg;
   
  //printf ("thread_id: %d; num_samples: %d; x: %f; y: %f\n", (int) pthread_self (), *samples_to_compute, x, y); 
  int k;
  for (k = 0; k < *samples_to_compute; k++)
  {  
    double x = pr_random_f(CIRCLE_RADIUS);
    double y = pr_random_f(CIRCLE_RADIUS);
 
    if (circle_miss(x,y))
    {
      pthread_mutex_lock(&circle_misses_mutex);
      //printf("begin circle_mutex: %d\n",(int) pthread_self ()); 
      circle_misses ++;
      //printf("end   circle_mutex: %d\n",(int) pthread_self ());
      pthread_mutex_unlock(&circle_misses_mutex);
    }
  }
  
  pthread_exit(NULL);
}

int 
main (int argc, char** argv)
{
  if (argc != 3)
   {
      printf ("Usage: pi <number_threads> <number_samples>\n");
      return 1;
   }
   
  unsigned int num_threads;
  unsigned int num_samples;
  
  if (sscanf (argv[1],"%u",&num_threads) != 1)
  {
    fprintf (stderr, "<number_threads> has to be a positive integer\n");
    return 1;
  }
  
  if (sscanf (argv[2],"%u",&num_samples) != 1)
  {
    fprintf (stderr, "<number_samples> has to be a positive integer\n");
    return 1;
  }
   
  printf ("start number_threads: %d ,number_samples: %d\n", num_threads, num_samples);
  
  // calculating number of samples per thread
  int samples_per_thread = num_samples / num_threads;
  int samples_per_thread_plus_one = samples_per_thread + 1;
  int remaining_samples = num_samples % num_threads;

  pthread_t tid[num_threads];
  
  int i;
  for (i = 0; i < num_threads; i++)
  {
    int *samples_to_compute = (remaining_samples > i) ? &samples_per_thread_plus_one : &samples_per_thread;
   
    pthread_create (&tid[i], NULL, &thread_routine, (void *) samples_to_compute);
  }

  for (i = 0; i < num_threads; i++)
  {
    pthread_join (tid[i], NULL);
  }
  
  printf ("main() reporting that all %d threads have terminated\n", num_threads);
  
  int circle_hits = num_samples - circle_misses;
  double pi = ((double) circle_hits / num_samples) * 4;
  printf ("estimation of pi: %f\n", pi);

  double relative_error = ((pi - M_PI) / M_PI);
  printf ("relative error: %f\n", relative_error);

  return 0;
}


