#include "random.h"
#include "math.h"
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>

#define THREAD_ARG 0
#define CIRCLE_RADIUS 1

static int circle_hits = 0;


typedef struct
{
  int samples_to_compute;
  int circle_hits_ret;
} thread_arg;


int
circle_hit(double x, double y)
{
  return ((x * x) + (y * y) <= 1);
}

static void *
thread_routine (void * arg)
{
  thread_arg *local_arg = (thread_arg *) arg;
  int samples_to_compute = local_arg->samples_to_compute;
  int local_circle_hits = 0;
   
  int k;
  for (k = 0; k < samples_to_compute; k++)
  {  
    double x = pr_random_f(CIRCLE_RADIUS);
    double y = pr_random_f(CIRCLE_RADIUS);
 
    if (circle_hit(x,y))
    {
      local_circle_hits ++;
    }
  }
  
  printf ("thread_id: %d; num_samples: %d; circle_hits: %d\n", (int) pthread_self (), samples_to_compute, local_circle_hits); 
  
  local_arg->circle_hits_ret = local_circle_hits;  

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
  int remaining_samples = num_samples % num_threads;

  pthread_t *tid;
  thread_arg *thread_args; 
  
  if ((tid = (pthread_t *) malloc(sizeof(pthread_t) * num_threads)) == NULL)
  {
    fprintf (stderr, "cannot allocate enough memory for pthread array\n");
    return 1;
  }

  if ((thread_args = (thread_arg *) malloc(sizeof(thread_arg) * num_threads)) == NULL)
  {
    fprintf (stderr, "cannot allocate enough memory for pthread arguments array\n");
    return 1;
  }
  
  int i;
  for (i = 0; i < num_threads; i++)
  {
    thread_arg *current_arg = &thread_args[i];
    current_arg->samples_to_compute = (remaining_samples > i) ? samples_per_thread + 1 : samples_per_thread;
   
    pthread_create (&tid[i], NULL, &thread_routine, (void *) current_arg);
  }

  for (i = 0; i < num_threads; i++)
  {
    pthread_join (tid[i], NULL);
    circle_hits += thread_args[i].circle_hits_ret;
  }
  
  printf ("main() reporting that all %d threads have terminated\n", num_threads);
  printf ("global number of circle hits: %d\n", circle_hits);
  
  double pi = ((double) circle_hits / num_samples) * 4;
  printf ("estimation of pi: %f\n", pi);

  double relative_error = ((pi - M_PI) / M_PI);
  printf ("relative error: %f\n", relative_error);

  // free heap 
  free(tid);
  free(thread_args);

  return 0;
}


