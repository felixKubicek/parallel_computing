#include "random.h"
#include "math.h"
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#define CIRCLE_RADIUS 1

static int circle_hits = 0;


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
  
  int k;
  #pragma omp parallel for schedule(static) private(k) reduction(+:circle_hits) num_threads(num_threads)
  for (k = 0; k < num_samples; k++)
  {  
    double x = pr_random_f(CIRCLE_RADIUS);
    double y = pr_random_f(CIRCLE_RADIUS);
 
    if (((x * x) + (y * y)) <= 1)
    {
      circle_hits++;
    } 
  }
  
  #pragma omp parallel num_threads(num_threads)
  {
    printf ("number of parallel threads: %d\n", omp_get_num_threads());
  }
  
  double pi = ((double) circle_hits / num_samples) * 4;
  printf ("estimation of pi: %f\n", pi);

  double relative_error = ((pi - M_PI) / M_PI);
  printf ("relative error: %f\n", relative_error);

  return 0;
}


