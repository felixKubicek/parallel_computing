#include <stdlib.h>
#include <stdio.h>
#include "random.h"
#include <pthread.h>
#include <math.h>

struct interimResult {
	int id;
	int circleCount;
	unsigned long long samples;
};

#define CIRCLERADIUS 1.0f

void *calculateSamples();
int isWithinCircle(double x, double y);

int main(int argc, char *argv[])
{

	if(argc != 3){
		printf("Wrong number of arguments! Usage: program [number_of_threads] [samples]\n");
		return 1;
	}

	// get number of threads and samples from program arguments
	char *p;
	int threadNum = atoi(argv[1]);
	unsigned long long samples = strtoull(argv[2], &p, 10);
	printf("number of threads: %d\n", threadNum);
	printf("number of samples: %d\n", samples);

	// thread ids
	pthread_t *threadIds = malloc(sizeof(pthread_t) * threadNum);

	// to save results from all threads
	struct interimResult *results = malloc(sizeof(struct interimResult) * threadNum);

	int samplesPerThread = samples/threadNum;
	int rest = samples % threadNum;


	int i;
	for(i = 0; i < threadNum; i++){		
     	//printf("In main: creating thread %d\n", i);

     	results[i] = (struct interimResult) {.id = i,
     			.circleCount = 0, 
     			.samples = samplesPerThread + rest};

     	if(i == threadNum -1){
     		results[i].samples += rest;
     	}
     	
     	pthread_create(&threadIds[i], NULL, calculateSamples, &results[i]);
    }

    int sum_circleCount = 0;
    for(i = 0; i < threadNum; i++){
    	pthread_join(threadIds[i], NULL );
    	sum_circleCount += results[i].circleCount;
    }

    double probability = (double) sum_circleCount / (double) samples;
	double pi = probability * 4;

	printf("PI: %f\n", pi);	
	printf("relative error: %f\n", (pi - M_PI)/M_PI);

	free(results);
	free(threadIds);

	return 0;
}

void *calculateSamples(void *s)
{
	struct interimResult *res = (struct interimResult *) s;
	//printf("id: %d samples: %llu\n", res->id, res->samples);
	
	unsigned long long i;
	for (i = 0; i < res->samples; i++){

		//generate coordinates for point
		double x = pr_random_f(1);
		double y = pr_random_f(1);		
		
		if(isWithinCircle(x,y)){
			res->circleCount += 1;
		}
	}
	pthread_exit(NULL);	
}

int isWithinCircle(double x, double y)
{
	return (x*x + y*y) <= CIRCLERADIUS;
}



