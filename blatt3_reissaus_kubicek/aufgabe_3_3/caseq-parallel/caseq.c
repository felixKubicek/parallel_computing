/* (c) 1996,1997 Peter Sanders, Ingo Boesnach */
/* simulate a cellular automaton (serial version)
 * periodic boundaries
 * 
 * #1: Number of lines 
 * #2: Number of iterations to be simulated
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <assert.h>

#include "random.h"
#include "md5tool.h"
#include <mpi.h>

/* size of ghostzone (one line for upper and lower region each) */
#define GHOSTZONE_SIZE 1

/* horizontal size of the configuration */
#define XSIZE 1024

/* "ADT" State and line of states (plus border) */
typedef char State;
typedef State Line[XSIZE + 2];

/* determine random integer between 0 and n-1 */
#define randInt(n) ((int)(nextRandomLEcuyer() * n))

/* --------------------- CA simulation -------------------------------- */


void print_line(Line line)
{
  int z;
  for (z=0; z < sizeof(Line); z++)
  {
    printf("%d", line[z]);
  }
  printf("\n");  
}


/* random starting configuration */
static void initConfig(Line *buf, int lines, int my_lines, int my_rank)
{ int x, y;
  
  initRandomLEcuyer(424243);

  int i; 
  long randResult = 0;
  for(i = 0; i< (my_rank * lines * XSIZE); i++)
  {
    randResult += (long) (randInt(100)>=50);
  }   
  
  if (my_rank == 2){
    printf("randResult: %ld\n", randResult);
  }
  
 
  for (y = 1;  y <= my_lines;  y++) {
    for (x = 1;  x <= XSIZE;  x++) {
       buf[y][x] = randInt(100) >= 50;
    }
  }
}

/* annealing rule from ChoDro96 page 34 
 * the table is used to map the number of nonzero
 * states in the neighborhood to the new state
 */
static State anneal[10] = {0, 0, 0, 0, 1, 0, 1, 1, 1, 1};

/* a: pointer to array; x,y: coordinates; result: n-th element of anneal,
      where n is the number of neighbors */
#define transition(a, x, y) \
   (anneal[(a)[(y)-1][(x)-1] + (a)[(y)][(x)-1] + (a)[(y)+1][(x)-1] +\
           (a)[(y)-1][(x)  ] + (a)[(y)][(x)  ] + (a)[(y)+1][(x)  ] +\
           (a)[(y)-1][(x)+1] + (a)[(y)][(x)+1] + (a)[(y)+1][(x)+1]])


/* treat torus like boundary conditions for left and right side */
static void boundary_left_right(Line *buf, int lines)
{  
  int y;
  for (y = 1;  y <= lines;  y++) {
      /* copy rightmost column to the buffer column 0 */
      buf[y][0      ] = buf[y][XSIZE];

      /* copy leftmost column to the buffer column XSIZE + 1 */
      buf[y][XSIZE+1] = buf[y][1    ];
   }

}



/* make one simulation iteration with lines lines.
 * old configuration is in from, new one is written to to.
 */
static void simulate(Line *from, Line *to, int lines, int my_rank, int world_size)
{
   int x,y;
   boundary_left_right(from, lines);
   
   char *hash1;
   char *hash2;
   char *hash3;
   char *hash4;
   
   MPI_Request reqs[4];
   
          
   int top_neighbour_rank = (my_rank == 0) ? world_size - 1 : my_rank - 1;
   int bottom_neighbour_rank = (my_rank + 1) % world_size;
   
   /* send top line */
   MPI_Isend(from[1], sizeof(Line), MPI_CHAR, top_neighbour_rank, 0, MPI_COMM_WORLD, &reqs[0]);   
   /* send bottom line */
   MPI_Isend(from[lines], sizeof(Line), MPI_CHAR, bottom_neighbour_rank, 1, MPI_COMM_WORLD, &reqs[1]);
   
   
   hash1 = getMD5DigestStr(from[1], sizeof(Line));
   //printf("sent top line: %s my_rank: %d \n", hash1, my_rank);
   
   hash2 = getMD5DigestStr(from[lines], sizeof(Line));
   //printf("sent bottom line: %s my_rank: %d \n", hash2, my_rank);
   
   
   /* calculate inner field if present (when more than 2 lines) */
   if (lines > 2 )
   {

     int inner_field_start = 2;
     int inner_field_end = lines - 1;

    
     for (y = inner_field_start;  y <= inner_field_end;  y++) {
        for (x = 1;  x <= XSIZE;  x++) {
           to[y][x  ] = transition(from, x  , y);
         }   
      }
   }
  
  
   /* receive ghost zones */
     
   /* receive top ghost zone */
   MPI_Irecv(from[0], sizeof(Line), MPI_CHAR, top_neighbour_rank, 1,  MPI_COMM_WORLD, &reqs[2]);
   /* receive bottom ghost zone */
   MPI_Irecv(from[lines+1], sizeof(Line), MPI_CHAR, bottom_neighbour_rank, 0,  MPI_COMM_WORLD, &reqs[3]);
   
   MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);
   
   
   
   
   
    
   hash3 = getMD5DigestStr(from[0], sizeof(Line));
   //printf("received top ghost zone: %s my_rank: %d\n", hash3, my_rank);
   
   hash4 = getMD5DigestStr(from[lines+1], sizeof(Line));
   //printf("received bottom ghost zone: %s my_rank: %d\n", hash4, my_rank);
   
   
   
   /* calculate outer field (bottom and top line with help of received ghost zones) */
   
   int top_line_index = 1;
   int bottom_line_index = lines;
   
   for (x = 1;  x <= XSIZE;  x++) {
      to[top_line_index][x  ] = transition(from, x  ,top_line_index);
    }

   if (bottom_line_index != top_line_index)
   {
     for (x = 1;  x <= XSIZE;  x++) {
        to[bottom_line_index][x  ] = transition(from, x  ,bottom_line_index);
      }
     
   } 
   
   /* free hash resources */
   free(hash1);
   free(hash2);
   free(hash3);
   free(hash4);
   
}


/* --------------------- measurement ---------------------------------- */

int main(int argc, char** argv)
{  
   int lines_global, its;
   int i;
   Line *from, *to, *temp, *result;
   char* hash = NULL;

   assert(argc == 3);
   
   
   MPI_Init(NULL, NULL);
   
   lines_global = atoi(argv[1]);
   its = atoi(argv[2]);
   
   // get number of processes
   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);

   // get own rank
   int my_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   
   
   
   
   int lines = lines_global / world_size;
   int rem_lines = lines_global % world_size;
   
   int my_lines = lines;
   
   if (my_lines <= 0)
   {
     MPI_Abort(MPI_COMM_WORLD, 1); 
   }
   else if (my_rank == world_size - 1)
   {
     my_lines += rem_lines; 
   }
   
   //printf("Start simulation: process %d; world_size %d; number_lines: %d \n", my_rank, world_size, my_lines);
  
   from = calloc(my_lines + (2 * GHOSTZONE_SIZE), sizeof(Line));
   to   = calloc(my_lines + (2 * GHOSTZONE_SIZE), sizeof(Line));

   initConfig(from, lines, my_lines, my_rank);
   
   
   
   
   

   for (i = 0;  i < its;  i++) {
      simulate(from, to, my_lines, my_rank, world_size);
    
      temp = from;  
	    from = to;  
	    to = temp;
      
   }
   
   if (my_rank == 0){ // gathers all results
    result = calloc(lines_global + 2, sizeof(Line));  
   } 
   
   // MPI_Gather(void* send_data, int send_count, MPI_Datatype send_datatype,
   //            void* recv_data, int recv_count, MPI_Datatype recv_datatype,
   //            int root, MPI_Comm communicator)

   //int size_my_area = my_lines * sizeof(Line);
   //MPI_Gather(to[1], size_my_area, MPI_CHAR, result + 1 + (my_rank * size_my_area), size_my_area, MPI_CHAR, 0,
   //   MPI_COMM_WORLD);
  
   if(my_rank == 0){
     hash3 = getMD5DigestStr(to[0], sizeof(Line));
     //printf("received top ghost zone: %s my_rank: %d\n", hash3, my_rank);
   }
   // free resources 
   
   free(hash);
   free(from);
   free(to);
     
   MPI_Finalize();

   return EXIT_SUCCESS;
}
