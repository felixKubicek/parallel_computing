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

void print_field(Line *buf, int lines, char* name);

/* --------------------- CA simulation -------------------------------- */


void print_line(Line line, int index)
{
    int z;
    for (z = 0; z < sizeof(Line); z++)
    {
        printf("%d", line[z]);
    }
    printf("   line %d \n", index);
}


/* random starting configuration */
static void initConfig(Line *buf, int lines, int rem_lines, int my_lines, int my_rank)
{
    int x, y;

    initRandomLEcuyer(424243);
    int no_rands = 0;

    /* calculate how often the random function was called before */
    if (my_rank <= rem_lines)
    {
        no_rands = my_rank * (lines + 1) * XSIZE;
    }
    else
    {
        no_rands = (rem_lines * (lines + 1) * XSIZE) +
                   (my_rank - rem_lines) * lines * XSIZE;
    }

    int i;
    long randResult = 0;
    for (i = 0; i < no_rands; i++)
    {
        randResult += (long) (randInt(100) >= 50);
    }

    for (y = 1;  y <= my_lines;  y++)
    {
        for (x = 1;  x <= XSIZE;  x++)
        {
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
    for (y = 0;  y <= lines + 1;  y++)
    {
        /* copy rightmost column to the buffer column 0 */
        buf[y][0      ] = buf[y][XSIZE];

        /* copy leftmost column to the buffer column XSIZE + 1 */
        buf[y][XSIZE + 1] = buf[y][1    ];
    }

}

/* make one simulation iteration with lines lines.
 * old configuration is in from, new one is written to to.
 */
static void simulate(Line *from, Line *to, int my_lines, int my_rank, int world_size)
{
    int x, y;
    boundary_left_right(from, my_lines);

    MPI_Request reqs[4];

    int top_neighbour_rank = (my_rank == 0) ? world_size - 1 : my_rank - 1;
    int bottom_neighbour_rank = (my_rank + 1) % world_size;

    /* send top line */
    MPI_Isend(from[1], sizeof(Line), MPI_CHAR, top_neighbour_rank, 0, MPI_COMM_WORLD, &reqs[0]);

    /* send bottom line */
    MPI_Isend(from[my_lines], sizeof(Line), MPI_CHAR, bottom_neighbour_rank, 1, MPI_COMM_WORLD, &reqs[1]);

    /* calculate inner field if present (when more than 2 my_lines) */
    if (my_lines > 2 )
    {
        int inner_field_start = 2;
        int inner_field_end = my_lines - 1;

        for (y = inner_field_start;  y <= inner_field_end;  y++)
        {
            for (x = 1;  x <= XSIZE;  x++)
            {
                to[y][x  ] = transition(from, x  , y);
            }
        }
    }

    /* receive ghost zones */

    /* receive top ghost zone */
    MPI_Irecv(from[0], sizeof(Line), MPI_CHAR, top_neighbour_rank, 1,  MPI_COMM_WORLD, &reqs[2]);
    /* receive bottom ghost zone */
    MPI_Irecv(from[my_lines + 1], sizeof(Line), MPI_CHAR, bottom_neighbour_rank, 0,  MPI_COMM_WORLD, &reqs[3]);

    MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);

    /* calculate outer field (bottom and top line with help of received ghost zones) */
    int top_line_index = 1;
    int bottom_line_index = my_lines;

    for (x = 1;  x <= XSIZE;  x++)
    {
        to[top_line_index][x  ] = transition(from, x  , top_line_index);
    }

    if (bottom_line_index != top_line_index)
    {
        for (x = 1;  x <= XSIZE;  x++)
        {
            to[bottom_line_index][x  ] = transition(from, x  , bottom_line_index);
        }

    }
}


/* --------------------- measurement ---------------------------------- */

int main(int argc, char **argv)
{
    int lines_global, its, i, *counts;
    Line *from, *to, *temp;
    Line *result;
    char *hash = NULL;

    assert(argc == 3);

    MPI_Init(&argc, &argv);

    lines_global = atoi(argv[1]);
    its = atoi(argv[2]);

    // get number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // get own rank
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // calculate number of lines per process
    int lines = lines_global / world_size;
    int rem_lines = lines_global % world_size;

    /*
    In order to have an equal distribution of lines, 'rem_lines' will be distributed
    among the first processes. For example if rem_lines = 3, then the 
    first three processes will get one more line than the others.
    */

    // gather different line counts for processes
    counts = calloc(world_size, sizeof(int));
    if (!counts)
    {
      printf("Error allocating requested memory.\n");
      exit(1);
    }

    for (i = 0; i < world_size; i++)
    {
        counts[i] = (rem_lines > i) ? lines + 1 : lines;    
    }
    int my_lines = counts[my_rank];

    
    if (my_lines <= 0)
    {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // create and initialize cellular automat fields
    from = calloc(my_lines + (2 * GHOSTZONE_SIZE), sizeof(Line));
    if (!from)
    {
      printf("Error allocating requested memory.\n");
      exit(1);
    }

    to   = calloc(my_lines + (2 * GHOSTZONE_SIZE), sizeof(Line));
    if (!to)
    {
      printf("Error allocating requested memory.\n");
      exit(1);
    }

    initConfig(from, lines, rem_lines, my_lines, my_rank);

    //simulate transition of cellular automat
    for (i = 0;  i < its;  i++)
    {
        simulate(from, to, my_lines, my_rank, world_size);
        
        temp = from;
        from = to;
        to = temp;
    }

    Line * start_buf = &from[1];

    int no_elements = counts[my_rank] * sizeof(Line);

    // send all results to process 0
    MPI_Send(*start_buf, no_elements, MPI_CHAR, 0, 1, MPI_COMM_WORLD);

    if (my_rank == 0) // process 0 collects all the results
    {
      result = calloc(lines_global,  sizeof(Line));
      
      MPI_Status status;
      for (i = 0; i < world_size; i++)
      {
        //check if a response has arrived
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        //calculate offset in result field depending on rank of sender process
        int offset = 0;
        if (status.MPI_SOURCE <= rem_lines)
        {
          offset = (status.MPI_SOURCE * (lines + 1));
        }
        else
        {
          offset = ((rem_lines * (lines + 1)) +
                   (status.MPI_SOURCE - rem_lines) * lines);
        }

        MPI_Recv(result[offset], counts[status.MPI_SOURCE] * sizeof(Line), MPI_CHAR, status.MPI_SOURCE,
                 status.MPI_TAG, MPI_COMM_WORLD, &status);
      }

      hash = getMD5DigestStr(result[0], sizeof(Line) * lines_global);
      printf("%s\n", hash);
      
      free(result);
    }

    free(hash);
    free(from);
    free(to);

    MPI_Finalize();

    return EXIT_SUCCESS;
}

void print_field(Line *buf, int lines, char* name){

  printf("field name: %s\n", name);
  int x, y;

  for (y = 0;  y <= lines + 1;  y++)
  {
    for (x = 0;  x <= XSIZE + 1;  x++)
    {
        printf("%d", buf[y][x]);
    }
    printf("   line %d\n", y);
  }
  printf("\n");
}
