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

typedef char charState;
typedef charState charLine[XSIZE + 2];

typedef int State;
typedef State Line[32 + 1];


/* determine random integer between 0 and n-1 */
#define randInt(n) ((int)(nextRandomLEcuyer() * n))

void print_field(Line *buf, int lines, char* name);
void print_line(Line * buf, int index);
void print_char_line(charLine * buf, int index);

/* --------------------- CA simulation -------------------------------- */

void setBit(State * line, int i)
{
  line[i/32] |= (1 << (i % 32));
}

void clearBit(State * line, int i)
{
  line[i/32] &= (~ (1 << (i % 32)));
}

int testBit(State * line, int i)
{
  return ((line[i/32] & (1 << (i % 32))) != 0);
}

void setBitTo(State * line, int i, int value)
{
  if (value > 0)
  {
    setBit(line, i);
  }
  else
  {
    clearBit(line,i);
  }
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
            setBitTo(buf[y], x, randInt(100) >= 50);
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
    (anneal[testBit((a)[(y)-1],(x)-1) + testBit((a)[(y)],(x)-1) + testBit((a)[(y)+1],(x)-1) +\
            testBit((a)[(y)-1],(x)) + testBit((a)[(y)],(x)) + testBit((a)[(y)+1],(x)) +\
            testBit((a)[(y)-1],(x)+1) + testBit((a)[(y)],(x)+1) + testBit((a)[(y)+1],(x)+1)])


/* treat torus like boundary conditions for left and right side */
static void boundary_left_right(Line *buf, int lines)
{
    int y;
    for (y = 0;  y <= lines + 1;  y++)
    {
        /* copy rightmost column to the buffer column 0 */
        setBitTo(buf[y], 0, testBit(buf[y], XSIZE));

        /* copy leftmost column to the buffer column XSIZE + 1 */
        setBitTo(buf[y], XSIZE + 1, testBit(buf[y], 1));
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

    /* receive ghost zones */
    /* receive top ghost zone */
    MPI_Irecv(from[0], sizeof(Line)/sizeof(int), MPI_INT, top_neighbour_rank, 1,  MPI_COMM_WORLD, &reqs[2]);
    /* receive bottom ghost zone */
    MPI_Irecv(from[my_lines + 1], sizeof(Line)/sizeof(int), MPI_INT, bottom_neighbour_rank, 0,  MPI_COMM_WORLD, &reqs[3]);

    /* send top line */
    MPI_Isend(from[1], sizeof(Line)/sizeof(int), MPI_INT, top_neighbour_rank, 0, MPI_COMM_WORLD, &reqs[0]);
    /* send bottom line */
    MPI_Isend(from[my_lines], sizeof(Line)/sizeof(int), MPI_INT, bottom_neighbour_rank, 1, MPI_COMM_WORLD, &reqs[1]);

    /* calculate inner field if present (when more than 2 my_lines) */
    if (my_lines > 2 )
    {
        int inner_field_start = 2;
        int inner_field_end = my_lines - 1;

        for (y = inner_field_start;  y <= inner_field_end;  y++)
        {
            for (x = 1;  x <= XSIZE;  x++)
            {
                setBitTo(to[y], x, transition(from, x  , y));
            }
        }
    }

    MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);

    /* calculate outer field (bottom and top line with help of received ghost zones) */
    int top_line_index = 1;
    int bottom_line_index = my_lines;

    for (x = 1;  x <= XSIZE;  x++)
    {
        setBitTo(to[top_line_index], x, transition(from, x  , top_line_index));
    }

    if (bottom_line_index != top_line_index)
    {
        for (x = 1;  x <= XSIZE;  x++)
        {
            setBitTo(to[bottom_line_index], x, transition(from, x  , bottom_line_index));
        }

    }
}


/* --------------------- measurement ---------------------------------- */

int main(int argc, char **argv)
{
    int lines_global, its, i, *line_counts, *line_displ;
    Line *from, *to, *temp;
    Line *result;
    
    charLine *char_result;
    
    char *hash = NULL;

    assert(argc == 3);

    MPI_Init(&argc, &argv);

    lines_global = atoi(argv[1]);
    its = atoi(argv[2]);
    
    //int print_line_index;
    //print_line_index = atoi(argv[3]);
    
    
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

    // gather different line line_counts for processes
    line_counts = calloc(world_size, sizeof(int));
    if (!line_counts)
    {
      printf("Error allocating requested memory.\n");
      exit(1);
    }
    
    // line_displacement for gather
    line_displ = calloc(world_size, sizeof(int));
    
    if (!line_displ)
    {
      printf("Error allocating requested memory.\n");
      exit(1);
    }
    
    for (i = 0; i < world_size; i++)
    {
        
      // set values of line_counts array  
      line_counts[i] = (rem_lines > i) ? lines + 1 : lines;
      
      // set values of line_displacement array  
      if (i <= rem_lines)
      {
        line_displ[i] = (i * (lines + 1));
      }
      else
      {
        line_displ[i] = ((rem_lines * (lines + 1)) + (i - rem_lines) * lines);
      }             
    }
    
    int my_lines = line_counts[my_rank];

    
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
    
    Line *start_buf = &from[1];
    
    if (my_rank == 0) // process 0 collects all the results
    {
      result = calloc(lines_global,  sizeof(Line));
      char_result = calloc(lines_global,  XSIZE + 2);
    }
    
    MPI_Datatype mpi_line_type;
  
    MPI_Type_contiguous(sizeof(Line)/sizeof(int), MPI_INT, &mpi_line_type);    
    MPI_Type_commit(&mpi_line_type);
     
    MPI_Gatherv(*start_buf, my_lines, mpi_line_type, result, line_counts, line_displ, mpi_line_type, 0, MPI_COMM_WORLD);
    
    
    if (my_rank == 0)
    {
      //print_line(result, print_line_index);
      //print_char_line(char_result, print_line_index);
      
      int y;
      for(y=0; y < lines_global; y++)
      {
        int x;
        for (x = 0; x < XSIZE + 2; x++)
        {
          char_result[y][x] = (char) testBit(result[y], x);
        }
      }
      
      hash = getMD5DigestStr(char_result[0], sizeof(charLine) * lines_global);
      printf("hash: %s\n", hash);
      
      // clean up
      free(char_result);
      free(result); 
      free(hash);
    }
    
    // clean up   
    free(line_counts);  
    free(line_displ);  
    free(from);
    free(to);
    
    
    
    
    MPI_Type_free(&mpi_line_type);
    MPI_Finalize();

    return EXIT_SUCCESS;
}

void print_line(Line * buf, int index)
{
    int z;
    int sum = 0;
    
    for (z = 0; z < (XSIZE + 2); z++)
    {
        int bit = testBit(buf[index], z);
        sum += bit;
      
        printf("%d", bit);
    }
    printf("   line %d \n", index);
    printf("sum of line %d: %d\n\n", index, sum);
}

void print_char_line(charLine * buf, int index)
{
    int z;
    int sum = 0;
    
    for (z = 0; z < (XSIZE + 2); z++)
    {
        int bit = buf[index][z];
        sum += bit;
      
        printf("%d", bit);
    }
    printf("   line %d \n", index);
    printf("sum of line %d: %d\n\n", index, sum);
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
