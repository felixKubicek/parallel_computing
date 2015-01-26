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


/* horizontal size of the configuration */
#define XSIZE 1024

/* "ADT" State and line of states (plus border) */
typedef char State;
typedef State Line[XSIZE + 2];

/* determine random integer between 0 and n-1 */
#define randInt(n) ((int)(nextRandomLEcuyer() * n))
void print_field(Line *buf, int lines, char* name);
void print_line(Line *buf, int index);

/* --------------------- CA simulation -------------------------------- */

/* random starting configuration */
static void initConfig(Line *buf, int lines)
{
    int x, y;
    initRandomLEcuyer(424243);
    for (y = 1;  y <= lines;  y++)
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

/* treat torus like boundary conditions */
static void boundary(Line *buf, int lines)
{
    int x, y;

    for (y = 0;  y <= lines + 1;  y++)
    {
        /* copy rightmost column to the buffer column 0 */
        buf[y][0      ] = buf[y][XSIZE];

        /* copy leftmost column to the buffer column XSIZE + 1 */
        buf[y][XSIZE + 1] = buf[y][1    ];
    }

    for (x = 0;  x <= XSIZE + 1;  x++)
    {
        /* copy bottommost row to buffer row 0 */
        buf[0][x      ] = buf[lines][x];

        /* copy topmost row to buffer row lines + 1 */
        buf[lines + 1][x] = buf[1][x    ];
    }
}

/* make one simulation iteration with lines lines.
 * old configuration is in from, new one is written to to.
 */
static void simulate(Line *from, Line *to, int lines)
{
    int x, y;

    boundary(from, lines);

    for (y = 1;  y <= lines;  y++)
    {
        for (x = 1;  x <= XSIZE;  x++)
        {
            to[y][x  ] = transition(from, x  , y);
        }
    }
}


/* --------------------- measurement ---------------------------------- */

int main(int argc, char **argv)
{
    int lines, its;
    int i;
    Line *from, *to, *temp;
    char *hash;

    assert(argc == 3);

    lines = atoi(argv[1]);
    its   = atoi(argv[2]);
    
    //int print_line_index;
    //print_line_index = atoi(argv[3]);

    from = (Line*) calloc((lines + 2), sizeof(Line));
    if (!from)
    {
      printf("Error allocating requested memory.\n");
      exit(1);
    }

    to = (Line*) calloc((lines + 2), sizeof(Line));
    if (!to)
    {
      printf("Error allocating requested memory.\n");
      exit(1);
    }

    initConfig(from, lines);
    
    for (i = 0;  i < its;  i++)
    {
        simulate(from, to, lines);

        temp = from;
        from = to;
        to = temp;


    }

    //print_line(from, print_line_index);
    
    hash = getMD5DigestStr(from[1], sizeof(Line) * (lines));
    printf("hash: %s\n", hash);
    
    free(from);
    free(to);
    free(hash);

    return EXIT_SUCCESS;
}

void print_line(Line * buf, int index)
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





