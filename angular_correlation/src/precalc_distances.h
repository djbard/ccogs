#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<unistd.h>

void define_ranges(float *vec_ranges_lo, float *vec_ranges_hi, float maxsep, int *vec_nbins, float *vec_binwidths)
{

    float r;
    for(int i=0;i<3;i++)
    {
        r = vec_ranges_hi[i] - vec_ranges_lo[i];
        vec_nbins[i] = (int)(r/maxsep);
        vec_binwidths[i] = r/vec_nbins[i];

    }

}

void assign_chunk_coordinate(float *coord, float *vec_ranges_lo, float *vec_binwidths, int *chunk_coord)
{

    int bin;
    for(int i=0;i<3;i++)
    {
        bin = (int)((coord[i]-vec_ranges_lo[i])/vec_binwidths[i]);
        chunk_coord[i] = bin;
        if(bin<0) {
            printf("LESS THAN 0!!!!");
            exit(-1);
        }
        //printf("%d\n",bin);
    }

}
