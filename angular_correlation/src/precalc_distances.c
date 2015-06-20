#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<unistd.h>

void define_ranges(float *vec_ranges_lo, float *vec_ranges_hi, float maxsep, int *vec_nbins, float *vec_binwidths)
{

    float r;
    int i = 0;
    for(i=0;i<3;i++)
    {
        r = vec_ranges_hi[i] - vec_ranges_lo[i];
        vec_nbins[i] = (int)(r/maxsep);
        vec_binwidths[i] = r/vec_nbins[i];

    }

}

void assign_chunk_coordinate(float *coord, float *vec_ranges_lo, float *vec_binwidths, int *chunk_coord)
{

    int bin;
    int i;
    for(i=0;i<3;i++)
    {
        bin = (int)((coord[i]-vec_ranges_lo[i])/vec_binwidths[i]);
        chunk_coord[i] = bin;
        if(bin<0) {
            printf("LESS THAN 0!!!!\n");
            printf("%d %f %f %f\n",i,coord[i],vec_ranges_lo[i],vec_binwidths[i]);
            exit(-1);
        }
        //printf("%d\n",bin);
    }

}

float distance(float x0, float y0, float z0, float x1, float y1, float z1)
{

    float diffx = x0-x1;
    float diffy = y0-y1;
    float diffz = z0-z1;

    //printf("diff: %f %f %f\n",diffx,x0,x1);
    //printf("diff: %f %f %f\n",diffy,y0,y1);
    //printf("diff: %f %f %f\n",diffz,z0,z1);
    //return sqrt(diffx*diffx + diffy*diffy + diffz*diffz);
    return diffx*diffx + diffy*diffy + diffz*diffz;
}
