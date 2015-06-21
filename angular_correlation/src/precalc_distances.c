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
    return sqrt(diffx*diffx + diffy*diffy + diffz*diffz);
    //return diffx*diffx + diffy*diffy + diffz*diffz;
}

void read_in_6_cols(FILE* infile, int max_num, float *x, float *y, float *z, int *n, float *lo, float *hi, int max_to_read_in) {

    int num=0;
    float temp0, temp1, temp2, dummy, tempdum0, tempdum1, tempdum2;

    while(fscanf(infile, "%e %e %e %e %e %e", &tempdum0, &tempdum1, &tempdum2, &temp0, &temp1, &temp2) != EOF) {
        x[num] = temp0;
        y[num] = temp1;
        z[num] = temp2;

        // Keep track of ranges of values
        if(x[num]<lo[0])
            lo[0]=x[num]-1;
        if(x[num]>hi[0])
            hi[0]=x[num]+1;

        if(y[num]<lo[1])
            lo[1]=y[num]-1;
        if(y[num]>hi[1])
            hi[1]=y[num]+1;

        if(z[num]<lo[2])
            lo[2]=z[num]-1;
        if(z[num]>hi[2])
            hi[2]=z[num]+1;

        if(num>=max_num) {
            printf("Exceeded max num galaxies: %d", max_num);
            exit(-1);
        }

        // Diagnostics
        if (num<10) {
            printf("%f %f %f\n", x[num],y[num],z[num]);
        }

        num += 1;

        if (num>=max_to_read_in)
            break;
    }

    *n = num;
}

