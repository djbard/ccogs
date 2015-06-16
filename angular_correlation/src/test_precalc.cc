//#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<unistd.h>

#include "precalc_distances.h"

#define NGALS 1000000
#define MAX_NGALS_IN_SUBVOLUME 600
#define NBINS_SUBVOLUME_X 15
#define NBINS_SUBVOLUME_Y 25
#define NBINS_SUBVOLUME_Z 12

int main()
{

    float maxsep = 150.;
    float vec_ranges_lo[3] = {-2000., -1720., 755.};
    float vec_ranges_hi[3] = {140., 1780., 2500.};
    float vec_ranges_width[3] = {0.0,0.0,0.0};
    int vec_nbins[3] = {0,0,0};
    float vec_binwidths[3] = {0.0,0.0,0.0};

    ////////////////////////////////////////////////////////////////////////////
    // Define the ranges for our galaxy chunks
    ////////////////////////////////////////////////////////////////////////////
    define_ranges(vec_ranges_lo, vec_ranges_hi, maxsep, vec_nbins, vec_binwidths);

    for(int i=0;i<3;i++)
    {
        printf("%d %d %f\n",i,vec_nbins[i],vec_binwidths[i]);
        vec_ranges_width[i] = vec_ranges_hi[i]-vec_ranges_lo[i];
    }

    //fflush(stdout); 
    //exit(1);

    ////////////////////////////////////////////////////////////////////////////
    // Gen some fake galaxies.
    ////////////////////////////////////////////////////////////////////////////
    //float galsx[NGALS];
    //float galsy[NGALS];
    //float galsz[NGALS];
    float *galsx, *galsy, *galsz;
    galsx = (float*)malloc(NGALS*sizeof(float));
    galsy = (float*)malloc(NGALS*sizeof(float));
    galsz = (float*)malloc(NGALS*sizeof(float));


    for(int i=0;i<NGALS;i++)
    {
        galsx[i] = vec_ranges_width[0]*rand()/(float)RAND_MAX + vec_ranges_lo[0];
        galsy[i] = vec_ranges_width[1]*rand()/(float)RAND_MAX + vec_ranges_lo[1];
        galsz[i] = vec_ranges_width[2]*rand()/(float)RAND_MAX + vec_ranges_lo[2];

    }

    //exit(1);

    ////////////////////////////////////////////////////////////////////////////
    // Assign a chunk coordinate to each
    ////////////////////////////////////////////////////////////////////////////
    //int gal_chunk_coordinate[NGALS][3];
    int *gal_chunk_coordinate[NGALS];
    for(int i=0;i<NGALS;i++) {
        gal_chunk_coordinate[i] = (int*)malloc(NGALS*sizeof(int));
    }

    float coord[3] = {0.0, 0.0, 0.0};
    int chunk_coord[3] = {0,0,0};
    for(int i=0;i<NGALS;i++) {
        coord[0] = galsx[i];
        coord[1] = galsy[i];
        coord[2] = galsz[i];
        //printf("coord: %f %f %f\n",coord[0],coord[1],coord[2]);
        assign_chunk_coordinate(coord, vec_ranges_lo, vec_binwidths, chunk_coord);
        for(int j=0;j<3;j++) {
            //printf("%d ",chunk_coord[j]);
            gal_chunk_coordinate[i][j] = chunk_coord[j];
        }
        //printf("\n");
        //printf("%f %f %f  (%d %d %d)\n",galsx[i],galsy[i],galsz[i], gal_chunk_coordinate[i][0],gal_chunk_coordinate[i][1],gal_chunk_coordinate[i][2]);
    }

    //exit(-1);

    ////////////////////////////////////////////////////////////////////////////
    // Divide up the galaxies into the chunked voxels.
    ////////////////////////////////////////////////////////////////////////////
    //float gal_chunksx[NBINS_SUBVOLUME_X][NBINS_SUBVOLUME_Y][NBINS_SUBVOLUME_Z][MAX_NGALS_IN_SUBVOLUME];
    //float gal_chunksy[NBINS_SUBVOLUME_X][NBINS_SUBVOLUME_Y][NBINS_SUBVOLUME_Z][MAX_NGALS_IN_SUBVOLUME];
    //float gal_chunksz[NBINS_SUBVOLUME_X][NBINS_SUBVOLUME_Y][NBINS_SUBVOLUME_Z][MAX_NGALS_IN_SUBVOLUME];

    float *gal_chunksx[NBINS_SUBVOLUME_X][NBINS_SUBVOLUME_Y][NBINS_SUBVOLUME_Z];
    float *gal_chunksy[NBINS_SUBVOLUME_X][NBINS_SUBVOLUME_Y][NBINS_SUBVOLUME_Z];
    float *gal_chunksz[NBINS_SUBVOLUME_X][NBINS_SUBVOLUME_Y][NBINS_SUBVOLUME_Z];

    for(int i=0;i<NBINS_SUBVOLUME_X;i++) {
        for(int j=0;j<NBINS_SUBVOLUME_Y;j++) {
            for(int k=0;k<NBINS_SUBVOLUME_Z;k++) {
                gal_chunksx[i][j][k] = (float*)malloc(MAX_NGALS_IN_SUBVOLUME*sizeof(float));
                gal_chunksy[i][j][k] = (float*)malloc(MAX_NGALS_IN_SUBVOLUME*sizeof(float));
                gal_chunksz[i][j][k] = (float*)malloc(MAX_NGALS_IN_SUBVOLUME*sizeof(float));
            }
        }
    }

    int num_in_gal_chunks[NBINS_SUBVOLUME_X][NBINS_SUBVOLUME_Y][NBINS_SUBVOLUME_Z];

    for(int i=0;i<NBINS_SUBVOLUME_X;i++) {
        for(int j=0;j<NBINS_SUBVOLUME_Y;j++) {
            for(int k=0;k<NBINS_SUBVOLUME_Z;k++) {
                num_in_gal_chunks[i][j][k] = 0;
                for(int n=0;n<MAX_NGALS_IN_SUBVOLUME;n++) {
                    gal_chunksx[i][j][k][n] = 0.0;
                    gal_chunksy[i][j][k][n] = 0.0;
                    gal_chunksz[i][j][k][n] = 0.0;
                    //printf("%f ",gal_chunksx[i][j][k][n]);
                }
            }
        }
    }

    //exit(1);

    int ii,jj,kk;
    int index = 0;

    for(int i=0;i<NGALS;i++) {
        ii = gal_chunk_coordinate[i][0];
        jj = gal_chunk_coordinate[i][1];
        kk = gal_chunk_coordinate[i][2];

        index = num_in_gal_chunks[ii][jj][kk];

        gal_chunksx[ii][jj][kk][index] = galsx[i];
        gal_chunksy[ii][jj][kk][index] = galsy[i];
        gal_chunksz[ii][jj][kk][index] = galsz[i];

        //printf("%f\n",galsx[i]);
        //printf("%f\n",gal_chunksx[ii][jj][kk][index]);

        num_in_gal_chunks[ii][jj][kk]++;
        if(num_in_gal_chunks[ii][jj][kk]>MAX_NGALS_IN_SUBVOLUME)
        {
            printf("Too many galaxies in a subvolume!!! %d",MAX_NGALS_IN_SUBVOLUME);
            exit(-1);
        }
    }

    //exit(1);

    for(int i=0;i<NBINS_SUBVOLUME_X;i++) {
        printf("----------- %d\n",i);
        for(int j=0;j<NBINS_SUBVOLUME_Y;j++) {
            for(int k=0;k<NBINS_SUBVOLUME_Z;k++) {
                printf("%d ",num_in_gal_chunks[i][j][k]);
            }
            printf("\n");
        }
    }

    //exit(1);
    // Test the chunking.
    float x,y,z;
    int ngals_in_chunk;
    for(int i=0;i<NBINS_SUBVOLUME_X;i++) 
        for(int j=0;j<NBINS_SUBVOLUME_Y;j++) 
            for(int k=0;k<NBINS_SUBVOLUME_Z;k++) {
                ngals_in_chunk = num_in_gal_chunks[i][j][k];
                    for(int n=0;n<ngals_in_chunk;n++) {
                        //printf("%d\n",n);
                        x = gal_chunksx[i][j][k][n];
                        y = gal_chunksy[i][j][k][n];
                        z = gal_chunksz[i][j][k][n];
                        printf("%d %d %d %d     %f %f %f\n",i,j,k,n,x,y,z);
                        //printf("%d %d %d %d %d\n",i,j,k,ngals_in_chunk,n);
                        //printf("%f %f %f\n",x,y,z);

                }
            }


    free(galsx);
    free(galsy);
    free(galsz);
    for(int i=0;i<NBINS_SUBVOLUME_X;i++) {
        for(int j=0;j<NBINS_SUBVOLUME_Y;j++) {
            for(int k=0;k<NBINS_SUBVOLUME_Z;k++) {
                free(gal_chunksx[i][j][k]);
                free(gal_chunksy[i][j][k]);
                free(gal_chunksz[i][j][k]);
            }
        }
    }

    return 0;
}
