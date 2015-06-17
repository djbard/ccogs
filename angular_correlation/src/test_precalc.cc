//#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<unistd.h>

#include "precalc_distances.h"

#define NGALS 100000
#define MAX_NGALS_IN_SUBVOLUME 600
#define MAX_DISTANCES (NGALS/1000)*NGALS
#define NBINS_SUBVOLUME_X 15
#define NBINS_SUBVOLUME_Y 25
#define NBINS_SUBVOLUME_Z 12

///////////////////////////////////////////////////////////////////////////////
int compare (const void * a, const void * b)
{
      float fa = *(const float*) a;
        float fb = *(const float*) b;
          return (fa > fb) - (fa < fb);
}
///////////////////////////////////////////////////////////////////////////////

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
    printf("Define the ranges...\n");
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
    printf("Gen fake galaxies....\n");
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
    printf("Assign chunk coordinates...\n");
    //int gal_chunk_coordinate[NGALS][3];
    int *gal_chunk_coordinate[NGALS];
    for(int i=0;i<NGALS;i++) {
        gal_chunk_coordinate[i] = (int*)malloc(30*sizeof(int));
    }

    float coord[3] = {0.0, 0.0, 0.0};
    int chunk_coord[3] = {0,0,0};
    for(int i=0;i<NGALS;i++) {
        if(i%100000==0)
        {
            printf("%d\n",i);
            fflush(stdout);
        }
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
    printf("Divide up the chunked voxels....\n");
    float *gal_chunksx[NBINS_SUBVOLUME_X][NBINS_SUBVOLUME_Y][NBINS_SUBVOLUME_Z];
    float *gal_chunksy[NBINS_SUBVOLUME_X][NBINS_SUBVOLUME_Y][NBINS_SUBVOLUME_Z];
    float *gal_chunksz[NBINS_SUBVOLUME_X][NBINS_SUBVOLUME_Y][NBINS_SUBVOLUME_Z];

    int *num_in_gal_chunks[NBINS_SUBVOLUME_X][NBINS_SUBVOLUME_Y]; //[NBINS_SUBVOLUME_Z];

    for(int i=0;i<NBINS_SUBVOLUME_X;i++) {
        for(int j=0;j<NBINS_SUBVOLUME_Y;j++) {
            num_in_gal_chunks[i][j] = (int*)malloc(NBINS_SUBVOLUME_Z*sizeof(int)); 
            for(int k=0;k<NBINS_SUBVOLUME_Z;k++) {
                gal_chunksx[i][j][k] = (float*)malloc(MAX_NGALS_IN_SUBVOLUME*sizeof(float));
                gal_chunksy[i][j][k] = (float*)malloc(MAX_NGALS_IN_SUBVOLUME*sizeof(float));
                gal_chunksz[i][j][k] = (float*)malloc(MAX_NGALS_IN_SUBVOLUME*sizeof(float));
            }
        }
    }


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

    //for(int i=0;i<NBINS_SUBVOLUME_X;i++) {
        //printf("----------- %d\n",i);
        //for(int j=0;j<NBINS_SUBVOLUME_Y;j++) {
            //for(int k=0;k<NBINS_SUBVOLUME_Z;k++) {
                //printf("%d ",num_in_gal_chunks[i][j][k]);
            //}
            //printf("\n");
        //}
    //}

    //exit(1);
    // Test the chunking.
    float x,y,z;
    int ngals_in_chunk;
    /*
    for(int i=0;i<NBINS_SUBVOLUME_X;i++) 
        for(int j=0;j<NBINS_SUBVOLUME_Y;j++) 
            for(int k=0;k<NBINS_SUBVOLUME_Z;k++) {
                ngals_in_chunk = num_in_gal_chunks[i][j][k];
                    for(int n=0;n<ngals_in_chunk;n++) {
                        //printf("%d\n",n);
                        x = gal_chunksx[i][j][k][n];
                        y = gal_chunksy[i][j][k][n];
                        z = gal_chunksz[i][j][k][n];
                        //printf("%d %d %d %d     %f %f %f\n",i,j,k,n,x,y,z);
                        //printf("%d %d %d %d %d\n",i,j,k,ngals_in_chunk,n);
                        //printf("%f %f %f\n",x,y,z);

                }
            }
    */

    printf("Calculating distances....\n");
    ///////////////////////////////////////////////////////////////////////////
    // Calc distances with only the chunks.
    ///////////////////////////////////////////////////////////////////////////
    float *gals_superchunk[30*MAX_NGALS_IN_SUBVOLUME];
    for(int i=0;i<30*MAX_NGALS_IN_SUBVOLUME;i++) {
        gals_superchunk[i] = (float*)malloc(3*sizeof(float));
    }


    int n0,n1,ntemp,nindex;
    int iindex,jindex,kindex;

    float *distances;
    printf("MAX_DISTANCES: %f\n",MAX_DISTANCES/1e6);
    distances = (float*)malloc(MAX_DISTANCES*sizeof(float));

    float *unique_distances;
    unique_distances = (float*)malloc(MAX_DISTANCES*sizeof(float));

    int ndist = 0;
    float disttemp = 0;

    for(int i=0;i<NBINS_SUBVOLUME_X;i++) {
        printf("i: %d\n",i);
        for(int j=0;j<NBINS_SUBVOLUME_Y;j++) {
            for(int k=0;k<NBINS_SUBVOLUME_Z;k++) {

                n0 = num_in_gal_chunks[i][j][k];

                n1 = 0;
                nindex = 0;
                for(int ii=-1;ii<2;ii++) {
                    for(int jj=-1;jj<2;jj++) {
                        for(int kk=-1;kk<2;kk++) {
                            ntemp = num_in_gal_chunks[i][j][k];
                            for (int nindex=0;nindex<ntemp;nindex++)
                            {
                                iindex = i+ii;
                                jindex = j+jj;
                                kindex = k+kk;
                                if (iindex>=0 && jindex>=0 && kindex>=0 && \
                                        iindex<NBINS_SUBVOLUME_X && \
                                        jindex<NBINS_SUBVOLUME_Y && \
                                        kindex<NBINS_SUBVOLUME_Z) {
                                    gals_superchunk[n1][0] = gal_chunksx[iindex][jindex][kindex][nindex];
                                    gals_superchunk[n1][1] = gal_chunksy[iindex][jindex][kindex][nindex];
                                    gals_superchunk[n1][2] = gal_chunksz[iindex][jindex][kindex][nindex];
                                    n1++;
                                }
                            }
                        }
                    }
                }
                //printf("n0/n1: %d %d\n",n0,n1);

                //printf("%d %d\n",n0,n1);
                for(int ii=0;ii<n0;ii++) {
                    for(int jj=0;jj<n1;jj++) {
                        //printf("%d %d %d %f %f %f\n",ii,jj,n0, gal_chunksx[i][j][k][ii],gal_chunksy[i][j][k][ii],gal_chunksz[i][j][k][ii]);
                        disttemp = distance(gal_chunksx[i][j][k][ii],gal_chunksy[i][j][k][ii],gal_chunksz[i][j][k][ii],\
                                        gals_superchunk[jj][0],gals_superchunk[jj][1],gals_superchunk[jj][2]);
                        //printf("disttemp: %f\n",disttemp);
                        if (disttemp>0 && disttemp<maxsep) {
                            //printf("here! %d\n",ndist);
                            distances[ndist] = disttemp;
                            //printf("there!\n");
                            ndist++;
                            if (ndist%100000==0)
                                printf("ndist: %d\n",ndist);
                        }
                    }
                }

            }
        }
    }

    // Sort the array
    float prev_dist=-999;
    float current_dist=-999;
    int final_ndist=0;
    qsort(distances,ndist,sizeof(float),compare);
    for(int i=0;i<ndist;i++)
    {
        //printf("%f ",distances[i]);
        current_dist=distances[i];
        if(current_dist!=prev_dist)
        {
            unique_distances[final_ndist] = current_dist;
            prev_dist = current_dist;
            final_ndist++;
        }
        //if(i%30==0)
            //printf("\n");
    }
    printf("ndist/unique: %d %d\n",ndist,final_ndist);

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
