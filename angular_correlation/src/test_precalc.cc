//#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<unistd.h>

#include "precalc_distances.h"

#define NGALS 1000000
#define MAX_NGALS_IN_SUBVOLUME 800
#define MAX_DISTANCES (NGALS/1000)*NGALS
#define NBINS_SUBVOLUME_X 16
#define NBINS_SUBVOLUME_Y 29
#define NBINS_SUBVOLUME_Z 16

///////////////////////////////////////////////////////////////////////////////
int compare (const void * a, const void * b)
{
    float fa = *(const float*) a;
    float fb = *(const float*) b;
    return (fa > fb) - (fa < fb);
}
///////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

    float maxsep = 150.;
    //float vec_ranges_lo[3] = {-2000., -1720., 755.};
    //float vec_ranges_hi[3] = {140., 1780., 2500.};
    float vec_ranges_lo[3] = {1e9,1e9,1e9};
    float vec_ranges_hi[3] = {-1e9,-1e9,-1e9};
    float vec_ranges_width[3] = {0.0,0.0,0.0};
    int vec_nbins[3] = {0,0,0};
    float vec_binwidths[3] = {0.0,0.0,0.0};

    extern char *optarg;
    extern int optind, optopt, opterr;

    ////////////////////////////////////////////////////////////////////////////
    // Gen some fake galaxies.
    ////////////////////////////////////////////////////////////////////////////
    //printf("Gen fake galaxies....\n");
    float *galsx, *galsy, *galsz;
    galsx = (float*)malloc(NGALS*sizeof(float));
    galsy = (float*)malloc(NGALS*sizeof(float));
    galsz = (float*)malloc(NGALS*sizeof(float));
    //
    //
    //for(int i=0;i<NGALS;i++)
    //{
    //galsx[i] = vec_ranges_width[0]*rand()/(float)RAND_MAX + vec_ranges_lo[0];
    //galsy[i] = vec_ranges_width[1]*rand()/(float)RAND_MAX + vec_ranges_lo[1];
    //galsz[i] = vec_ranges_width[2]*rand()/(float)RAND_MAX + vec_ranges_lo[2];
    //
    //}

    //exit(1);

    ////////////////////////////////////////////////////////////////////////////
    // Read in galaxies.
    ////////////////////////////////////////////////////////////////////////////
    printf("Reading in the galaxies....\n");
    float temp0, temp1, temp2, dummy, tempdum;

    FILE *infile0, *infile1, *outfile ;

    infile0 = fopen(argv[optind],"r");
    infile1 = fopen(argv[optind+1],"r");

    int NUM_GALAXIES[2] = {0,0};

    int j = 0;
    //while(fscanf(infile[i], "%d %f %f %f %f %f %f", &idummy, &dummy, &dummy, &dummy, &temp0, &temp1, &temp2) != EOF)
    while(fscanf(infile0, "%e %e %e %e %e %e", &tempdum, &tempdum, &tempdum, &temp0, &temp1, &temp2) != EOF)
    {
        galsx[j] = temp0;///scale_factor;
        galsy[j] = temp1;///scale_factor;
        galsz[j] = temp2;///scale_factor;

        // Keep track of ranges of values
        if(galsx[j]<vec_ranges_lo[0])
            vec_ranges_lo[0]=galsx[j]-1;
        if(galsx[j]>vec_ranges_hi[0])
            vec_ranges_hi[0]=galsx[j]+1;

        if(galsy[j]<vec_ranges_lo[1])
            vec_ranges_lo[1]=galsy[j]-1;
        if(galsy[j]>vec_ranges_hi[1])
            vec_ranges_hi[1]=galsy[j]+1;

        if(galsz[j]<vec_ranges_lo[2])
            vec_ranges_lo[2]=galsz[j]-1;
        if(galsz[j]>vec_ranges_hi[2])
            vec_ranges_hi[2]=galsz[j]+1;


        if (NUM_GALAXIES[0]>=NGALS)
        {
            printf("Exceeded max num galaxies");
            exit(-1);
        }
        if (j<10)
        {
            printf("%f %f %f\n", galsx[j],galsy[j],galsz[j]);
        }
        NUM_GALAXIES[0] += 1;
        j += 1;
    }

    printf("NUM_GALAXIES[0]: %d\n",NUM_GALAXIES[0]);

    for(int i=0;i<3;i++) {
        printf("lo/hi: %d %f %f\n",i,vec_ranges_lo[i],vec_ranges_hi[i]);
    }

    ////////////////////////////////////////////////////////////////////////////
    // Define the ranges for our galaxy chunks
    ////////////////////////////////////////////////////////////////////////////
    printf("Define the ranges...\n");
    define_ranges(vec_ranges_lo, vec_ranges_hi, maxsep, vec_nbins, vec_binwidths);

    for(int i=0;i<3;i++) {
        printf("%d %d %f\n",i,vec_nbins[i],vec_binwidths[i]);
        vec_ranges_width[i] = vec_ranges_hi[i]-vec_ranges_lo[i];
    }

    //fflush(stdout); 
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
    for(int i=0;i<NUM_GALAXIES[0];i++) {
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

    // Allocate memory.
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

    printf("Allocated the memory....\n");


    // Set everything equal to 0
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

    printf("Everything set to 0.\n");

    //exit(1);

    int ii,jj,kk;
    int index = 0;

    for(int i=0;i<NUM_GALAXIES[0];i++) {

        ii = gal_chunk_coordinate[i][0];
        jj = gal_chunk_coordinate[i][1];
        kk = gal_chunk_coordinate[i][2];


        index = num_in_gal_chunks[ii][jj][kk];

        //printf("%d %d %d %d %f %f %f\n",ii,jj,kk,index,galsx[i],galsy[i],galsz[i]);

        gal_chunksx[ii][jj][kk][index] = galsx[i];
        gal_chunksy[ii][jj][kk][index] = galsy[i];
        gal_chunksz[ii][jj][kk][index] = galsz[i];

        //printf("%f\n",galsx[i]);
        //printf("%f %f %f (%d,%d,%d) %d\n",gal_chunksx[ii][jj][kk][index],gal_chunksy[ii][jj][kk][index],gal_chunksz[ii][jj][kk][index],ii,jj,kk,index);

        num_in_gal_chunks[ii][jj][kk]++;

        if(num_in_gal_chunks[ii][jj][kk]>MAX_NGALS_IN_SUBVOLUME)
        {
            printf("Too many galaxies in a subvolume!!! %d\n",MAX_NGALS_IN_SUBVOLUME);
            printf("%d %d %d %d %d\n",i,ii,jj,kk,index);
            exit(-1);
        }
    }

    printf("Chunked up the galaxies.\n");

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
    printf("MAX_DISTANCES: %d M\n",MAX_DISTANCES/1000000);
    distances = (float*)malloc(MAX_DISTANCES*sizeof(float));

    float *unique_distances;
    unique_distances = (float*)malloc(MAX_DISTANCES*sizeof(float));

    int ndist = 0;
    float disttemp = 0;

    for(int i=0;i<NBINS_SUBVOLUME_X;i++) {
        printf("i: %d\n",i);
        for(int j=0;j<NBINS_SUBVOLUME_Y;j++) {
            for(int k=0;k<NBINS_SUBVOLUME_Z;k++) {

                //if(i==15)
                    //printf("%d %d %d\n",i,j,k);
                    //printf("HERE!!!!!");

                n0 = num_in_gal_chunks[i][j][k];

                //if(i==15)
                    //printf("n0: %d     %d %d %d\n",n0,i,j,k);

                // Only do the calculation if there are galaxies in the chunk.
                if(n0>0)
                {

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
                        if(i==15)
                            printf("%d %d %d %f %f %f\n",ii,jj,n0, gal_chunksx[i][j][k][ii],gal_chunksy[i][j][k][ii],gal_chunksz[i][j][k][ii]);
                        disttemp = distance(gal_chunksx[i][j][k][ii],gal_chunksy[i][j][k][ii],gal_chunksz[i][j][k][ii],\
                                gals_superchunk[jj][0],gals_superchunk[jj][1],gals_superchunk[jj][2]);
                        //printf("disttemp: %f\n",disttemp);
                        if (disttemp>0 && disttemp<maxsep) {
                            //printf("here! %d\n",ndist);
                            distances[ndist] = disttemp;
                            //printf("there!\n");
                            ndist++;
                            if (ndist%1000000==0)
                                printf("ndist: %d M\n",ndist/1000000);
                        }
                    }
                }
                }

            }
        }
    }

    printf("Calculated the distances.\n");

    /*
    printf("Copying over the distances,\n");
    float *temp_distances;
    temp_distances = (float*)malloc(ndist*sizeof(float));
    for(int i=0;i<ndist;i++) {
        temp_distances[i] = 0.;
    }

    for(int i=0;i<ndist;i++) {
        temp_distances[i] = distances[i];
        if(distances[i]!=distances[i]){
            printf("BAD!!!! distances[i]: %f\n",distances[i]);
            exit(-1);
        }
        if(i%10000==0)
            printf("%f ",temp_distances[i]);
    }
    printf("\n");
    printf("Copied over the distances,\n");
    */
    printf("ndist: %d\n",ndist);
    
    // Free up some memory.
    //delete distances;

    // Sort the array
    float prev_dist=-999;
    float current_dist=-999;
    int final_ndist=0;
    printf("About to qsort....\n");
    qsort(distances,ndist,sizeof(float),compare);
    printf("Sorted!!!!\n");
    printf("%f %f %f %f\n",distances[0],distances[10000],distances[1000000],distances[10000000]);

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
        else
        {
            printf("%f %f\n",current_dist,prev_dist);
        }
        //if(i%30==0)
        //printf("\n");
    }
    printf("ndist/unique: %d %d\n",ndist,final_ndist);

    printf("Histogramming!\n");
    float lo=0.;
    float hi=256.;
    int nbins=256;
    float binwidth=(hi-lo)/nbins;
    int *histogram;
    histogram = (int*)malloc(nbins*sizeof(int));
    // Initialize histogram
    for(int i=0;i<nbins;i++){
        histogram[i]=0;
    }

    int bin_index;
    float dist;
    for(int i=0;i<final_ndist;i++){
        dist=unique_distances[i];
        if(dist>=lo && dist<=hi)
            bin_index = int((dist-lo)/binwidth);
        histogram[bin_index]++;
    }
    
    for(int i=0;i<nbins;i++){
        printf("%.1f %d\n",i*binwidth + lo + (binwidth/2.),histogram[i]);
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
