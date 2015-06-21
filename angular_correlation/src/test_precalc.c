#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<unistd.h>

#include "precalc_distances.h"

#define NGALS 600000
#define MAX_NGALS_IN_SUBVOLUME 1000
#define MAX_DISTANCES (NGALS/1000)*NGALS
#define NBINS_SUBVOLUME_X 16
#define NBINS_SUBVOLUME_Y 29
#define NBINS_SUBVOLUME_Z 16

int compare (const void * a, const void * b);

int main(int argc, char **argv)
{

    float maxsep;
    maxsep = 150;
    //float vec_ranges_lo[3] = {-2000., -1720., 755.};
    //float vec_ranges_hi[3] = {140., 1780., 2500.};
    float vec_ranges_lo[3] = {1000000.,1000000.,1000000.};
    float vec_ranges_hi[3] = {-1e9,-1e9,-1e9};
    float vec_ranges_width[3] = {0.0,0.0,0.0};
    int vec_nbins[3] = {0,0,0};
    float vec_binwidths[3] = {0.0,0.0,0.0};

    extern char *optarg;
    extern int optind, optopt, opterr;

    int i=0,j=0,k=0,ii=0,jj=0,kk=0,n=0;

    ////////////////////////////////////////////////////////////////////////////
    // Allocate some memory.
    ////////////////////////////////////////////////////////////////////////////
    printf("Allocating memory for the galaxy coordinates...\n");
    float *galsx0, *galsy0, *galsz0;
    float *galsx1, *galsy1, *galsz1;
    galsx0 = (float*)malloc(NGALS*sizeof(float));
    galsy0 = (float*)malloc(NGALS*sizeof(float));
    galsz0 = (float*)malloc(NGALS*sizeof(float));

    galsx1 = (float*)malloc(NGALS*sizeof(float));
    galsy1 = (float*)malloc(NGALS*sizeof(float));
    galsz1 = (float*)malloc(NGALS*sizeof(float));

    printf("Allocated memory for the galaxy coordinates.\n");

    for(i=0;i<3;i++) {
        vec_ranges_width[i] = vec_ranges_hi[i]-vec_ranges_lo[i];
    }
    for(i=0;i<NGALS;i++) {
        galsx0[i]=0;
        galsy0[i]=0;
        galsz0[i]=0;
        galsx1[i]=0;
        galsy1[i]=0;
        galsz1[i]=0;
    }

    int NUM_GALAXIES[2] = {NGALS,NGALS};

    ////////////////////////////////////////////////////////////////////////////
    // Read in galaxies.
    ////////////////////////////////////////////////////////////////////////////
    ///*
    printf("Reading in the galaxies....\n");
    float temp0, temp1, temp2, dummy, tempdum0, tempdum1, tempdum2;

    FILE *infile[2], *outfile ;
        for (i=0;i<2;i++) {
            infile[i] = fopen(argv[optind+i],"r");
            printf("Opening input file %d: %s\n",i,argv[optind+i]);
        }

    NUM_GALAXIES[0] = 0;
    NUM_GALAXIES[1] = 0;

    int max = 800000;

    read_in_6_cols(infile[0], NGALS, galsx0, galsy0, galsz0, &NUM_GALAXIES[0], vec_ranges_lo, vec_ranges_hi, max);
    read_in_6_cols(infile[1], NGALS, galsx1, galsy1, galsz1, &NUM_GALAXIES[1], vec_ranges_lo, vec_ranges_hi, max);

    printf("NUM_GALAXIES[0]: %d\n",NUM_GALAXIES[0]);
    printf("NUM_GALAXIES[1]: %d\n",NUM_GALAXIES[1]);

    for(i=0;i<3;i++) {
        printf("lo/hi: %d %f %f\n",i,vec_ranges_lo[i],vec_ranges_hi[i]);
    }

    ////////////////////////////////////////////////////////////////////////////
    // Define the ranges for our galaxy chunks
    ////////////////////////////////////////////////////////////////////////////
    printf("Define the ranges...\n");
    define_ranges(vec_ranges_lo, vec_ranges_hi, maxsep, vec_nbins, vec_binwidths);

    for(i=0;i<3;i++) {
        printf("%d %d %f\n",i,vec_nbins[i],vec_binwidths[i]);
    }

    ////////////////////////////////////////////////////////////////////////////
    // Assign a chunk coordinate to each
    ////////////////////////////////////////////////////////////////////////////
    printf("Assign chunk coordinates...\n");
    
    int *gal_chunk_coordinate0;
    gal_chunk_coordinate0 = (int*)malloc(NUM_GALAXIES[0]*3*sizeof(int));

    int *gal_chunk_coordinate1;
    gal_chunk_coordinate1 = (int*)malloc(NUM_GALAXIES[0]*3*sizeof(int));


    float coord[3] = {0.0, 0.0, 0.0};
    int chunk_coord[3] = {0,0,0};

    for(i=0;i<NUM_GALAXIES[0];i++) {
        if(i%100000==0)
        {
            printf("%d\n",i);
            fflush(stdout);
        }
        coord[0] = galsx0[i];
        coord[1] = galsy0[i];
        coord[2] = galsz0[i];
        
        assign_chunk_coordinate(coord, vec_ranges_lo, vec_binwidths, chunk_coord);
        for(j=0;j<3;j++) {
            gal_chunk_coordinate0[i*3 + j] = chunk_coord[j];
        }
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

    for(i=0;i<NBINS_SUBVOLUME_X;i++) {
        for(j=0;j<NBINS_SUBVOLUME_Y;j++) {
            num_in_gal_chunks[i][j] = (int*)malloc(NBINS_SUBVOLUME_Z*sizeof(int)); 
            for(k=0;k<NBINS_SUBVOLUME_Z;k++) {
                gal_chunksx[i][j][k] = (float*)malloc(MAX_NGALS_IN_SUBVOLUME*sizeof(float));
                gal_chunksy[i][j][k] = (float*)malloc(MAX_NGALS_IN_SUBVOLUME*sizeof(float));
                gal_chunksz[i][j][k] = (float*)malloc(MAX_NGALS_IN_SUBVOLUME*sizeof(float));
            }
        }
    }

    printf("Allocated the memory....\n");


    // Set everything equal to 0
    for(i=0;i<NBINS_SUBVOLUME_X;i++) {
        for(j=0;j<NBINS_SUBVOLUME_Y;j++) {
            for(k=0;k<NBINS_SUBVOLUME_Z;k++) {
                num_in_gal_chunks[i][j][k] = 0;
                for(n=0;n<MAX_NGALS_IN_SUBVOLUME;n++) {
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

    ii,jj,kk;
    int index = 0;

    for(i=0;i<NUM_GALAXIES[0];i++) {

        ii = gal_chunk_coordinate0[i*3 + 0];
        jj = gal_chunk_coordinate0[i*3 + 1];
        kk = gal_chunk_coordinate0[i*3 + 2];


        index = num_in_gal_chunks[ii][jj][kk];

        //printf("%d %d %d %d %f %f %f\n",ii,jj,kk,index,galsx0[i],galsy0[i],galsz0[i]);

        gal_chunksx[ii][jj][kk][index] = galsx0[i];
        gal_chunksy[ii][jj][kk][index] = galsy0[i];
        gal_chunksz[ii][jj][kk][index] = galsz0[i];

        //printf("%f\n",galsx0[i]);
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

    //for(i=0;i<NBINS_SUBVOLUME_X;i++) {
    //printf("----------- %d\n",i);
    //for(j=0;j<NBINS_SUBVOLUME_Y;j++) {
    //for(k=0;k<NBINS_SUBVOLUME_Z;k++) {
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
       for(i=0;i<NBINS_SUBVOLUME_X;i++) 
       for(j=0;j<NBINS_SUBVOLUME_Y;j++) 
       for(k=0;k<NBINS_SUBVOLUME_Z;k++) {
       ngals_in_chunk = num_in_gal_chunks[i][j][k];
       for(n=0;n<ngals_in_chunk;n++) {
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
    for(i=0;i<30*MAX_NGALS_IN_SUBVOLUME;i++) {
        gals_superchunk[i] = (float*)malloc(3*sizeof(float));
    }


    int n0,n1,ntemp,nindex;
    int iindex,jindex,kindex;

    float *distances;
    printf("MAX_DISTANCES: %d M\n",MAX_DISTANCES/1000000);
    distances = (float*)malloc(MAX_DISTANCES*sizeof(float));
    if(distances)
        printf("Success!\n");
    else
        printf("Booooo!\n");

    float *unique_distances;
    unique_distances = (float*)malloc(MAX_DISTANCES*sizeof(float));

    int ndist = 0;
    float disttemp = 0;

    float testdist = 16383.5068359375;
    int x0,y0,z0,x1,y1,z1;
    int igal0,igal1;
    int igal1min = 0;

    for(i=0;i<NBINS_SUBVOLUME_X;i++) {
        printf("i: %d\n",i);
        for(j=0;j<NBINS_SUBVOLUME_Y;j++) {
            for(k=0;k<NBINS_SUBVOLUME_Z;k++) {

                n0 = num_in_gal_chunks[i][j][k];

                // Only do the calculation if there are galaxies in the chunk.
                if(n0>0)
                {
                    for(ii=0;ii<2;ii++) {
                        for(jj=0;jj<2;jj++) {
                            for(kk=0;kk<2;kk++) {

                                iindex = i+ii;
                                jindex = j+jj;
                                kindex = k+kk;

                                if (iindex>=0 && jindex>=0 && kindex>=0 && \
                                        iindex<NBINS_SUBVOLUME_X && \
                                        jindex<NBINS_SUBVOLUME_Y && \
                                        kindex<NBINS_SUBVOLUME_Z) {

                                    n1 = num_in_gal_chunks[iindex][jindex][kindex];

                                    if(n1>0)
                                    {

                                        for(igal0=0;igal0<n0;igal0++) {

                                            x0 = gal_chunksx[i][j][k][igal0];
                                            y0 = gal_chunksy[i][j][k][igal0];
                                            z0 = gal_chunksz[i][j][k][igal0];

                                            igal1min = 0;
                                            if (ii==0 && jj==0 && kk==0)
                                                igal1min = igal0+1;

                                            //printf("igal1min: %d\n",igal1min);
                                            for(igal1=igal1min;igal1<n1;igal1++) {

                                                x1 = gal_chunksx[iindex][jindex][kindex][igal1];
                                                y1 = gal_chunksy[iindex][jindex][kindex][igal1];
                                                z1 = gal_chunksz[iindex][jindex][kindex][igal1];

                                                if(i==15)
                                                    printf("%d %d %d %d %f %f %f\n",igal0,igal1,n0,n1,x0,y0,z0);

                                                disttemp = distance(x0,y0,z0,x1,y1,z1);

                                                if (disttemp==testdist) {
                                                    printf("%.18f\n",disttemp);
                                                    printf("%d %d %d %f %f %f %f %f %f\n",igal0,jj,n0,x0,y0,z0,x1,y1,z1);
                                                }
                                                //printf("disttemp: %32.32f\n",disttemp);
                                                if (disttemp>0 && disttemp<maxsep) {
                                                    if(ndist>=MAX_DISTANCES) {
                                                        printf("The number of calculated distances (%d) had exceeded the memory allotted (%d).\n",ndist,MAX_DISTANCES);
                                                        exit(-1);
                                                    }
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
                                    {
                                    }
                                }
                            }
                        }
                    }
                    //printf("n0/n1: %d %d\n",n0,n1);

                    //printf("%d %d\n",n0,n1);
                }

            }
        }
    }

    printf("Calculated the distances.\n");

    /*
       printf("Copying over the distances,\n");
       float *temp_distances;
       temp_distances = (float*)malloc(ndist*sizeof(float));
       for(i=0;i<ndist;i++) {
       temp_distances[i] = 0.;
       }

       for(i=0;i<ndist;i++) {
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
    /*
       printf("About to qsort....\n");
       qsort(distances,ndist,sizeof(float),compare);
       printf("Sorted!!!!\n");
       printf("%f %f %f %f\n",distances[0],distances[10000],distances[1000000],distances[10000000]);

       for(i=0;i<ndist;i++)
       {
    //printf("%f ",distances[i]);
    current_dist=distances[i];
    //if(current_dist>100 && current_dist<101)
    //{
    //printf("%.18f\n",current_dist);
    //}
    if(current_dist>16383.5 && current_dist<16384.0)
    printf("current_dist: %.18f\n",current_dist);

    if(current_dist!=prev_dist)
    {
    unique_distances[final_ndist] = current_dist;
    prev_dist = current_dist;
    final_ndist++;
    }
    else
    {
    //printf("%.10f %.10f\n",current_dist,prev_dist);
    }
    //if(i%30==0)
    //printf("\n");
    }
    printf("ndist/unique: %d %d\n",ndist,final_ndist);
     */
    printf("ndist: %d\n",ndist);

    printf("Histogramming!\n");
    float lo=0.;
    float hi=256.;
    int nbins=256;
    float binwidth=(hi-lo)/nbins;
    long int *histogram;
    histogram = (long int*)malloc(nbins*sizeof(long int));
    // Initialize histogram
    for(i=0;i<nbins;i++){
        histogram[i]=0;
    }

    int bin_index;
    float dist;
    //for(i=0;i<final_ndist;i++){
    for(i=0;i<ndist;i++){
        //dist=unique_distances[i];
        dist=distances[i];
        if(dist>=lo && dist<=hi)
        {
            bin_index = (int)((dist-lo)/binwidth);
            //printf("%f %d\n",dist,bin_index);
            histogram[bin_index]++;
        }
    }

    outfile = fopen(argv[optind+2],"w");
    char *output;
    output = (char*)malloc(1024*sizeof(char));

    fprintf(outfile,"%d\n",NUM_GALAXIES[0]);
    for(i=0;i<nbins;i++){
        //printf("----------------- %d\n",i);
        //printf("%.1f %.1f   %d\n",i*binwidth + lo , (i+1)*binwidth + lo,histogram[i]);
        sprintf(output,"%5.1f %5.1f %7d\n",i*binwidth + lo , (i+1)*binwidth + lo,histogram[i]);
        printf(output);
        fprintf(outfile,output);
    }
    fclose(outfile);

    ////////////////////////////////////////////////////////////////////////////
    // Free up everything. 
    ////////////////////////////////////////////////////////////////////////////
    free(galsx0);
    free(galsy0);
    free(galsz0);
    free(galsx1);
    free(galsy1);
    free(galsz1);
    free(distances);
    free(unique_distances);
    free(histogram);
    free(gal_chunk_coordinate0);
    free(gal_chunk_coordinate1);
    //for(i=0;i<NGALS;i++) 
    //free(gal_chunk_coordinate[i]);

    for(i=0;i<NBINS_SUBVOLUME_X;i++) {
        for(j=0;j<NBINS_SUBVOLUME_Y;j++) {

            free(num_in_gal_chunks[i][j]);

            for(k=0;k<NBINS_SUBVOLUME_Z;k++) {
                free(gal_chunksx[i][j][k]);
                free(gal_chunksy[i][j][k]);
                free(gal_chunksz[i][j][k]);
            }
        }
    }

    for(i=0;i<30*MAX_NGALS_IN_SUBVOLUME;i++)
        free(gals_superchunk[i]);

    fclose(infile[0]);
    fclose(infile[1]);
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
int compare (const void * a, const void * b)
{
    float fa = *(const float*) a;
    float fb = *(const float*) b;
    return (fa > fb) - (fa < fb);
}
///////////////////////////////////////////////////////////////////////////////
