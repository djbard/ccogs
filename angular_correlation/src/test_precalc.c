#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<unistd.h>

#include "precalc_distances.h"

#define NGALS 1000000
#define MAX_NGALS_IN_SUBVOLUME 2000
#define MAX_DISTANCES (NGALS/1000)*NGALS
#define NBINS_SUBVOLUME_X 16
#define NBINS_SUBVOLUME_Y 29
#define NBINS_SUBVOLUME_Z 16

int compare (const void * a, const void * b);

int main(int argc, char **argv)
{

    float maxsep;
    maxsep = 200;
    //float vec_ranges_lo[3] = {-2000., -1720., 755.};
    //float vec_ranges_hi[3] = {140., 1780., 2500.};
    float vec_ranges_lo[3] = {1000000.,1000000.,1000000.};
    float vec_ranges_hi[3] = {-1e9,-1e9,-1e9};
    float vec_ranges_width[3] = {0.0,0.0,0.0};
    int vec_nbins[3] = {0,0,0};
    float vec_binwidths[3] = {0.0,0.0,0.0};

    extern char *optarg;
    extern int optind, optopt, opterr;

    int i=0,j=0,k=0,ii=0,jj=0,kk=0,n=0,nfile=0;

    ////////////////////////////////////////////////////////////////////////////
    // Allocate some memory.
    ////////////////////////////////////////////////////////////////////////////
    printf("Allocating memory for the galaxy coordinates...\n");
    //float *galsx0, *galsy0, *galsz0;
    //float *galsx1, *galsy1, *galsz1;
    float *galsx[2], *galsy[2], *galsz[2];

    for(i=0;i<2;i++) {
        galsx[i] = (float*)malloc(NGALS*sizeof(float));
        galsy[i] = (float*)malloc(NGALS*sizeof(float));
        galsz[i] = (float*)malloc(NGALS*sizeof(float));
    }
    //galsx0 = (float*)malloc(NGALS*sizeof(float));
    //galsy0 = (float*)malloc(NGALS*sizeof(float));
    //galsz0 = (float*)malloc(NGALS*sizeof(float));

    //galsx1 = (float*)malloc(NGALS*sizeof(float));
    //galsy1 = (float*)malloc(NGALS*sizeof(float));
    //galsz1 = (float*)malloc(NGALS*sizeof(float));

    printf("Allocated memory for the galaxy coordinates.\n");

    for(i=0;i<3;i++) {
        vec_ranges_width[i] = vec_ranges_hi[i]-vec_ranges_lo[i];
    }
    for(i=0;i<2;i++) {
        for(j=0;j<NGALS;j++) {
            galsx[i][j]=0;
            galsy[i][j]=0;
            galsz[i][j]=0;
        }
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

    for(i=0;i<2;i++) {
        read_in_6_cols(infile[i], NGALS, galsx[i], galsy[i], galsz[i], &NUM_GALAXIES[i], vec_ranges_lo, vec_ranges_hi, max);
    }

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

    //int *gal_chunk_coordinate0;
    //gal_chunk_coordinate0 = (int*)malloc(NUM_GALAXIES[0]*3*sizeof(int));

    //int *gal_chunk_coordinate1;
    //gal_chunk_coordinate1 = (int*)malloc(NUM_GALAXIES[0]*3*sizeof(int));

    int *gal_chunk_coordinate[2];
    for(i=0;i<2;i++) {
        gal_chunk_coordinate[i] = (int*)malloc(NUM_GALAXIES[i]*3*sizeof(int));
    }


    float coord[3] = {0.0, 0.0, 0.0};
    int chunk_coord[3] = {0,0,0};

    for(i=0;i<2;i++) {
        for(j=0;j<NUM_GALAXIES[i];j++) {
            if(j%100000==0) {
                printf("%d\n",j);
                fflush(stdout);
            }
            coord[0] = galsx[i][j];
            coord[1] = galsy[i][j];
            coord[2] = galsz[i][j];

            assign_chunk_coordinate(coord, vec_ranges_lo, vec_binwidths, chunk_coord);
            for(k=0;k<3;k++) {
                gal_chunk_coordinate[i][j*3 + k] = chunk_coord[k];
            }
        }
    }

    //exit(-1);

    ////////////////////////////////////////////////////////////////////////////
    // Divide up the galaxies into the chunked voxels.
    ////////////////////////////////////////////////////////////////////////////
    printf("Divide up the chunked voxels....\n");
    float *gal_chunksx[2][NBINS_SUBVOLUME_X][NBINS_SUBVOLUME_Y][NBINS_SUBVOLUME_Z];
    float *gal_chunksy[2][NBINS_SUBVOLUME_X][NBINS_SUBVOLUME_Y][NBINS_SUBVOLUME_Z];
    float *gal_chunksz[2][NBINS_SUBVOLUME_X][NBINS_SUBVOLUME_Y][NBINS_SUBVOLUME_Z];

    // Allocate memory.
    int *num_in_gal_chunks[2][NBINS_SUBVOLUME_X][NBINS_SUBVOLUME_Y]; //[NBINS_SUBVOLUME_Z];

    for(nfile=0;nfile<2;nfile++) {
        for(i=0;i<NBINS_SUBVOLUME_X;i++) {
            for(j=0;j<NBINS_SUBVOLUME_Y;j++) {
                num_in_gal_chunks[nfile][i][j] = (int*)malloc(NBINS_SUBVOLUME_Z*sizeof(int)); 
                for(k=0;k<NBINS_SUBVOLUME_Z;k++) {
                    gal_chunksx[nfile][i][j][k] = (float*)malloc(MAX_NGALS_IN_SUBVOLUME*sizeof(float));
                    gal_chunksy[nfile][i][j][k] = (float*)malloc(MAX_NGALS_IN_SUBVOLUME*sizeof(float));
                    gal_chunksz[nfile][i][j][k] = (float*)malloc(MAX_NGALS_IN_SUBVOLUME*sizeof(float));
                }
            }
        }
    }

    printf("Allocated the memory....\n");


    // Set everything equal to 0
    for(nfile=0;nfile<2;nfile++) {
        for(i=0;i<NBINS_SUBVOLUME_X;i++) {
            for(j=0;j<NBINS_SUBVOLUME_Y;j++) {
                for(k=0;k<NBINS_SUBVOLUME_Z;k++) {
                    num_in_gal_chunks[nfile][i][j][k] = 0;
                    for(n=0;n<MAX_NGALS_IN_SUBVOLUME;n++) {
                        gal_chunksx[nfile][i][j][k][n] = 0.0;
                        gal_chunksy[nfile][i][j][k][n] = 0.0;
                        gal_chunksz[nfile][i][j][k][n] = 0.0;
                        //printf("%f ",gal_chunksx[i][j][k][n]);
                    }
                }
            }
        }
    }

    printf("Everything set to 0.\n");

    //exit(1);

    ii,jj,kk;
    int index = 0;

    printf("Building the chunks of galaxies.\n");

    for(nfile=0;nfile<2;nfile++) {
        for(i=0;i<NUM_GALAXIES[nfile];i++) {

            ii = gal_chunk_coordinate[nfile][i*3 + 0];
            jj = gal_chunk_coordinate[nfile][i*3 + 1];
            kk = gal_chunk_coordinate[nfile][i*3 + 2];

            index = num_in_gal_chunks[nfile][ii][jj][kk];

            //printf("%d %d %d %d %f %f %f\n",ii,jj,kk,index,galsx0[i],galsy0[i],galsz0[i]);

            gal_chunksx[nfile][ii][jj][kk][index] = galsx[nfile][i];
            gal_chunksy[nfile][ii][jj][kk][index] = galsy[nfile][i];
            gal_chunksz[nfile][ii][jj][kk][index] = galsz[nfile][i];

            //printf("%f\n",galsx0[i]);
            //printf("%f %f %f (%d,%d,%d) %d\n",gal_chunksx[ii][jj][kk][index],gal_chunksy[ii][jj][kk][index],gal_chunksz[ii][jj][kk][index],ii,jj,kk,index);

            num_in_gal_chunks[nfile][ii][jj][kk]++;

            if(num_in_gal_chunks[nfile][ii][jj][kk]>MAX_NGALS_IN_SUBVOLUME)
            {
                printf("Too many galaxies in a subvolume!!! %d\n",MAX_NGALS_IN_SUBVOLUME);
                printf("%d %d %d %d %d %d\n",nfile,i,ii,jj,kk,index);
                exit(-1);
            }
        }
    }

    printf("Chunked up the galaxies.\n");
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    // Calc distances with only the chunks.
    ///////////////////////////////////////////////////////////////////////////
    printf("Calculating distances....\n");

    float x,y,z;
    int ngals_in_chunk;

    int n0,n1,ntemp,nindex;
    int iindex,jindex,kindex;

    float *distances;
    printf("MAX_DISTANCES: %d M\n",MAX_DISTANCES/1000000);
    distances = (float*)malloc(MAX_DISTANCES*sizeof(float));
    if(distances)
        printf("Success!\n");
    else {
        printf("Booooo!\n");
        exit(-1);
    }

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

                n0 = num_in_gal_chunks[0][i][j][k];

                // Only do the calculation if there are galaxies in the chunk.
                if(n0>0) {
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

                                    n1 = num_in_gal_chunks[1][iindex][jindex][kindex];

                                    if(n1>0) {
                                        for(igal0=0;igal0<n0;igal0++) {

                                            // Gals from first file.
                                            x0 = gal_chunksx[0][i][j][k][igal0];
                                            y0 = gal_chunksy[0][i][j][k][igal0];
                                            z0 = gal_chunksz[0][i][j][k][igal0];

                                            igal1min = 0;
                                            if (ii==0 && jj==0 && kk==0)
                                                igal1min = igal0+1;

                                            //printf("igal1min: %d\n",igal1min);
                                            for(igal1=igal1min;igal1<n1;igal1++) {

                                                x1 = gal_chunksx[1][iindex][jindex][kindex][igal1];
                                                y1 = gal_chunksy[1][iindex][jindex][kindex][igal1];
                                                z1 = gal_chunksz[1][iindex][jindex][kindex][igal1];

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
    fprintf(outfile,"%d\n",NUM_GALAXIES[1]);
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
    for(nfile=0;nfile<2;nfile++) {
        free(galsx[nfile]);
        free(galsy[nfile]);
        free(galsz[nfile]);
        free(gal_chunk_coordinate[nfile]);
    }
    free(distances);
    free(histogram);
    //for(i=0;i<NGALS;i++) 
    //free(gal_chunk_coordinate[i]);

    for(nfile=0;nfile<2;nfile++) {
        for(i=0;i<NBINS_SUBVOLUME_X;i++) {
            for(j=0;j<NBINS_SUBVOLUME_Y;j++) {

                free(num_in_gal_chunks[nfile][i][j]);

                for(k=0;k<NBINS_SUBVOLUME_Z;k++) {
                    free(gal_chunksx[nfile][i][j][k]);
                    free(gal_chunksy[nfile][i][j][k]);
                    free(gal_chunksz[nfile][i][j][k]);
                }
            }
        }
    }

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
