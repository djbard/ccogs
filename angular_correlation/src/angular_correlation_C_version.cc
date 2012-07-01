#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<cmath>

using namespace std;

#define SUBMATRIX_SIZE 10000
//#define NUM_BIN 5000
//#define HIST_MIN 0.0
//#define HIST_MAX 3.5 
#define NUM_BIN 27 // for log binning
//#define NUM_BIN 37 // for log binning
#define HIST_MIN 0.0 // for degrees
#define HIST_MAX 100.0 // for degrees

#define CONV_FACTOR 57.2957795 // 180/pi

////////////////////////////////////////////////////////////////////////
float distance(float a0, float d0, float a1, float d1) 
{

    float alpha = a0, delta0 = d0;
    float cos_d0 = cos(delta0), sin_d0 = sin(delta0), dist;

    float a_diff, sin_a_diff, cos_a_diff;
    float cos_d1, sin_d1, numer, denom, mult1, mult2;    

    a_diff = a1 - alpha;

    sin_a_diff = sin(a_diff);
    cos_a_diff = cos(a_diff);

    sin_d1 = sin(d1);
    cos_d1 = cos(d1);

    mult1 = cos_d1 * cos_d1 * sin_a_diff * sin_a_diff;
    mult2 = cos_d0 * sin_d1 - sin_d0 * cos_d1 * cos_a_diff;
    mult2 = mult2 * mult2;

    numer = sqrt(mult1 + mult2); 

    denom = sin_d0 *sin_d1 + cos_d0 * cos_d1 * cos_a_diff;

    //dist = atan(num);  
    dist = atan2(numer,denom);  
    dist *= CONV_FACTOR;  // Convert to degrees

    return dist;

}

////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{

    float *h_alpha0, *h_delta0;
    float *h_alpha1, *h_delta1;

    int NUM_PARTICLES;

    if (argc < 4)
    {

        printf("\nMust pass in cluster_data file  on command line!\n");
        printf("\nUsage: ", argv[0] );
        printf(" <cluster_data file> <distances file> \n\n");
        exit(1);
    }

    FILE *infile0, *infile1, *outfile ;
    infile0 = fopen(argv[1],"r");
    infile1 = fopen(argv[2],"r");
    outfile = fopen(argv[3], "w");

    bool two_different_files = 1;
    if (strcmp(argv[1],argv[2])==0)
    {
        two_different_files = 0;
        printf("Using the same file!\n");
    }

    //////////////////////////////////////////////////////////////////////
    // Read in the cluster_data file
    ////////////////////////////////////////////////////////////////////////////

    char axis_titles[256];
    char dummy[256];

    ////////////////////////////////////////////////////////////////////////////
    // Read in the first file
    ////////////////////////////////////////////////////////////////////////////

    fscanf(infile0, "%s %s %s", &axis_titles, &dummy, &axis_titles);
    fscanf(infile0, "%d", &NUM_PARTICLES);

    int size = NUM_PARTICLES * sizeof(float);    
    printf("SIZE0 # particles: %d\n",NUM_PARTICLES);

    h_alpha0 = (float*)malloc(size);
    h_delta0 = (float*)malloc(size);

    for(int i=0; i<NUM_PARTICLES; i++)
    {
        fscanf(infile0, "%f %s %f %s ", &h_alpha0[i], &dummy, &h_delta0[i], &dummy);
    }

    ////////////////////////////////////////////////////////////////////////////
    // Read in the second file
    ////////////////////////////////////////////////////////////////////////////

    fscanf(infile1, "%s %s %s", &axis_titles, &dummy, &axis_titles);
    fscanf(infile1, "%d", &NUM_PARTICLES);
    printf("SIZE1 # particles: %d\n",NUM_PARTICLES);

    h_alpha1 = (float*)malloc(size);
    h_delta1 = (float*)malloc(size);

    for(int i=0; i<NUM_PARTICLES; i++)
    {
        fscanf(infile1, "%f %s %f %s ", &h_alpha1[i], &dummy, &h_delta1[i], &dummy);
    }

    ////////////////////////////////////////////////////////////////////////////
    //allocation of histogram
    ///////////////////////////////////////////////////////////////////////////

    int *hist, *dev_hist;

    // Log binning
    //float h_bin_edges[30] = {0.001000,0.001585,0.002512,0.003981,0.006310,0.010000,0.010000,0.015849,0.025119,0.039811,0.063096,0.100000,0.100000,0.158489,0.251189,0.398107,0.630957,1.000000,1.000000,1.584893,2.511886,3.981072,6.309573,10.000000,10.000000,15.848932,25.118864,39.810717,63.095734,100.000000};

    // For 27 bins
    float h_bin_edges[NUM_BIN] = {0.0000,0.001000,0.001585,0.002512,0.003981,0.006310,0.010000,0.015849,0.025119,0.039811,0.063096,0.100000,0.158489,0.251189,0.398107,0.630957,1.000000,1.584893,2.511886,3.981072,6.309573,10.000000,15.848932,25.118864,39.810717,63.095734,100.000000};

    // For 37 bins
    //float h_bin_edges[NUM_BIN] = {0.0000,0.001000,0.001389,0.001931,0.002683,0.003728,0.005179,0.007197,0.010000,0.013895,0.019307,0.026827,0.037276,0.051795,0.071969,0.100000,0.138950,0.193070,0.268270,0.372759,0.517947,0.719686,1.000000,1.389495,1.930698,2.682696,3.727594,5.179475,7.196857,10.000000,13.894955,19.306977,26.826958,37.275937,51.794747,71.968567,100.000000};

    /*
       for (int i=0;i<NUM_BIN;i++)
       {
       printf("%d %f\n",i,h_bin_edges[i]);
       }
       printf("\n");
     */

    int size_hist = (NUM_BIN+2);
    int size_hist_bytes = size_hist*sizeof(int);

    hist = (int*)malloc(size_hist_bytes);
    memset(hist, 0, size_hist_bytes);

    int x, y;
    float dist = 0;

    int bin_index = 0;
    for(int i = 0; i < NUM_PARTICLES; i++)
    {
        if (i%1000==0)
        {
        printf("%d\n",i);
        }

        for(int j = 0; j < NUM_PARTICLES; j++)
        {
            //printf("----\n");
            //printf("%d %d\t\t%d %d\n",k,y,j,x);
            //printf("----\n");


            bool do_calc = 1;
            if (two_different_files)
            {
                do_calc = 1;
            }
            else // Doing the same file
            {
                if(i > j)
                    do_calc=1;
                else
                    do_calc=0;
            }
            //if(idx > i) ///////// CHECK THIS
            if (do_calc)
            {
                dist = distance(h_alpha0[i],h_delta0[i],h_alpha1[j],h_delta1[j]);

                    if(dist < HIST_MIN)
                        bin_index = 0; 
                    else if(dist >= HIST_MAX)
                        bin_index = NUM_BIN + 1;
                    else
                    {
                        //bin_index = int(((dist - HIST_MIN) * NUM_BIN / HIST_MAX) +1);    
                        bin_index = 0;
                        for (int k=0;k<NUM_BIN-1;k++)
                        {
                            //bin_index = 5;
                            //if (dist>=0.1*j && dist<0.1*(j+1))
                            //if (dist>=dev_bin_edges[j] && dist<dev_bin_edges[j+1])
                            if (dist>=h_bin_edges[k] && dist<h_bin_edges[k+1])
                            {
                                bin_index = k+1;
                                break;
                            }
                        }
                    }

                hist[bin_index]++;
            }
        }
    }  

    unsigned long total = 0;
    float  bin_width = (HIST_MAX - HIST_MIN) / NUM_BIN;
    float bins_mid = 0;

    fprintf(outfile, "%s %s\n", "Angular Distance(radians)","Number of Entries");      
    for(int k=0; k<NUM_BIN+1; k++)
    {
        //bins_mid = bin_width*(k - 0.5);

        float lo = h_bin_edges[k];
        float hi = h_bin_edges[k+1];

        bins_mid = (hi+lo)/2.0;

        fprintf(outfile, "%.3e %s %lu \n", bins_mid, ",",  hist[k]);
        total += hist[k];

    }
    printf("total: %lu \n", total);

    fclose(infile0);
    fclose(infile1);
    fclose(outfile);

    free(h_alpha0);
    free(h_delta0);
    free(h_alpha1);
    free(h_delta1);
    free(hist);

    return 0;
}  
//////////////////////////////////////////////////////////////////////
