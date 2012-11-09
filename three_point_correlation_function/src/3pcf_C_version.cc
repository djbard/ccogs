#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<cmath>
#include <unistd.h>

using namespace std;

#define SUBMATRIX_SIZE 10000
//#define NUM_BIN 5000
//#define HIST_MIN 0.0
//#define HIST_MAX 3.5 
//#define NUM_BIN 27 // for log binning
//#define NUM_BIN 37 // for log binning
#define HIST_MIN 0.0 // for degrees
#define HIST_MAX 100.0 // for degrees

#define DEFAULT_NBINS 254 // for log binning
//#define DEFAULT_NBINS 126 // for log binning
//#define DEFAULT_NBINS 62 // for log binning

#define CONV_FACTOR 57.2957795 // 180/pi

////////////////////////////////////////////////////////////////////////
float distance(float a0, float d0, float a1, float d1,float conv_factor_angle) 
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
    dist *= conv_factor_angle;  // Convert to degrees

    return dist;

}

////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{

    // Needed for parsing command-line arguments.
    extern char *optarg;
    extern int optind, optopt, opterr;
    int c;
    char *filename;
    char *outfilename = NULL;
    char defaultoutfilename[256];
    sprintf(defaultoutfilename,"default_out.dat");
    FILE *binning_file = NULL;

    float hist_lower_range = 0.0000001;
    float hist_upper_range = 0;
    int nbins = DEFAULT_NBINS;
    float hist_bin_width = 0.05;
    int log_binning_flag = 0; // False

    float scale_factor = 1.0; // For if we need to convert input to arcsec or arcmin
    float conv_factor_angle = 57.2957795; // 180/pi // For if we need to convert arcdistance to arcsec or arcmin

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////


    while ((c = getopt(argc, argv, "ao:L:l:w:sm")) != -1) {
        switch(c) {
            case 'L':
                printf("L is set\n");
                hist_lower_range = atof(optarg);
                break;
            case 'w':
                hist_bin_width = atof(optarg);
                printf("Histogram bin width: %f\n",hist_bin_width);
                break;
            case 'l':
                log_binning_flag = atoi(optarg);
                printf("Will use log binning.\n");
                break;
            case 's':
                scale_factor = 206264.0; // To convert arcseconds to radians.
                conv_factor_angle *= 3600.0; // convert radians to arcseconds.
                printf("Reading in values assuming they are arcseconds.\n");
                printf("scale_factor: %f\n",scale_factor);
                printf("conv_factor_angle: %f\n",conv_factor_angle);
                break;
            case 'm':
                scale_factor = 3437.74677; // To convert arcminutes to radians.
                conv_factor_angle *= 60.0; // convert radians to arcminutes.
                printf("scale_factor: %f\n",scale_factor);
                printf("conv_factor_angle: %f\n",conv_factor_angle);
                printf("Reading in values assuming they are arcminutes.\n");
                break;
            case 'o':
                outfilename = optarg;
                printf("Output filename is %s\n", outfilename);
                break;
            case '?':
                printf("unknown arg %c\n", optopt);
                break;
        }
    }

    if (argc < 2)
    {

        printf("\nMust pass in at least two input files on command line!\n");
        printf("\nUsage: ", argv[0] );
        //printf(" <cluster_data file> <distances file> \n\n");
        exit(1);
    }

    // Set a default output file name, if none was passed in on the 
    // command line.
    if (outfilename == NULL) 
    {
        outfilename = defaultoutfilename;
        printf("Output filename is %s\n", outfilename);
    }

    float temp_lo = hist_lower_range;
    if (hist_upper_range == 0)
    {
        if (log_binning_flag==0)
        {
            for (int i=0;i<nbins;i++)
            {
                hist_upper_range = temp_lo + hist_bin_width;
                temp_lo = hist_upper_range;
            }
        }
        else if (log_binning_flag==1)
        {
            for (int i=0;i<nbins;i++)
            {
                hist_upper_range = exp(log(temp_lo) + hist_bin_width);
                temp_lo = hist_upper_range;
            }
        }
        else if (log_binning_flag==2)
        {
            for (int i=0;i<nbins;i++)
            {
                hist_upper_range = pow(10,(log10(temp_lo) + hist_bin_width));
                temp_lo = hist_upper_range;
            }
        }
    }
    printf("hist_upper_range: %f\n",hist_upper_range);

    //FILE *infile0, *infile1, *outfile ;
    //infile0 = fopen(argv[optind],"r");
    //infile1 = fopen(argv[optind+1],"r");


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
    //infile0 = fopen(argv[1],"r");
    //infile1 = fopen(argv[2],"r");
    infile0 = fopen(argv[optind],"r");
    infile1 = fopen(argv[optind+1],"r");

    outfile = fopen(outfilename, "w");

    printf("Opening input file 0: %s\n",argv[optind]);
    printf("Opening input file 1: %s\n",argv[optind+1]);


    bool two_different_files = 1;
    if (strcmp(argv[1],argv[2])==0)
    {
        two_different_files = 0;
        printf("Using the same file!\n");
    }

    // Set a default output file name, if none was passed in on the
    // command line.
    /*
    if (outfilename == NULL)
    {
        outfilename = defaultoutfilename;
        printf("Output filename is %s\n", outfilename);
    }
    */



    //////////////////////////////////////////////////////////////////////
    // Read in the cluster_data file
    ////////////////////////////////////////////////////////////////////////////

    //char axis_titles[256];
    //char dummy[256];

    ////////////////////////////////////////////////////////////////////////////
    // Read in the first file
    ////////////////////////////////////////////////////////////////////////////

    int NUM_GALAXIES0;
    int NUM_GALAXIES1;
    //fscanf(infile0, "%s %s %s", &axis_titles, &dummy, &axis_titles);
    fscanf(infile0, "%d", &NUM_GALAXIES0);

    int size_of_galaxy_array = NUM_GALAXIES0 * sizeof(float);    
    printf("SIZE 0 # GALAXIES: %d\n",NUM_GALAXIES0);

    h_alpha0 = (float*)malloc(size_of_galaxy_array);
    h_delta0 = (float*)malloc(size_of_galaxy_array);
    float temp0, temp1;

    for(int i=0; i<NUM_GALAXIES0; i++)
    {
        fscanf(infile0, "%f %f", &temp0, &temp1);
        h_alpha0[i] = temp0/scale_factor;
        h_delta0[i] = temp1/scale_factor;
        //fscanf(infile0, "%f %f", &h_alpha0[i]*scale_factor, &h_delta0[i]*scale_factor);
        //if (i<10)
        //printf("%e %e\n", h_alpha0[i], h_delta0[i]);
    }

    ////////////////////////////////////////////////////////////////////////////
    // Read in the second file
    ////////////////////////////////////////////////////////////////////////////

    //fscanf(infile1, "%s %s %s", &axis_titles, &dummy, &axis_titles);
    fscanf(infile1, "%d", &NUM_GALAXIES1);
    printf("SIZE 1 # GALAXIES: %d\n",NUM_GALAXIES1);

    size_of_galaxy_array = NUM_GALAXIES1 * sizeof(float);    
    h_alpha1 = (float*)malloc(size_of_galaxy_array);
    h_delta1 = (float*)malloc(size_of_galaxy_array);

    for(int i=0; i<NUM_GALAXIES1; i++)
    {
        fscanf(infile1, "%f %f", &temp0, &temp1);
        h_alpha1[i] = temp0/scale_factor;
        h_delta1[i] = temp1/scale_factor;
        //fscanf(infile1, "%f %f", &h_alpha1[i]*scale_factor, &h_delta1[i]*scale_factor);
        //if (i<10)
        //printf("%e %e\n", h_alpha1[i], h_delta1[i]);
    }

    ////////////////////////////////////////////////////////////////////////////
    // Read in the first file
    ////////////////////////////////////////////////////////////////////////////
    /*

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
     */


    ////////////////////////////////////////////////////////////////////////////
    //allocation of histogram
    ///////////////////////////////////////////////////////////////////////////

    unsigned long *hist;
    float hist_min = hist_lower_range;
    float hist_max = hist_upper_range;
    //int nbins;
    float bin_width = hist_bin_width;
    int log_binning=log_binning_flag;

    int size_hist = (nbins+2);
    int size_hist_bytes = size_hist*sizeof(unsigned long);

    hist = (unsigned long*)malloc(size_hist_bytes);
    memset(hist, 0, size_hist_bytes);

    int x, y;
    float dist = 0;

    int bin_index = 0;
    for(int i = 0; i < NUM_GALAXIES0; i++)
    {
        if (i%1000==0)
        {
            printf("%d\n",i);
        }

        for(int j = 0; j < NUM_GALAXIES1; j++)
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
            
            if (do_calc)
            {
                dist = distance(h_alpha0[i],h_delta0[i],h_alpha1[j],h_delta1[j],conv_factor_angle);

                if(dist < hist_min)
                    bin_index = 0;
                else if(dist >= hist_max)
                    bin_index = nbins + 1;
                else
                {
                    if (log_binning==0)
                    {
                        bin_index = int((dist-hist_min)/bin_width) + 1;
                    }
                    else if (log_binning==1)// log binning
                    {
                        bin_index = int((log(dist)-log(hist_min))/bin_width) + 1;
                    }
                    else if (log_binning==2)// log 10 binning
                    {
                        bin_index = int((log10(dist)-log10(hist_min))/bin_width) + 1;
                    }
                }


                hist[bin_index]++;
            }
        }
    }  

    unsigned long total = 0;
    float bins_mid = 0;

    float lo = hist_lower_range;
    float hi = 0;
    for(int k=0; k<nbins+1; k++)
    {
        if (k==0)
        {
            //fprintf(outfile, "Underflow below %.3e %s %lu \n", lo, ",",  hist[k]);
        }
        else
        {
            if (log_binning_flag==0)
            {
                hi = lo + hist_bin_width;
            }
            else if (log_binning_flag==1)
            {
                //printf("lo: %f\t\tlog(lo): %f\n",lo,log(lo));
                hi = exp(log(lo) + hist_bin_width);
            }
            else if (log_binning_flag==2)
            {
                //printf("lo: %f\t\tlog10(lo): %f\n",lo,log10(lo));
                hi = pow(10,(log10(lo) + hist_bin_width));
            }

            fprintf(outfile, "%.3e %.3e %lu \n",lo,hi,hist[k]);
            total += hist[k];

            lo = hi;
        }
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
