#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include <unistd.h>

#include <cuda_runtime.h>

using namespace std;

#define SUBMATRIX_SIZE 16384

////////////////////////////////////////////////////////////////////////
// Number of histogram bins has to be edited by hand, prior to
// copmilation.
////////////////////////////////////////////////////////////////////////

#define DEFAULT_NBINS 254 
//#define DEFAULT_NBINS 126 
//#define DEFAULT_NBINS 62 
//#define DEFAULT_NBINS 30 

#define CONV_FACTOR 57.2957795 // 180/pi

////////////////////////////////////////////////////////////////////////
// Kernel to calculate angular distances between galaxies and histogram
// the distances.
////////////////////////////////////////////////////////////////////////
__global__ void distance(volatile float *a0, volatile float *d0, volatile float *a1, volatile float *d1, int xind, int yind, int max_xind, int max_yind, volatile int *dev_hist, float hist_min, float hist_max, int nbins, float bin_width, int log_binning=0, bool two_different_files=1, float conv_factor_angle=57.2957795)
{

    ////////////////////////////////////////////////////////////////////////////
    // Idx will keep track of which thread is being calculated within a given 
    // warp.
    ////////////////////////////////////////////////////////////////////////////
    int idx = blockIdx.x * blockDim.x + threadIdx.x; // This should range to SUBMATRIX_SIZE

    idx += xind;

    ////////////////////////////////////////////////////////////////////////
    // Shared memory stuff.
    ////////////////////////////////////////////////////////////////////////
    __shared__ int shared_hist[DEFAULT_NBINS+2];
    // Note that we only clear things out for the first thread on each block.
    if(threadIdx.x==0)
    {
        for (int i=0;i<nbins+2;i++)
            shared_hist[i] = 0;
    }
    __syncthreads();
    ////////////////////////////////////////////////////////////////////////

    if (idx<max_xind)
    {
        int i=0;

        float alpha_rad = a0[idx];
        float delta0 = d0[idx];
        float cos_d0 = cos(delta0);
        float sin_d0 = sin(delta0);
        float dist;

        int bin_index = 0; 

        float a_diff, sin_a_diff, cos_a_diff;
        float cos_d1, sin_d1, numer, denom, mult1, mult2;    
        float d1_rad;

        bool do_calc = 1;

        int ymax = yind + SUBMATRIX_SIZE;

        if (ymax>max_yind)
        {
            ymax = max_yind;
        }

        for(i=yind; i<ymax; i++)
        {
            if (two_different_files)
            {
                do_calc = 1;
            }
            else // Doing the same file
            {
                if(idx > i)
                    do_calc=1;
                else
                    do_calc=0;
            }
            //if(idx > i) ///////// CHECK THIS
            if (do_calc)
            {
                a_diff = a1[i] - alpha_rad;
                d1_rad = d1[i];

                sin_a_diff = sin(a_diff);
                cos_a_diff = cos(a_diff);

                sin_d1 = sin(d1_rad);
                cos_d1 = cos(d1_rad);

                mult1 = cos_d1 * cos_d1 * sin_a_diff * sin_a_diff;
                mult2 = cos_d0 * sin_d1 - sin_d0 * cos_d1 * cos_a_diff;
                mult2 = mult2 * mult2;

                numer = sqrt(mult1 + mult2); 

                denom = sin_d0 *sin_d1 + cos_d0 * cos_d1 * cos_a_diff;

                dist = atan2(numer,denom);  
                dist *= conv_factor_angle;  // Convert to degrees or what have you.

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

                atomicAdd(&shared_hist[bin_index],1);

            }
        }
    }

    __syncthreads();

    if(threadIdx.x==0)
    {
        for(int i=0;i<nbins+2;i++)
            dev_hist[i+(blockIdx.x*(nbins+2))]=shared_hist[i];
    }

}

////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
    // Needed for parsing command-line arguments.
    extern char *optarg;
    extern int optind, optopt, opterr;
    int c;
    char *outfilename = NULL;
    char defaultoutfilename[256];
    sprintf(defaultoutfilename,"default_out.dat");

    float hist_lower_range = 0.0000001;
    float hist_upper_range = 0;
    int nbins = DEFAULT_NBINS;
    float hist_bin_width = 0.05;
    int log_binning_flag = 0; // False

    float scale_factor = 1.0; // For if we need to convert input to arcsec or arcmin
    float conv_factor_angle = 57.2957795; // 180/pi // For if we need to convert arcdistance to arcsec or arcmin

    bool silent_on_GPU_testing = false;

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////

    while ((c = getopt(argc, argv, "ao:L:l:w:smS")) != -1) {
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
            case 'S':
                printf("Silent mode - don't run the GPU test (suppresses some output)");
                silent_on_GPU_testing = true;
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

    FILE *infile0, *infile1, *outfile ;
    infile0 = fopen(argv[optind],"r");
    infile1 = fopen(argv[optind+1],"r");

    printf("Opening input file 0: %s\n",argv[optind]);
    printf("Opening input file 1: %s\n",argv[optind+1]);
    outfile = fopen(outfilename, "w");

    ////////////////////////////////////////////////////////////////////////////
    // Check to see if the two files are actually the same file.
    // This is the case for the DD and RR calculations and change slightly
    // the exact calculations being performed.
    ////////////////////////////////////////////////////////////////////////////
    bool two_different_files = 1;
    if (strcmp(argv[optind],argv[optind+1])==0)
    {
        two_different_files = 0;
        printf("Using the same file!\n");
    }
    printf("\n");
    ////////////////////////////////////////////////////////////////////////////
    // Now get the info from the device.
    ////////////////////////////////////////////////////////////////////////////
    if (!silent_on_GPU_testing)
    {
        printf("\n------ CUDA device diagnostics ------\n\n");

        int tot_gals = 100000;
        int nx = SUBMATRIX_SIZE;
        int ncalc = nx * nx;
        int gpu_mem_needed = int(tot_gals * sizeof(float)) * 2; // need to allocate ra, dec.
        printf("Requirements: %d calculations and %d bytes memory on the GPU \n\n", ncalc, gpu_mem_needed);

        int deviceCount = 0;
        cudaError_t error_id = cudaGetDeviceCount(&deviceCount);
        if (error_id != cudaSuccess) {
            printf( "cudaGetDeviceCount returned %d\n-> %s\n", (int)error_id, cudaGetErrorString(error_id) );
        }
        // This function call returns 0 if there are no CUDA capable devices.
        if (deviceCount == 0)
            printf("There is no device supporting CUDA\n");
        else
            printf("Found %d CUDA Capable device(s)\n", deviceCount);


        int dev=0;
        for (dev = 0; dev < deviceCount; ++dev) {
            cudaDeviceProp deviceProp;
            cudaGetDeviceProperties(&deviceProp, dev);
            printf("\nDevice %d: \"%s\"\n", dev, deviceProp.name);

            printf("  Total amount of global memory:                 %.0f MBytes (%llu bytes)\n",
                    (float)deviceProp.totalGlobalMem/1048576.0f, (unsigned long long) deviceProp.totalGlobalMem);


            printf("  Warp size:                                     %d\n", deviceProp.warpSize);
            printf("  Maximum number of threads per block:           %d\n", deviceProp.maxThreadsPerBlock);
            printf("  Maximum sizes of each dimension of a block:    %d x %d x %d\n",
                    deviceProp.maxThreadsDim[0],
                    deviceProp.maxThreadsDim[1],
                    deviceProp.maxThreadsDim[2]);
            printf("  Maximum sizes of each dimension of a grid:     %d x %d x %d\n",
                    deviceProp.maxGridSize[0],
                    deviceProp.maxGridSize[1],
                    deviceProp.maxGridSize[2]);

            // does this device have enough capcacity for the calculation?
            printf("\n*************\n");

            // check memory
            if((unsigned long long) deviceProp.totalGlobalMem < gpu_mem_needed) printf(" FAILURE: Not eneough memeory on device for this calculation! \n");
            else
            {
                printf("Hurrah! This device has enough memory to perform this calculation\n");

                // check # threads

                int threadsPerBlock = deviceProp.maxThreadsPerBlock; // maximal efficiency exists if we use max # threads per block.
                int blocksPerGrid = int(ceil(ncalc / threadsPerBlock)); // need nx*nx threads total
                if(deviceProp.maxThreadsDim[0] >blocksPerGrid) printf("FAILURE: Not enough threads on the device to do this calculation!\n");
                else
                {
                    printf("Hurrah! This device supports enough threads to do this calculation\n");
                    // how many kernels can we run at once on this machine?
                    int n_mem = floor(deviceProp.totalGlobalMem / float(gpu_mem_needed));
                    int n_threads = floor(threadsPerBlock * deviceProp.maxThreadsDim[0]*deviceProp.maxThreadsDim[1] / float(ncalc) ); // max # threads possible?

                    printf("%d %d  \n",  n_threads, deviceProp.maxThreadsDim[0]);

                    int max_kernels = 0;
                    n_mem<n_threads ? max_kernels = n_mem : max_kernels = n_threads;

                    printf(" you can run %d kernels at a time on this device without overloading the resources \n", max_kernels);
                }
            }

        }

        printf("\n------ End CUDA device diagnostics ------\n\n");
    }
    ////////////////////////////////////////////////////////////////////////////

    float *d_alpha0, *d_delta0;
    float *h_alpha0, *h_delta0;

    float *d_alpha1, *d_delta1;
    float *h_alpha1, *h_delta1;

    int NUM_GALAXIES0;
    int NUM_GALAXIES1;

    //////////////////////////////////////////////////////////////////////
    // Read in the galaxy files.
    ////////////////////////////////////////////////////////////////////////////
    // Read in the first file
    ////////////////////////////////////////////////////////////////////////////

    fscanf(infile0, "%d", &NUM_GALAXIES0);

    int size_of_galaxy_array0 = NUM_GALAXIES0 * sizeof(float);    
    printf("SIZE 0 # GALAXIES: %d\n",NUM_GALAXIES0);

    h_alpha0 = (float*)malloc(size_of_galaxy_array0);
    h_delta0 = (float*)malloc(size_of_galaxy_array0);
    float temp0, temp1;

    for(int i=0; i<NUM_GALAXIES0; i++)
    {
        fscanf(infile0, "%f %f", &temp0, &temp1);
        h_alpha0[i] = temp0/scale_factor;
        h_delta0[i] = temp1/scale_factor;
        //if (i<10)
        //printf("%e %e\n", h_alpha0[i], h_delta0[i]);
    }

    ////////////////////////////////////////////////////////////////////////////
    // Read in the second file
    ////////////////////////////////////////////////////////////////////////////

    fscanf(infile1, "%d", &NUM_GALAXIES1);

    int size_of_galaxy_array1 = NUM_GALAXIES1 * sizeof(float);    
    printf("SIZE 1 # GALAXIES: %d\n",NUM_GALAXIES1);

    h_alpha1 = (float*)malloc(size_of_galaxy_array1);
    h_delta1 = (float*)malloc(size_of_galaxy_array1);

    for(int i=0; i<NUM_GALAXIES1; i++)
    {
        fscanf(infile1, "%f %f", &temp0, &temp1);
        h_alpha1[i] = temp0/scale_factor;
        h_delta1[i] = temp1/scale_factor;
        //if (i<10)
        //printf("%e %e\n", h_alpha1[i], h_delta1[i]);
    }

    ////////////////////////////////////////////////////////////////////////////
    // Allocation of histogram
    ///////////////////////////////////////////////////////////////////////////

    int *hist, *dev_hist;

    int size_hist = SUBMATRIX_SIZE * (nbins+2);
    int size_hist_bytes = size_hist*sizeof(int);

    hist = (int*)malloc(size_hist_bytes);
    memset(hist, 0, size_hist_bytes);

    printf("Size of histogram: %d bytes\n",size_hist_bytes);
    cudaMalloc((void **) &dev_hist, (size_hist_bytes));
    cudaMemset(dev_hist, 0, size_hist_bytes);

    unsigned long  *hist_array;

    int hist_array_size = (nbins+2) * sizeof(unsigned long);
    hist_array =  (unsigned long*)malloc(hist_array_size);
    printf("Size of histogram array: %d bytes\n",hist_array_size);
    memset(hist_array,0,hist_array_size); 

    ////////////////////////////////////////////////////////////////////////////
    // Define the grid and block size
    ////////////////////////////////////////////////////////////////////////////
    dim3 grid, block;
    // 128*4 = 512, the amount of memory needed for one histogram.
    // 8192*4 = 32768 is max memory to ask for for the histograms.
    // 8192/128 = 64, is is the right number of blocks?
    grid.x = 8192/(DEFAULT_NBINS+2); // Is this the number of blocks?
    block.x = SUBMATRIX_SIZE/grid.x; // Is this the number of threads per block? NUM_GALAXIES/block.x;
    // SUBMATRIX is the number of threads per warp? Per kernel call?
    ////////////////////////////////////////////////////////////////////////////

    cudaMalloc((void **) &d_alpha0, size_of_galaxy_array0 );
    cudaMalloc((void **) &d_delta0, size_of_galaxy_array0 );

    cudaMalloc((void **) &d_alpha1, size_of_galaxy_array1 );
    cudaMalloc((void **) &d_delta1, size_of_galaxy_array1 );

    // Check to see if we allocated enough memory.
    if (0==d_alpha0 || 0==d_delta0 || 0==d_alpha1 || 0==d_delta1 || 0==dev_hist)
    {
        printf("couldn't allocate memory\n");
        return 1;
    }

    // Initialize array to all 0's
    cudaMemset(d_alpha0,0,size_of_galaxy_array0);
    cudaMemset(d_delta0,0,size_of_galaxy_array0);
    cudaMemset(d_alpha1,0,size_of_galaxy_array1);
    cudaMemset(d_delta1,0,size_of_galaxy_array1);

    cudaMemcpy(d_alpha0, h_alpha0, size_of_galaxy_array0, cudaMemcpyHostToDevice );
    cudaMemcpy(d_delta0, h_delta0, size_of_galaxy_array0, cudaMemcpyHostToDevice );
    cudaMemcpy(d_alpha1, h_alpha1, size_of_galaxy_array1, cudaMemcpyHostToDevice );
    cudaMemcpy(d_delta1, h_delta1, size_of_galaxy_array1, cudaMemcpyHostToDevice );

    int x, y;
    int num_submatrices_x = NUM_GALAXIES0 / SUBMATRIX_SIZE;
    int num_submatrices_y = NUM_GALAXIES1 / SUBMATRIX_SIZE;
    // Take care of edges of matrix.
    if (NUM_GALAXIES0%SUBMATRIX_SIZE != 0)
    {
        num_submatrices_x += 1;
    }
    if (NUM_GALAXIES1%SUBMATRIX_SIZE != 0)
    {
        num_submatrices_y += 1;
    }


    printf("Breaking down the calculations.\n");
    printf("Number of submatrices: %dx%d\n",num_submatrices_x,num_submatrices_y);
    printf("Number of calculations per submatrices: %dx%d\n",SUBMATRIX_SIZE,SUBMATRIX_SIZE);

    int bin_index = 0;
    for(int k = 0; k < num_submatrices_y; k++)
    {
        y = k*SUBMATRIX_SIZE;
        //printf("%d %d\n",k,y);
        for(int j = 0; j < num_submatrices_x; j++)
        {
            x = j*SUBMATRIX_SIZE; 

            //printf("----\n");
            //printf("%d %d\t\t%d %d\n",k,y,j,x);
            //printf("----\n");

            // Set the histogram to all zeros each time.
            cudaMemset(dev_hist,0,size_hist_bytes);

            int max_x = NUM_GALAXIES0;
            int max_y = NUM_GALAXIES1;

            distance<<<grid,block>>>(d_alpha0, d_delta0,d_alpha1, d_delta1, x, y, max_x, max_y, dev_hist, hist_lower_range, hist_upper_range, nbins, hist_bin_width, log_binning_flag, two_different_files,conv_factor_angle);
            cudaMemcpy(hist, dev_hist, size_hist_bytes, cudaMemcpyDeviceToHost);

            ////////////////////////////////////////////////////////////////////
            // Sum up the histograms from each thread (hist).
            ////////////////////////////////////////////////////////////////////
            for(int m=0; m<size_hist; m++)
            {
                bin_index = m%(nbins+2);
                hist_array[bin_index] += hist[m];
            }    
        }  
    }

    unsigned long total = 0;

    float lo = hist_lower_range;
    float hi = 0;
    for(int k=0; k<nbins+1; k++)
    {
        if (k==0)
        {
            //fprintf(outfile, "Underflow below %.3e %s %lu \n", lo, ",",  hist_array[k]);
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

            fprintf(outfile, "%.3e %.3e %lu \n",lo,hi,hist_array[k]);
            total += hist_array[k];

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

    cudaFree(d_alpha0);
    cudaFree(d_delta0);  
    cudaFree(d_alpha1);
    cudaFree(d_delta1);  
    cudaFree(dev_hist);

    return 0;
}  
//////////////////////////////////////////////////////////////////////
