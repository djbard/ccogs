#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <assert.h>
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>

#include <unistd.h>



#include <cuda.h>
#include <cuda_runtime_api.h>
#include "cutil.h"


///////////////////////////////////// this version calculates the aperture mass at the galaxy positions. 



void checkCUDAerror(const char *msg);

int checkDeviceSpecs(int number_of_galaxies);



__global__ void mApKernel(float* rgamma1, float* rgamma2, float* ra, float* dec, float* mAp_rgamma, float* var_rgamma, float* SN_rgamma,  int tot_gals, float theta_max)
{
  int idx = blockDim.x * blockIdx.x + threadIdx.x;
   
  float this_ra = ra[idx];
  float this_dec = dec[idx];
  
  // want to include the tails outside the halo radius to which our filter is tuned......
  int kernel_radius = 1.5*theta_max;
  float ang = 0.0;
  float xc = 0.15; // could be passed to fn
  float x = 0, Q = 0;
  
  float rgammaMap = 0;
  float rgammaVar=0;
  float radiff, decdiff, dist;

  float npoints = 0;

  for(int i=0; i<tot_gals; i++){

    radiff = (float)this_ra-ra[i]; 
    if(abs(radiff)>kernel_radius) continue;
    decdiff=(float)this_dec-dec[i];
    if(abs(decdiff)>kernel_radius || (radiff==0 && decdiff==0)) continue;
    
    dist = sqrtf(radiff*radiff + decdiff*decdiff);
    if(abs(dist)>kernel_radius) continue;
    
    // have to do something a bit complicated for teh angle - make sure it's getting the correct range. 
    // based on Jan's get_angle(radiff, decdiff)
    if(radiff==0 && decdiff>0) ang = M_PI/2.0;
    else if(radiff==0 && decdiff<0) ang = -1.0 * M_PI/2.0;
    else if(radiff>0) ang = atanf(decdiff/radiff);
    else if(radiff<0 && decdiff>0) ang = atanf(decdiff/radiff) + M_PI;
    else if(radiff<0 and decdiff<0) ang = atanf(decdiff/radiff)-M_PI;
    
    x = dist / theta_max;
    Q = (1.0 / (1.0 + exp(6.0 - 150.0*x) + exp(-47.0 + 50.0*x))) * (tanh(x/xc) / (x/xc));
    
    rgammaMap+=Q* (-1* (rgamma1[i]*cos(2*ang) + rgamma2[i]*sin(2*ang) ));
    rgammaVar+= Q*Q* (rgamma1[i]*rgamma1[i] + rgamma2[i]*rgamma2[i]);
    
    
     npoints++;

}
  
  
  mAp_rgamma[idx] = rgammaMap/npoints;// got to normalise by the # gals I did the sum over. 
  var_rgamma[idx] = rgammaVar /(2*npoints*npoints); 
  SN_rgamma[idx] = sqrtf(2) * rgammaMap / sqrtf(rgammaVar);

}




////////////////////////////////////////////////////////////////////////////////////
//////////  ********** aperture mass for ellip ***************  ////////////////////
////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv){
  

  char* input_filename; char* output_filename;
  int number_of_galaxies; 
  float filter_rad; 
  if (argc>1)
    {
      input_filename = argv[1];
      output_filename = argv[2];
      number_of_galaxies = atoi(argv[3]);
      filter_rad = atof(argv[4]);
    }

 

     
  // CPU memory
  size_t sizeneeded = number_of_galaxies*sizeof(float);
  float *h_rgamma1 = 0;
  float *h_rgamma2 = 0;
  float *h_ra = 0;
  float *h_dec = 0;
  h_rgamma1 = (float*) malloc(sizeneeded);
  h_rgamma2 = (float*) malloc(sizeneeded);
  h_ra = (float*) malloc(sizeneeded);
  h_dec = (float*) malloc(sizeneeded);
  

  ifstream infile;
  infile.open(input_filename);
  
  int i=0;
  float x, y, g1, g2;
  while(1)
    {
      infile>>x>>y>>g1>>g2;
      h_ra[i] = x;
      h_dec[i] = y;
      h_rgamma1[i] = g1;
      h_rgamma2[i] = g2;   
      i += 1;
      if(!infile.good()) break;       

    }
             
    
  // check whether the device has the capacity to do this calculation. 
  int max_threads = checkDeviceSpecs(number_of_galaxies);

  
  /// first, I need to test whether the device is busy. If so, it can wait a little while.
    // of course, even when this loop finished the device may still be busy. So it goes. 
    while(1){
      size_t testsize = 1*sizeof(float); 
      float *d_test;
      cudaMalloc(&d_test, testsize);
      cudaError_t err = cudaGetLastError();
      if( cudaSuccess != err){
	printf("gotta wait for a bit!: %s\n",  cudaGetErrorString( err) );
	sleep(10);
      }
      else break;
    }
    
    
    // GPU memory
    float *d_rgamma1, *d_rgamma2, *d_ra, *d_dec;
    cudaMalloc(&d_rgamma1, sizeneeded);
    cudaMalloc(&d_rgamma2, sizeneeded);
    cudaMalloc(&d_ra, sizeneeded);
    cudaMalloc(&d_dec, sizeneeded);
    
    // output_mAp vector is going to be the calculated value for each point
    size_t sizeneeded_out = number_of_galaxies*sizeof(float);
    float *h_mAp_rgamma,*d_mAp_rgamma, *h_var_rgamma, *d_var_rgamma, *h_SN_rgamma, *d_SN_rgamma;
    
    h_mAp_rgamma = (float*)malloc(sizeneeded_out);
    cudaMalloc(&d_mAp_rgamma, sizeneeded_out);
    h_var_rgamma = (float*)malloc(sizeneeded_out);
    cudaMalloc(&d_var_rgamma, sizeneeded_out);
    h_SN_rgamma = (float*)malloc(sizeneeded_out);
    cudaMalloc(&d_SN_rgamma, sizeneeded_out);
    
    
    //copy vectors from host to device memory
    cudaMemcpy(d_rgamma1, h_rgamma1, sizeneeded, cudaMemcpyHostToDevice);
    cudaMemcpy(d_rgamma2, h_rgamma2, sizeneeded, cudaMemcpyHostToDevice);
    cudaMemcpy(d_ra, h_ra, sizeneeded, cudaMemcpyHostToDevice);
    cudaMemcpy(d_dec, h_dec, sizeneeded, cudaMemcpyHostToDevice);
    cudaMemcpy(d_mAp_rgamma, h_mAp_rgamma, sizeneeded_out, cudaMemcpyHostToDevice);
    cudaMemcpy(d_var_rgamma, h_var_rgamma, sizeneeded_out, cudaMemcpyHostToDevice);
    cudaMemcpy(d_SN_rgamma, h_SN_rgamma, sizeneeded_out, cudaMemcpyHostToDevice);
    
    //check memory is alright
    if (0==h_rgamma1 || 0==h_rgamma2  || 0==h_ra || 0==h_dec || 0==h_mAp_rgamma || 0==h_var_rgamma|| 0==h_SN_rgamma ) printf("can't allocate memory on host \n");
    if (0==d_rgamma1 || 0==d_rgamma2  || 0==d_ra || 0==d_dec  || 0==d_mAp_rgamma || 0==d_var_rgamma || 0==d_SN_rgamma ) printf("can't allocate memory on device \n");
    checkCUDAerror("memory");
    
    
    
    // set up kernel params
    int threadsPerBlock = max_threads; // 
    int blocksPerGrid = int(ceil(number_of_galaxies / float(max_threads)) ); // need nx*nx threads total
    printf(" theads per block: %d and blocks per grid: %d for a total of: %d\n", threadsPerBlock, blocksPerGrid, threadsPerBlock*blocksPerGrid);
    
    
    mApKernel<<<blocksPerGrid, threadsPerBlock >>>(d_rgamma1, d_rgamma2, d_ra, d_dec, d_mAp_rgamma, d_var_rgamma, d_SN_rgamma,  number_of_galaxies, filter_rad); 
    
    checkCUDAerror("kernel");
    
    
    //get the output_mAp back off the device
    cudaMemcpy(h_mAp_rgamma, d_mAp_rgamma, sizeneeded_out, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_var_rgamma, d_var_rgamma, sizeneeded_out, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_SN_rgamma, d_SN_rgamma, sizeneeded_out, cudaMemcpyDeviceToHost);
    
    
    
    // finally, write out to the output file! 
    
   
    FILE *output_file;
    i=0;
    double sq2=sqrt(2.0);
    output_file = fopen(output_filename, "w");
    fprintf(output_file, " # ra  dec  mAp  Var S/N \n");
    for(i=0 ; i<number_of_galaxies; i++){
      fprintf(output_file, "%f %f %f %f %f \n", h_ra[i], h_dec[i], h_mAp_rgamma[i], h_var_rgamma[i],  h_SN_rgamma[i]);
    }
    fclose(output_file);
    
    printf("successfuly completed!\n");
}





//////////////////////////////////////////////////////////////////////////////////////////////////
//simple function to check for errors. From Dr Dobbs. 
void checkCUDAerror(const char *msg)
{
  cudaError_t err = cudaGetLastError();
  if( cudaSuccess != err) 
    {
      fprintf(stderr, "Cuda error: %s: %s.\n", msg, 
	      cudaGetErrorString( err) );
        exit(EXIT_FAILURE);
    }                         
}



//////////////////////////////////////////////////////////////////////////////////

int checkDeviceSpecs( int number_of_galaxies){

  /// what size simulation are we trying to use?

  /// in case we decide to calcuate aperture mass at more points than the number of galaxies...

  int nxny = number_of_galaxies;

  int gpu_mem_needed = int(number_of_galaxies * sizeof(float)) * 7; // need to allocate gamma1, gamma2, ra, dec and output mAp and var and SN. 
  printf("Requirements: %d calculations and %d bytes memory on the GPU \n\n", nxny, gpu_mem_needed);  

  int threadsPerBlock=0;
  // now get the info from the device. 
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
  
  
  int dev, driverVersion = 0, runtimeVersion = 0;     
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
        
    // check memory
    if((unsigned long long) deviceProp.totalGlobalMem < gpu_mem_needed) {
      printf(" FAILURE: Not eneough memeory on device for this calculation! \n");
      exit(1);
    }    
    else
      { 
	printf("Hurrah! This device has enough memory to perform this calculation\n");
	
	// check # threads
	
	threadsPerBlock = deviceProp.maxThreadsPerBlock; // maximal efficiency exists if we use max # threads per block. 
	int blocksPerGrid = int(ceil(nxny / threadsPerBlock)); // need nx*nx threads total

	if( int(deviceProp.maxThreadsDim[1])*int(deviceProp.maxThreadsDim[2]) <blocksPerGrid) {
	  printf("FAILURE: Not enough threads on teh device to do this calculation!\n");
	    exit(1);
	  }
	else 
	  {
	    printf("Hurrah! This device supports enough threads to do this calculation\n");
	  }
      }

  }// loop over devices
  
  return threadsPerBlock;
}
