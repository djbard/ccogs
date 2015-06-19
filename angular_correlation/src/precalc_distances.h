#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<unistd.h>

void define_ranges(float *vec_ranges_lo, float *vec_ranges_hi, float maxsep, int *vec_nbins, float *vec_binwidths);

void assign_chunk_coordinate(float *coord, float *vec_ranges_lo, float *vec_binwidths, int *chunk_coord);

float distance(float x0, float y0, float z0, float x1, float y1, float z1);
