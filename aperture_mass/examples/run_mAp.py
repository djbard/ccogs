
import sys, time, subprocess

from mAp_tools import *

## Edit this file to include the paths to your input and output files,
## the number of galaxies in your sample and the filter radius 

## input catalogue path. The catalogue must be in the format:
<<<<<<< HEAD
### ra dec reduced_gamma1  reduced_gamma 2
### where ra and dec shoudl be in arcminutes.
input_catname = "example_cat_small.cat"
=======
##   ra dec reduced_gamma1  reduced_gamma 2
## where ra and dec should be in arcminutes.
## There are two example catalogues included,
## one containing 100,000 objects (sample_data/example_cat_small.cat)
## one containing 1,500,000 objects (sample_data/example_cat_large.cat)
input_catname = "../sample_data/example_cat_small.cat"
>>>>>>> 9a2ce61e6fa25835d555f7705d216b8a9e2808ab

## output catalogue path
output_catname = "./out.txt"


## Where is your executable?
exe = "../src/mAp_grid.x"


## critical radius of NFW filter in arcminutes
theta_max = 8


## What dimesnsions should the grid be?
## This should depend on the memory of your GPU card - if
## you find it fails to run due to insufficient memory then
## reduce the grid size! 
grid_size = 1024




## get the # galaxies in your catalogue:
number_of_galaxies = 0
for line in open(input_catname):
    if "#" in line:
        continue
    number_of_galaxies+=1


## run the code. 

## for this, we also need the min and max ra and dec,
## so we can get the grid in the correct coordinates
ra, dec = [], []
for line in open(input_catname):
    cols = line.split()
    ra.append(float(cols[0]))
    dec.append(float(cols[1]))
min_ra = min(ra)
max_ra = max(ra)
min_dec = min(dec)
max_dec = max(dec)


subcmd = [exe, input_catname, output_catname, str(number_of_galaxies), str(theta_max), str(grid_size), str(min_ra), str(max_ra), str(min_dec), str(max_dec)]
subprocess.call(subcmd)


### Uncomment this section if you wish to turn the output into a fits file.
#fitsname = "output-grid.fits"
#makeFits_grid(output_catname, fitsname, grid_size)
