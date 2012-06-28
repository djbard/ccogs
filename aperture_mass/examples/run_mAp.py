
import sys, time, subprocess

from mAp_tools import *

## Edit this file to include the paths to your input and output files,
## the number of galaxies in your sample and the filter radius 

## input catalogue path. The catalogue must be in the format:
### ra dec reduced_gamma1  reduced_gamma 2
### where ra and dec shoudl be in arcminutes.
input_catname = "example_cat.txt"

## output catalogue path
output_catname = "./out.txt"

## Where ismAp_gals.x" your executable?
exe = "../src/mAp_gals.x"
## critical radius of NFW filter in arcminutes
theta_max = 8

## Do you want to calculate the aperture mass at the position
## of each galaxy, or at a regular grid of positions?
#output_type = "galaxy"
output_type = "grid"


## If a grid, what dimesnsions should it be?
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
if output_type == "galaxy":
    subcmd = [exe, input_catname, output_catname, str(number_of_galaxies), str(theta_max)]
    subprocess.call(subcmd)
    ### now turn the output into a fits file.
    fitsname = "output-gals.fits"
    makeFits_gal(output_catname, fitsname)

    
elif output_type == "grid":
    ## for this, we also need the min and max ra and dec,
    ## so we can get the grid in teh correct coordinates
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
    ### now turn the output into a fits file.
    fitsname = "output-grid.fits"
    makeFits_grid(output_catname, fitsname, grid_size)

else:
    print "You need to specify whether the aperture mass will be evaulated at the positions of the input galaxies or on a grid. Please set output_type. " 
