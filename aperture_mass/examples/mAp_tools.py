import math
import numpy as np


## create a fits file for the case where the aperture mass
## has been calculated at the galaxy positions
def makeFits_gal(incat, outfile):

    ra, dec, sn=[], [], []
    for line in open(incat ,"r"):
        cols = line.split()
        if "#" in line:
            continue
        ra.append(float(cols[0]))
        dec.append(float(cols[1]))
        sn.append(float(cols[4]))


    min_ra = min(ra)
    max_ra = max(ra)
    min_dec = min(dec)
    max_dec = max(dec)

    rough_size = math.sqrt(len(ra))

    ### The galaxies are proably close to randomly distributed
    ### but to preserve any clustering effects we put the
    ### number of pixels in the output fits file to 2*# galaxies.
    ### Feel free to change this if you want.
    
    side = math.floor(2*rough_size)
    arr = np.empty( (side, side) )
    pix_size = (max_ra - min_ra) /side 

    for i in range(0, len(ra)):
        this_ra = math.floor( (ra[i]-min_ra) / pix_size)
        this_dec = math.floor( (dec[i]-min_dec) / pix_size)
        if ra[i]==max_ra or dec[i]==max_dec:
            continue
        arr[this_ra][this_dec] = sn[i]
        
    out_hdu = pyfits.PrimaryHDU(arr)
    out_hdu.writeto(outfile)

    return 1


## create a fits file for the case where the aperture mass
## has been calculated on a grid
## The ouput fits file has the same coords as the input file
## only the points are on a grid. 
def makeFits_grid(incat, outfile, grid_size):

    ra, dec, sn=[], [], []
    for line in open(incat ,"r"):
        cols = line.split()
        if "#" in line:
            continue
        ra.append(float(cols[0]))
        dec.append(float(cols[1]))
        sn.append(float(cols[4]))

    min_ra = min(ra)
    max_ra = max(ra)
    min_dec = min(dec)
    max_dec = max(dec)

    ra_pixsize = (max_ra - min_ra)/float(grid_size);
    dec_pixsize = (max_dec - min_dec)/float(grid_size);

    min_ra = min(ra)
    max_ra = max(ra)
    min_dec = min(dec)
    max_dec = max(dec)
    arr = np.empty( (grid_size, grid_size) )
    print len(ra), grid_size
    for i in range(0, len(ra)):
        tempra = math.floor(i/grid_size);
        tempdec = math.floor(i - grid_size*tempra);

        if ra[i]==max_ra or dec[i]==max_dec:
            continue
        
        arr[tempra][tempdec] = sn[i]
        
        
    #out_hdu = pyfits.PrimaryHDU(arr)
    #out_hdu.writeto(outfile)

    return 1

