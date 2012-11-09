import os, math, sys
### If you're not using a fits table, and you don't have pyfits, you shoudl comment out this import. 
import pyfits 
import numpy as np
from ROOT import *

###############################################
### Fit the redshift - Mpc relation taken from CosmoCalc reference points. 
def getRedshiftMpcFit():
    ### first, I'm gonna plot my cosmo-dependent redshift/distance relation
    redshift_model = [.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
    Mpc_model = [209.0, 413.4, 613.5, 808.8, 999.4, 1185.3, 1366.5, 1542.9, 1714.6, 1881.7]

    gr_model = TGraph(len(redshift_model), np.array(redshift_model), np.array(Mpc_model))
    gr_model.Fit("pol2")
    fitresults = gr_model.GetFunction('pol2')

    p0 = fitresults.GetParameter(0)
    p1 = fitresults.GetParameter(1)
    p2 = fitresults.GetParameter(2)

    print "polynomial fit params: p0 =", p0, ", p1 =", p1, ", p2 =", p2
    
    return p0, p1, p2




###############################################
### convert redshift z into Mpc
def getMpcFromRedshift(redshift, p0, p1, p2):
    Mpc = []
    for i in range(0, len(redshift)):
        Mpc.append(p0 + p1*redshift[i] + p2*redshift[i]*redshift[i])
    
    return Mpc


###############################################
### Convert ra/dec/Mpc coords into x/y/z coords. 
def convertRaDecMpcToXYZ(ra, dec, Mpc):
    x, y, z, = [], [], []
    rad = math.pi/180.0
    
    for i in range(0, len(ra)):
        x.append(Mpc[i]*sin(rad*(-1.0*dec[i]+90))*cos(rad*(ra[i])))
        y.append(Mpc[i]*sin(rad*(-1.0*dec[i]+90))*sin(rad*(ra[i])))
        z.append(Mpc[i]*cos(rad*(-1.0*dec[i]+90)))
                 
            
    
    return x, y, z

###############################################
### read info out of the fits file and into arrays. 
### You'll want to edit this to match the data in your fits table
def readFitsFile(hdulist):
    data = hdulist[1].data
    
    ra = data.field('RA')
    dec = data.field('DEC')
    redshift = data.field('ZRED')
    return ra, dec, redshift


###############################################
### main bit of code!

### input file:
inputfile = sys.argv[1]

### Fit 2nd order polynomial to Redshift-Mpc relation (based on standard cosmology) and get fit function params. 
p0, p1, p2 = getRedshiftMpcFit()


### Get ra/dec/redshift from the input files
### Is it a fits table? 
if ".fit" in inputfile:
    gallist = pyfits.open(inputfile)
    gal_ra, gal_dec, gal_redshift = readFitsFile(gallist)
### Is it a .txt file?  Assume format is: ra dec z
else:
    gal_ra, gal_Dec, gal_z = [], [], []
    for line in open(inputfile, 'r'):
        cols = line.split()
        if len(cols)<3:
            continue
        gal_ra.append(float(cols[0]))
        gal_dec.append(float(cols[0]))
        gal_z.append(float(cols[0]))
        

### convert redshift to Mpc
gal_Mpc = getMpcFromRedshift(gal_redshift, p0, p1, p2)


### convert ra,dec, Mpc to x,y,z in Mpc
gal_x, gal_y, gal_z = convertRaDecMpcToXYZ(gal_ra, gal_dec, gal_Mpc)


### Now, print out the files.
gal_datfile = open("default_xyz.dat", "w")
## first line is # gals in the file.
ngals = 0
for i in range(0, len(gal_ra)):
    ngals+=1
    
print >> gal_datfile, ngals

for i in range(0, len(gal_x)):
    print >> gal_datfile, gal_x[i], gal_y[i], gal_z[i]
gal_datfile.close()

