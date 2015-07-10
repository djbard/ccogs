import sys
import numpy as np
import matplotlib.pyplot as plt

################################################################################
# main
################################################################################
def main():

    ############################################################################
    # Open the files and read the data.
    # 
    # You may prefer to manually edit the filename and number of galaxies
    # in the datasets.
    ############################################################################
    
    # Read in off the command line some string to look for in the 
    # input files.
    '''
    tag = 'logbinning_GPU_100k'
    if len(sys.argv)>=2:
        tag = sys.argv[1]
    '''

    filenames = [None,None,None]
    filenames[0] = sys.argv[1] # DD
    filenames[1] = filenames[0].replace("data","flat") # RR
    filenames[2] = filenames[0].replace("data_data","data_flat") # DR

    tag = filenames[0].split("_data")[0]

    # Pull the number of galaxies out of the file name.
    ngal_in_file = tag.split('_')[-1][0:-1]
    #ngalaxies = float(ngal_in_file)*1000.0
    #nd = 401926.
    #nr = 213277.
    #nr = 1075634.
    #nr = 1075634.
    #nr = 4299990.
    #nr = 5375624.
    #nr = 216323.

    #nd = 10000.
    #nr = 100000.

    #nd = 401927.
    #nr = 216323.

    #nd = 401927.
    #nr = 213277.

    nd = 101292.
    nr = 57997.

    ############################################################################
    ############################################################################

    dd = None
    rr = None
    dr = None
    bin_lo = None
    bin_hi = None

    # Loop over the files and pull out the necessary info.
    for i,name in enumerate(filenames):

        print "Opening: ",name 
        infile = open(name)

        # Parse the entire contents of the file into a big array of floats.
        content = np.array(infile.read().split()).astype('float')

        # We know there are three columns of numbers, so we can pull out what
        # we want using an array of the indices.
        nentries = len(content)
        nbins = nentries/3
        index = np.arange(0,nentries,3)

        if i==0:
            bin_lo = content[index]
            bin_hi = content[index+1]
            dd = content[index+2]
        elif i==1:
            rr = content[index+2]
        elif i==2:
            dr = content[index+2]

    print dd
    ############################################################################

    # Calculate the normalization.
    dd_norm = ((nd*nd)-nd)/2.0
    rr_norm = ((nr*nr)-nr)/2.0
    dr_norm = (nd*nr)/1.0

    print "DD normalization:",dd_norm 
    print "RR normalization:",rr_norm 
    print "DR normalization:",dr_norm 

    # Normalize the data appropriately.
    dd /= dd_norm
    rr /= rr_norm
    dr /= dr_norm
     
    # Calculate the angular correlation function here.
    w = (dd-(2.0*dr)+rr)/rr

    bin_mid = (bin_hi+bin_lo)/2.0
    bin_width = (bin_hi-bin_lo)

    # Divide out the bin width.
    #w /= bin_width

    ############################################################################
    # Write out the function to a file.
    ############################################################################
    outfile = open('default_acf.dat','w+')
    for lo,hi,wval in zip(bin_lo,bin_hi,w):
        if wval==wval: # Check for nans and infs
            output = "%.3e %.3e %f\n" % (lo,hi,wval)
            outfile.write(output)
    outfile.close()
    ############################################################################

    ################################################################################
    # Make a figure on which to plot the angular correlation function.
    ################################################################################
    fig0 = plt.figure(figsize=(9,6),dpi=100,facecolor='w',edgecolor='k')
    ax0 = fig0.add_subplot(1,1,1)
    fig0.subplots_adjust(top=0.95,bottom=0.15,right=0.95)
    ################################################################################
    
    ############################################################################
    # Format the plot.
    ############################################################################
    ax0.set_xlabel(r"$r$ (Mpc)", fontsize=24, weight='bold')
    ax0.set_ylabel(r"w($r$)", fontsize=24, weight='bold')
    plt.xticks(fontsize=24,weight='bold')
    plt.yticks(fontsize=24,weight='bold')

    #bin_mid *= 0.7

    #ax0.scatter(bin_mid,bin_mid*w,s=30)
    #ax0.scatter(bin_mid,np.log10(w),s=30)
    ax0.scatter(bin_mid,w,s=30)
    #ax0.set_xlabel(r"$w$ (r)",fontsize=24, weight='bold')
    #ax0.set_ylabel(r"w($\theta$)",fontsize=24, weight='bold')

    plt.xticks(fontsize=24,weight='bold')
    plt.yticks(fontsize=24,weight='bold')

    #ax0.set_xscale('log')
    #ax0.set_yscale('log')
   
    ax0.set_xlim(10,200)
    #ax0.set_xlim(-10,130)
    #ax0.set_ylim(0.001,0.4)
    #ax0.set_ylim(0.01,100)
    #ax0.set_ylim(0.01,5)

    ################################################################################
    # Make a figure on which to plot the DD,RR,DR
    ################################################################################
    '''
    fig1 = plt.figure(figsize=(15,4),dpi=100,facecolor='w',edgecolor='k')
    fig1.add_subplot(1,3,1)
    plt.scatter(bin_mid,dd,s=30)
    plt.xlabel(r"$\theta$ (degrees)", fontsize=24, weight='bold')
    plt.ylabel(r"DD normalized counts", fontsize=18, weight='bold')

    fig1.add_subplot(1,3,2)
    plt.scatter(bin_mid,rr,s=30)
    plt.xlabel(r"$\theta$ (degrees)", fontsize=24, weight='bold')
    plt.ylabel(r"RR normalized counts", fontsize=18, weight='bold')

    fig1.add_subplot(1,3,3)
    plt.scatter(bin_mid,dr,s=30)
    plt.xlabel(r"$\theta$ (degrees)", fontsize=24, weight='bold')
    plt.ylabel(r"DR normalized counts", fontsize=18, weight='bold')
    '''

    plt.tight_layout()


    plt.show()

################################################################################
# Top-level script evironment
################################################################################
if __name__ == "__main__":
    main()

