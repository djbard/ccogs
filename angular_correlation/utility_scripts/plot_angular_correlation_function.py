import numpy as np
import matplotlib.pyplot as plt

################################################################################
# main
################################################################################
def main():

    ############################################################################
    # Open the files and read the data.
    # 
    # You will have to manually edit the filename and number of galaxies
    # in the datasets.
    ############################################################################
    
    ngalaxies = 10000.0

    filenames = [None,None,None]
    filenames[0] = "logbinning_GPU_10k_data_data_arcmin.dat" # DD
    filenames[1] = "logbinning_GPU_10k_flat_flat_arcmin.dat" # RR
    filenames[2] = "logbinning_GPU_10k_data_flat_arcmin.dat" # DR

    ############################################################################
    ############################################################################

    dd = None
    rr = None
    dr = None
    bin_lo = None
    bin_hi = None
    for i,name in enumerate(filenames):
        infile = open(name)
        content = np.array(infile.read().split()).astype('float')
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

    ############################################################################

    # Calculate the normalization.
    dd_norm = ((ngalaxies*ngalaxies)-ngalaxies)/2.0
    rr_norm = ((ngalaxies*ngalaxies)-ngalaxies)/2.0
    dr_norm = (ngalaxies*ngalaxies)

    print dd_norm 
    print rr_norm 
    print dr_norm 

    dd /= dd_norm
    rr /= rr_norm
    dr /= dr_norm
     
    w = (dd-(2*dr)+rr)/rr

    bin_mid = (bin_hi+bin_lo)/2.0
    bin_width = (bin_hi-bin_lo)

    # Divide out the bin width.
    w /= bin_width

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
    ax0.set_xlabel(r"$\theta$ (arcmin)", fontsize=24, weight='bold')
    ax0.set_ylabel(r"w($\theta$)", fontsize=24, weight='bold')
    plt.xticks(fontsize=24,weight='bold')
    plt.yticks(fontsize=24,weight='bold')

    ax0.scatter(bin_mid,w,s=30)
    ax0.set_xlabel(r"$\theta$ (arcminutes)",fontsize=24, weight='bold')
    ax0.set_ylabel(r"w($\theta$)",fontsize=24, weight='bold')

    plt.xticks(fontsize=24,weight='bold')
    plt.yticks(fontsize=24,weight='bold')

    #ax0.set_xscale('log')
    #ax0.set_yscale('log')
   
    #ax0.set_xlim(-100,5000)
    ax0.set_xlim(-10,130)
    #ax0.set_ylim(-0.7,2.8)
    ax0.set_ylim(0.01,100)

    plt.show()

################################################################################
# Top-level script evironment
################################################################################
if __name__ == "__main__":
    main()

