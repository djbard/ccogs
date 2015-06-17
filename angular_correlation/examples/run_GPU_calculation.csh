#!/bin/csh

set BIN_DIR = '../bin/'
set executable = $BIN_DIR/'angular_correlation'

set ngals = 10 # In thousands (10 = 10k)
if ( $1 != '' ) then
    set ngals = $1
endif

#set flat = '../sample_data/flat_'$ngals'k_arcmin.dat'
#set data = '../sample_data/data_'$ngals'k_arcmin.dat'

set data = '~/Downloads/dr10dat.txt'
set flat = '~/Downloads/mockdat.txt'
#set data = '~/Downloads/adr10dat.txt'
#set flat = '~/Downloads/amockdat.txt'
#set data = '~/Downloads/dr10-10k.dat'
#set flat = '~/Downloads/dr10-100k-mock.dat'
#set flat = '~/Downloads/mockdat.txt'
#set flat = '~/Downloads/allrandom.txt'
#set flat = '~/Downloads/abigrandom.txt'

################################################################################
# Read in data assuming arc minutes. (-m)
# Even-spaced binning (-l 0)
# Bin width of 1.0 (-w 1.0)
# Low-edge of 1st bin is 1 arg min. (-L 1.00)
################################################################################
#set global_params = '-w 0.5 -L 0.00 -l 0 -D'
#set tag = 'evenbinning_GPU'

# Mpc
set global_params = '-w 1.0 -L 0.00 -l 0 -p'
#set tag = 'evenbinning_GPU_cartesian_allrandom_width1Mpc'
#set tag = 'evenbinning_GPU_cartesian_abigrandom_width2Mpc'
#set tag = 'evenbinning_GPU_cartesian_amockdat_width2Mpc'
set tag = 'evenbinning_GPU_cartesian_mockdat_width1Mpc'
#set tag = 'evenbinning_GPU_cartesian_DJB_files_10k100k'

################################################################################
# Read in data assuming arc minutes. (-m)
# Log binning (base e) (-l 1)
# Bin width of 0.05 (-w 0.05)
# Low-edge of 1st bin is 1 arg min. (-L 1.00)
################################################################################
#set global_params = '-w 0.05 -L 1.00 -l 1 -m'
#set tag = 'logbinning_GPU'

################################################################################
# Read in data assuming arc minutes. (-m)
# Log10 binning (base 10) (-l 2)
# Bin width of 0.02 (-w 0.02)
# Low-edge of 1st bin is 1 arg min. (-L 1.00)
################################################################################
#set global_params = '-w 0.02 -L 1.00 -l 2 -m'
#set tag = 'log10binning_GPU'


echo "#####################"
time $executable $data $data $global_params -o "$tag"_"$ngals"k_data_data_arcmin.dat 
echo "#####################"
time $executable $flat $flat $global_params -o "$tag"_"$ngals"k_flat_flat_arcmin.dat 
echo "#####################"
time $executable $data $flat $global_params -o "$tag"_"$ngals"k_data_flat_arcmin.dat 

