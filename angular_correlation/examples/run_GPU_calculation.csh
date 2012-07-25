#!/bin/csh

set BIN_DIR = '../bin/'
set executable = $BIN_DIR/'angular_correlation'
set tag = 'GPU'

set ngals = 10 # In thousands (10 = 10k)
if ( $1 != '' ) then
    set ngals = $1
endif

set flat = '../sample_data/flat_'$ngals'k_arcmin.dat'
set data = '../sample_data/data_'$ngals'k_arcmin.dat'

################################################################################
# Read in data assuming arc minutes. (-m)
# Log binning (base e) (-l 1)
# Bin width of 0.05 (-w 0.05)
# Low-edge of 1st bin is 1 arg min. (-L 1.00)
################################################################################

set global_params = '-w 0.05 -L 1.00 -l 1 -m'

echo "\n#####################\n"
time $executable $data $data $global_params -o logbinning_"$tag"_"$ngals"k_data_data_arcmin.dat 
echo "\n#####################\n"
time $executable $flat $flat $global_params -o logbinning_"$tag"_"$ngals"k_flat_flat_arcmin.dat 
echo "\n#####################\n"
time $executable $data $flat $global_params -o logbinning_"$tag"_"$ngals"k_data_flat_arcmin.dat 

