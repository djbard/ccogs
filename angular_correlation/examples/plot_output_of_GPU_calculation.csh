#!/usr/bin/env csh

set pyscript = '../utility_scripts/plot_angular_correlation_function.py'

set tag = 'evenbinning_GPU_10k'
#set tag = 'evenbinning_GPU_100k'
#set tag = 'log10binning_GPU_10k'
#set tag = 'log10binning_GPU_100k'
#set tag = 'logbinning_GPU_10k'
#set tag = 'logbinning_GPU_100k'

python $pyscript $tag
