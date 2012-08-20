# Makefile
#
# author      : Deborah Bard
# date        : Jan., 2012
#

################################################################################
# This section contains some environment variables that you may want to edit, or
# perhaps need to edit, depending on your CUDA installation.
################################################################################

# INSTALL_DIR is where the CUDA and C executables will be installed. You may find it 
# useful to install in some system-wide directory (e.g. /usr/local/bin, 
# $(HOME)/bin, etc.).
INSTALL_DIR       := ./

BIN_GRID           := $(INSTALL_DIR)/mAp_grid.x


# SDK_INSTALL_PATH and CUDA_INSTALL_PATH are probably the more significant
# environment variables and will depend on your system and where the CUDA
# libraries and header files are installed.

### You will need to edit this line to point to your own SDK installation
CUDA_SDK = /path/to/my/cuda/sdk/installation

### You will need to edit this line to point to your own cuda installation
CUDA_INSTALL_PATH := /opt/cuda




# While the INCLUDES and LIBS flags are dependent on the install path, there may 
# still be some subdirectories that are machine-dependent.

INCLUDES          += -I. -I$(CUDA_INSTALL_PATH)/include  -I$(CUDA_SDK)/C/common/inc 
LIBS		  := -L$(CUDA_INSTALL_PATH)/lib64  -L$(CUDA_SDK)/C/common/lib/ -L$(CUDA_SDK)/C/lib/ 
COMMONFLAGS       := -m64
CXXFLAGS          := $(COMMONFLAGS)
LDFLAGS           := -lm -lcudart

# compilers
NVCC              := $(CUDA_INSTALL_PATH)/bin/nvcc -arch sm_11 $(COMMONFLAGS)


CU_SOURCES_GRID   := mAp_gpu_grid.cu
CU_OBJS_GRID      := $(patsubst %.cu, %.cu_o, $(CU_SOURCES_GRID))


###############################################################################

################################################################################
# Rules for Makefile
################################################################################


%.cu_o : %.cu 
	$(NVCC) -c $(INCLUDES) -o $@ $^
all:  ${BIN_GRID}


$(BIN_GRID): $(CPP_OBJS) $(CU_OBJS_GRID)
	$(CXX) -fPIC -o $(BIN_GRID) $(CU_OBJS_GRID) $(LDFLAGS) $(INCLUDES) $(LIBS)


clean:
	rm -f $(BIN) *.o *.cu_o *.cubin



