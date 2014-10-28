# Makefile
#
# author      : Matt Bellis
# date        : Jan., 2012
#

################################################################################
# This section contains some environment variables that you may want to edit, or
# perhaps need to edit, depending on your CUDA installation.
################################################################################

# INSTALL_DIR is where the CUDA and C executables will be installed. You may find it 
# useful to install in some system-wide directory (e.g. /usr/local/bin, 
# $(HOME)/bin, etc.).
INSTALL_DIR       := ../bin/

# SDK_INSTALL_PATH and CUDA_INSTALL_PATH are probably the more significant
# environment variables and will depend on your system and where the CUDA
# libraries and header files are installed.
#SDK_INSTALL_PATH  := ${HOME}/CUDA/SDK_4.0/
#CUDA_INSTALL_PATH := /opt/cuda/
SDK_INSTALL_PATH  := $(HOME)/NVIDIA_GPU_Computing_SDK
CUDA_INSTALL_PATH := /usr/local/cuda/

# While the INCLUDES and LIBS flags are dependent on the install path, there may 
# still be some subdirectories that are machine-dependent. 
INCLUDES          += -I. -I$(SDK_INSTALL_PATH)/C/common/inc/ -I$(CUDA_INSTALL_PATH)/include/
LIBS              += -L. -L$(SDK_INSTALL_PATH)/C/lib/ -L$(CUDA_INSTALL_PATH)//usr/lib/x86_64-linux-gnu/
LIBS              += -L$(CUDA_INSTALL_PATH)/lib64

LDFLAGS := -lcudart

COMMONFLAGS       := -m64 
OBJ_DIR           := ./objects/
CXXFLAGS          := $(COMMONFLAGS)

# compilers
CC              := g++
NVCC              := $(CUDA_INSTALL_PATH)/bin/nvcc -arch sm_20 $(COMMONFLAGS)

################################################################################

################################################################################
# Rules for Makefile
################################################################################
all: 3pcf 3pcf_C_version

3pcf : $(OBJ_DIR)/3pcf.cu_o
	mkdir -p $(INSTALL_DIR)
	$(CXX) -fPIC -o $(INSTALL_DIR)/3pcf $(OBJ_DIR)//3pcf.cu_o $(LDFLAGS) $(INCLUDES) $(LIBS)

$(OBJ_DIR)/3pcf.cu_o : 3pcf.cu
	mkdir -p $(OBJ_DIR)
	$(NVCC) -c $(INCLUDES) -o $(OBJ_DIR)/3pcf.cu_o $<

3pcf_C_version : 3pcf_C_version.cc
	mkdir -p $(INSTALL_DIR)
	$(CC) -O1 -o $(INSTALL_DIR)/$@ $@.cc -pg
	#$(CC) -pg -O1 -o $(INSTALL_DIR)/$@ $@.cc

clean:
	rm -f $(INSTALL_DIR)/3pcf_C_version 
	rm -f $(INSTALL_DIR)/3pcf 
	rm -f $(OBJ_DIR)/*.o 
	rm -f $(OBJ_DIR)/*.cu_o 

################################################################################
