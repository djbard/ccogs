# Makefile
#
# author      : Matt Bellis
# date        : Jan., 2012
# discription : generic Makefile for making cuda programs
#
INSTALL_DIR       := ../bin/
OBJ_DIR           := ./objects/

#SDK_INSTALL_PATH  := /home/bellis/NVIDIA_GPU_Computing_SDK
#CUDA_INSTALL_PATH := /usr/
SDK_INSTALL_PATH  := ${HOME}/CUDA/SDK_4.0/
CUDA_INSTALL_PATH := /opt/cuda/

# flags

INCLUDES          += -I. -I$(SDK_INSTALL_PATH)/C/common/inc/ -I$(CUDA_INSTALL_PATH)/include/
LIBS              += -L. -L$(SDK_INSTALL_PATH)/C/lib/ -L$(CUDA_INSTALL_PATH)//usr/lib/x86_64-linux-gnu/
LIBS              += -L$(CUDA_INSTALL_PATH)/lib64

CXXFLAGS          := $(COMMONFLAGS)

#LDFLAGS := -lSystem.B -lm -lcudart -lcutil_x86_64 
#LDFLAGS := -lm -lcudart 
LDFLAGS := -lcudart

# compilers
CC              := g++
#NVCC              := nvcc -arch sm_11 $(COMMONFLAGS) --device-emulation
NVCC              := $(CUDA_INSTALL_PATH)/bin/nvcc -arch sm_20 $(COMMONFLAGS)

all: angular_correlation angular_correlation_C_version

angular_correlation : $(OBJ_DIR)/angular_correlation.cu_o
	mkdir -p $(INSTALL_DIR)
	$(CXX) -fPIC -o $(INSTALL_DIR)/angular_correlation $(OBJ_DIR)//angular_correlation.cu_o $(LDFLAGS) $(INCLUDES) $(LIBS)

$(OBJ_DIR)/angular_correlation.cu_o : angular_correlation.cu
	mkdir -p $(OBJ_DIR)
	$(NVCC) -c $(INCLUDES) -o $(OBJ_DIR)/angular_correlation.cu_o $<

angular_correlation_C_version : angular_correlation_C_version.cc
	mkdir -p $(INSTALL_DIR)
	$(CC) -o $(INSTALL_DIR)/$@ $@.cc

clean:
	rm -f $(INSTALL_DIR)/angular_correlation_C_version 
	rm -f $(INSTALL_DIR)/angular_correlation 
	rm -f $(OBJ_DIR)/*.o 
	rm -f $(OBJ_DIR)/*.cu_o 
