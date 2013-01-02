# Makefile for FITS version of Debbie Bard's Aperture Mass Code
BIN_GALS           := mAp_gals.x
BIN_GRID           := mAp_grid.x

### You will need to edit this line to point to your own SDK installation
#CUDA_SDK = /path/to/my/cuda/sdk/installation
CUDA_SDK  := ${HOME}/CUDA/SDK_4.0/

### You will need to edit this line to point to your own cuda installation
CUDA_INSTALL_PATH := /opt/cuda


INCLUDES          += -I. -I$(CUDA_INSTALL_PATH)/include  -I$(CUDA_SDK)/C/common/inc 
LIBS		  := -L$(CUDA_INSTALL_PATH)/lib64  -L$(CUDA_SDK)/C/common/lib/ -L$(CUDA_SDK)/C/lib/ 
COMMONFLAGS       := -m64
CXXFLAGS          := $(COMMONFLAGS)
LDFLAGS           := -lrt -lm -lcudart

# compilers
NVCC              := $(CUDA_INSTALL_PATH)/bin/nvcc -arch sm_11 $(COMMONFLAGS)


CU_SOURCES_GALS   := mAp_gpu_gals.cu
CU_SOURCES_GRID   := mAp_gpu_grid.cu
CU_OBJS_GALS      := $(patsubst %.cu, %.cu_o, $(CU_SOURCES_GALS))
CU_OBJS_GRID      := $(patsubst %.cu, %.cu_o, $(CU_SOURCES_GRID))


%.cu_o : %.cu 
	$(NVCC) -c $(INCLUDES) -o $@ $^
all: ${BIN_GALS} ${BIN_GRID}

$(BIN_GALS): $(CPP_OBJS) $(CU_OBJS_GALS)
	$(CXX) -fPIC -o $(BIN_GALS) $(CU_OBJS_GALS)  $(LDFLAGS) $(INCLUDES) $(LIBS)

$(BIN_GRID): $(CPP_OBJS) $(CU_OBJS_GRID)
	$(CXX) -fPIC -o $(BIN_GRID) $(CU_OBJS_GRID) $(LDFLAGS) $(INCLUDES) $(LIBS)


clean:
	rm -f $(BIN) *.o *.cu_o *.cubin



