CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17
TARGET = main
SRCS = main.cpp cpt_sam_adaptive_left.cpp qr_functions.cpp find_stuff.cpp cpt_resRho.cpp transpose.cpp
OBJS = $(SRCS:.cpp=.o)

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S), Darwin)
    OPENBLAS_INC = -I/opt/homebrew/opt/openblas/include
    OPENBLAS_LIB = -L/opt/homebrew/opt/openblas/lib
    LIBOMP_INC   = -I/opt/homebrew/opt/libomp/include
    LIBOMP_LIB   = -L/opt/homebrew/opt/libomp/lib
    
    INCLUDES = $(OPENBLAS_INC) $(LIBOMP_INC)
    LIBS = $(OPENBLAS_LIB) $(LIBOMP_LIB) -llapack -lblas -lomp

    M1_FLAGS = -mcpu=apple-m1
    CXXFLAGS += -Xpreprocessor -fopenmp -fPIC $(M1_FLAGS)
endif

ifeq ($(UNAME_S), Linux)
    INCLUDES := -I/usr/include/
    LIBS     := -llapack -lblas
    CXXFLAGS += -fopenmp -march=native
endif

all: CXXFLAGS += -O3 -ffast-math -flto
all: LDFLAGS += -flto
all: $(TARGET)

xtrace: CXXFLAGS += -g3 -O0 -fno-omit-frame-pointer -fno-inline
ifeq ($(UNAME_S), Darwin)
xtrace: GEN_DSYM = dsymutil $(TARGET)
endif
xtrace: $(TARGET)

debug: CXXFLAGS += -g -O0
debug: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS) -o $(TARGET) $(LIBS)
	$(GEN_DSYM)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -rf $(OBJS) $(TARGET) $(TARGET).dSYM

.PHONY: all xtrace debug clean
