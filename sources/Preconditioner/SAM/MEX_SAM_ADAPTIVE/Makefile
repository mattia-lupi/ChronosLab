# Variables
CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17 -ffast-math -flto
TARGET = main
SRCS = main.cpp cpt_sam_adaptive_left.cpp qr_functions.cpp find_stuff.cpp cpt_resRho.cpp
OBJS = $(SRCS:.cpp=.o)

# Detect Operating System
UNAME_S := $(shell uname -s)

# --- OS Conditional Logic ---
ifeq ($(UNAME_S), Darwin)
    # macOS Configuration
    OPENBLAS_INC = -I/opt/homebrew/Cellar/openblas/0.3.33/include/
    LIBOMP_INC   = -I/opt/homebrew/Cellar/libomp/22.1.7/include/
    LIBOMP_LIB   = -L/opt/homebrew/Cellar/libomp/22.1.7/lib/

    INCLUDES = $(OPENBLAS_INC) $(LIBOMP_INC)
    # Append the library path and the explicit OpenMP library flag
    LIBS = -llapack -lblas $(LIBOMP_LIB) -lomp

	 M1_FLAGS = -mcpu=apple-m1
    
    CXXFLAGS += -Xpreprocessor -fopenmp -fPIC -mmacosx-version-min=15.0 $(M1_FLAGS)
    LDFLAGS  += -mmacosx-version-min=15.0 -flto
endif

ifeq ($(UNAME_S), Linux)
    # Linux Configuration
    INCLUDES := -I/usr/include/
    LIBS     := -llapack -lblas

    CXXFLAGS += -fopenmp -march=native
endif

# --- Build Rules ---

# Default Target: Fully Optimized Build
all: CXXFLAGS += -O3
all: $(TARGET)

# Debug Target: Debugging symbols and no optimization
debug: CXXFLAGS += -g -O0
debug: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS) -o $(TARGET) $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all debug clean
