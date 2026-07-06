# Variables
CXX = g++
CXXFLAGS = -g -Wall -Wextra -std=c++17
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
    
    INCLUDES = $(OPENBLAS_INC) $(LIBOMP_INC)
    LIBS = -llapack -lblas
endif

ifeq ($(UNAME_S), Linux)
    # Linux Configuration
    INCLUDES := -I/usr/include/
    LIBS     := -llapack -lblas
    
    CXXFLAGS += -fopenmp
endif

# --- Build Rules ---

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET) $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
