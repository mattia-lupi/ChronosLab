mex -g -R2018a sam_adaptive_left_mex.cpp cpt_sam_adaptive_left.cpp qr_functions.cpp find_stuff.cpp cpt_resRho.cpp CXXFLAGS="-fopenmp -fPIC" LDFLAGS="-fopenmp -z noexecstack" -lmwblas -lmwlapack
