mex -R2018a SYM_aFSAI_compute.cpp compute_local_fsai.cpp cpt_afsai_coef.cpp CXXFLAGS="-fopenmp -fPIC" LDFLAGS="-fopenmp -z noexecstack" -lmwblas -lmwlapack
