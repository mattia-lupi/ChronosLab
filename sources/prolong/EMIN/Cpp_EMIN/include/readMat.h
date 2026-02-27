//int readBinaryMatrix(char *filename, int *nr, int *nc, int *nt, int **iat, int **ja, double **coef);
//int readASCIIMatrix(char *filename, int *nr, int *nc, int *nt, int **iat, int **ja, double **coef);
int readCSRmat(int *nrows, int *ncols, int *nterms, int **ia, int **ja, double **coef, const char *fname, bool BINREAD);
