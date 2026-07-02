

void computeFirstQR(double *Ahat, int sizeI, int sizeJ, double *R, double *Rtriang, double *tau, double *work, int lwork, int &info);
void applyFirstQt(double *Ahat, int sizeI, int sizeJ, double *tau, double *a0k, double *work, int lwork, int &info);
void applyR(int sizeJ, double *R, double *a0k, int &info);
void applyQt(int t, int *sizeJ, int *sizeI, int *qStart, double *Ahat, 
             double *tau, double *a0k, int nrowsA0k, int ncolsA0k, double *work, int lwork, int &info);
void computeNewQR(int t, int *sizeI, int *sizeJ, int *qStart, double *Ahat, double *tau, double *R, 
                  double *Rtriang, double *work, int lwork, int &info);
