

void computeFirstQR(double *Ahat, int sizeI, int sizeJ, double *R, double *Rtriang, double *tau, double *work, int lwork, int &info);
void applyFirstQt(double *Ahat, int sizeI, int sizeJ, double *tau, double *a0k, double *work, int lwork, int &info);
void applyR(int sizeJ, double *R, double *a0k, int &info);
void applyOldQ(int t, int oldSizeJ, int sizeJ, int oldSizeI, int sizeI, 
               double *Ahat, double *tau, double *a0k, double *work, int lwork, int &info);
void computeNewQR(int rowSizeB2, int colSizeB2, int startB2, int sizeB, int oldSizeTau, double *Ahat,
                  double *tau, double *R, int sqrtStartR, double *Rtriang, double *work, int lwork, int &info);
