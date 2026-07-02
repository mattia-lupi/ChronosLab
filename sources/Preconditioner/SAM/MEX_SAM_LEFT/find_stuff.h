void fullA0k(int nn_A, int *iat0, int *ja0, double *coef0, int k, double *A0k);
void findNonZeroInColJ(int *J, int *iatk, int *jak, int n2, int *I, int &sizeI);
void getA0k(double *a0k,  int *I, int sizeI, int oldSizeI, int *iat0,  int *ja0, double *coef0,  int k);
void getAhat( int *I,  int sizeI,  int *J,  int Jstart,  int Jend,
             int *iatk,  int *jak, double *coefk, double *Ahat,  int &Astart);
void getAJ(int *J, int Jsize, int nn_A, int *iatk, int *jak, double *coefk, double *AJ);

void fillL(int *L, double *res, int nn_A, int &usedL);
void findJtilde(int *Jtilde, int &JtildeSize, int *L, int sizeL, int *iatk, int *jak, int *J, int sizeJ);
void fullAJtilde(int nn_A, int *iatk, int *jak, double *coefk, int *Jtilde, int JtildeSize, double *AJtilde);
