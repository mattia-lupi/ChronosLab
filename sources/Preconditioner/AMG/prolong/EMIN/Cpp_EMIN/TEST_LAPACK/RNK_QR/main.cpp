#include <iostream>  // to use: cout,endl
#include <iomanip>
#include <stdlib.h>  // to use: exit
#include <fstream>   // to use: ifstream,ofstream
#include <sstream>   // to use: stringstream
#include <math.h>
#include <libgen.h>  // to use: basename
#include <chrono>    // to use: system_clock,duration
#include "cblas.h"    
#include "lapacke.h"  // to use: dpotrf,dpotrs
#include "lapacke_utils.h"

using namespace std;

#define charStrLen 2000
#define CONDMAX 1.e+5

// MAIN PROGRAM
int main(int argc, const char* argv[]){

// Check arguments
if (argc < 2) {
   printf("Too few arguments.\n Usage: [parm file]\n");
   exit(1);
}

// Input matrix
int nn, mm, kk;

// Open the input file
ifstream fileBB(argv[1]);
// Read header
int k = 0;
fileBB >> nn >> mm;
kk = min(mm,nn);
// Allocate and read BB matrix (BB is nn x mm but it is stored by columns)
double *g       = (double*) malloc( mm*sizeof(double) );
double *zn      = (double*) malloc( nn*sizeof(double) );
double *res     = (double*) malloc( mm*sizeof(double) );
double *p0      = (double*) malloc( nn*sizeof(double) );
double *delta   = (double*) malloc( nn*sizeof(double) );
double **BB     = (double**) malloc( mm*sizeof(double*) );
double *BBbuf   = (double*) malloc( (mm*nn)*sizeof(double) );
double *TTbuf   = (double*) malloc( (mm*nn)*sizeof(double) );
double *mat_R   = (double*) malloc( (mm*nn)*sizeof(double) );
double *mat_RRT = (double*) malloc( (mm*mm)*sizeof(double) );
for (int i = 0; i < mm; i++){
   BB[i] = &(BBbuf[k]);
   k += nn;
}
for (int i = 0; i < nn; i++){
   for (int j = 0; j < mm; j++) fileBB >> BB[j][i];
}
for (int i = 0; i < mm; i++) fileBB >> g[i];
for (int i = 0; i < nn; i++) fileBB >> p0[i];
for (int i = 0; i < nn; i++) fileBB >> zn[i];

cout << "Tot memory: " << ( 2*mm+3*nn+mm*(3*nn+mm) )*sizeof(double) << " bytes " << endl;
cout << sizeof(double) << endl;

// Close the input fileBB
fileBB.close();
cout << "MATRIX read " << nn << " " << mm  << endl;

/*
FILE *of = fopen("XXX","w");
for (int i = 0; i < nn; i++){
   for (int j = 0; j < mm; j++) fprintf(of," %15.6e",BB[j][i]);
   fprintf(of,"\n");
}
fclose(of);
*/
cout << g[0] << endl;

// BB is nn x mm
// TT is mm x nn, so its factorization is Q mm x mm and R mm x nn


double query_work_1;
double query_work_2;
double *work;
double *tau;
lapack_int ierr;
lapack_int l_nn = static_cast<lapack_int>(nn);
lapack_int l_mm = static_cast<lapack_int>(mm);
lapack_int l_kk = static_cast<lapack_int>(kk);
lapack_int *JPVT;

// Query workspace
lapack_int lwork = -1;
cout << l_mm << endl;

// Workspace for QR of transpose(BB)
ierr = LAPACKE_dgeqp3_work(LAPACK_COL_MAJOR,l_mm,l_nn,BBbuf,l_mm,JPVT,tau,&query_work_1,lwork);
ierr = LAPACKE_dormqr_work(LAPACK_COL_MAJOR,'L','T',l_mm,1,l_kk,BBbuf,l_mm,tau,res,l_mm,
                           &query_work_2,lwork);
cout << "query_work_1 " << query_work_1 << endl;
cout << "query_work_2 " << query_work_2 << endl;
lwork = static_cast<lapack_int>(max(query_work_1,query_work_2));
cout << "lwork " << lwork << endl;

// Transpose BB
LAPACKE_dge_trans(LAPACK_COL_MAJOR,l_nn,l_mm,BBbuf,l_nn,TTbuf,l_mm);

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
FILE *of = fopen("XXX","w");
for (int i = 0; i < nn; i++){
   for (int j = 0; j < mm; j++) fprintf(of," %15.6e",BB[j][i]);
   fprintf(of,"\n");
}
fprintf(of,"\n");
fprintf(of,"\n");
for (int i = 0; i < mm; i++){
   for (int j = 0; j < nn; j++) fprintf(of," %15.6e",TTbuf[j*mm+i]);
   fprintf(of,"\n");
}
fprintf(of,"\n");
fprintf(of,"\n");
for (int j = 0; j < nn*mm; j++) fprintf(of," %15.6e\n",TTbuf[j]);
fclose(of);
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// Compute residual: res = g - BT*p0
for (int i = 0; i < mm; i++) res[i] = g[i];
cblas_dgemv(CblasColMajor,CblasNoTrans,mm,nn,-1.0,TTbuf,mm,p0,1,1.0,res,1);

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
cout << "RES: " << endl;
for (int i = 0; i < mm; i++) cout << res[i] << endl;
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// Allocate workspace
JPVT = (lapack_int*) malloc( nn*sizeof(lapack_int) );
tau = (double*) malloc( min(nn,mm)*sizeof(double) );
work = (double*) malloc( lwork*sizeof(double) );

// Perform QR with column pivoting
for (int i = 0; i < nn; i++) JPVT[i] = 0;
ierr = LAPACKE_dgeqp3_work(LAPACK_COL_MAJOR,l_mm,l_nn,TTbuf,l_mm,JPVT,tau,work,lwork);
cout << "IERR DGEQP3: " << ierr << endl;

// Find rank of B
int rank = 0;
while (fabs( TTbuf[0] / TTbuf[rank*mm+rank] ) < CONDMAX && rank < kk) rank++;

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2
cout << "RANK: " << rank << endl;
double fac = TTbuf[0];
cout << fac << endl;
for (int i = 0; i < l_kk; i++) cout << fabs(TTbuf[i*mm+i] / fac) << endl;
for (int i = 0; i < l_nn; i++) cout << JPVT[i] << endl;
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2

// Multiply res by QT: res = QT*res
lapack_int l_rank = static_cast<lapack_int>(rank);
ierr = LAPACKE_dormqr_work(LAPACK_COL_MAJOR,'L','T',l_mm,1,l_mm,TTbuf,l_mm,tau,res,l_mm,
                           work,lwork);
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2
cout << "QT RES: " << endl;
for (int i = 0; i < l_rank; i++) cout << res[i] << endl;
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2

// Compute R*RT
for (int j = 0; j < nn; j++){
   int col = l_mm*j;
   for (int i = 0; i < rank; i++) mat_R[col+i] = (j>=i) ? TTbuf[col+i]:0.0;
}
cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,rank,rank,nn,1.0,mat_R,l_mm,mat_R,l_mm,
            0.0,mat_RRT,l_mm);
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2
cout << "mat_RRT: " << endl;
for (int i = 0; i < rank; i++){
   for (int j = 0; j < rank; j++) cout << "  " << mat_RRT[j*l_mm+i];
   cout << endl;
}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2

// Compute res = inv(mat_RRT)*res;
ierr = LAPACKE_dpotrf_work(LAPACK_COL_MAJOR,'U',l_rank,mat_RRT,l_mm);
ierr = LAPACKE_dpotrs_work(LAPACK_COL_MAJOR,'U',l_rank,1,mat_RRT,l_mm,res,l_mm);
//ierr = LAPACKE_dposv_work(LAPACK_COL_MAJOR,'U',l_rank,1,mat_RRT,l_mm,res,l_mm);
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2
cout << "SOLV RES: " << endl;
for (int i = 0; i < l_rank; i++) cout << res[i] << endl;
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2
//FINO QUIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

// Compute delta = mat_R^T * res
cblas_dgemv(CblasColMajor,CblasTrans,rank,nn,1.0,mat_R,mm,res,1,0.0,delta,1);
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2
cout << "DELTA: " << endl;
for (int i = 0; i < nn; i++) cout << delta[i] << endl;
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2

// Update p0
for (int i = 0; i < nn; i++) p0[JPVT[i]-1] += delta[i];
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2
cout << "p2: " << endl;
for (int i = 0; i < nn; i++) cout <<p0[i] << endl;
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2

// mat_RRT = U'*U ---> inv(U)*inv(U') ---> 

// Create the orthogonal projector
cblas_dtrsm(CblasColMajor,CblasLeft,CblasUpper,CblasTrans,CblasNonUnit,
            rank,nn,1.0,mat_RRT,mm,mat_R,mm);
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2
cout << "mat_Q: " << endl;
for (int i = 0; i < rank; i++){
   for (int j = 0; j < nn; j++) cout << "  " << mat_R[j*l_mm+i];
   cout << endl;
}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2

// Permute columns of Q with PVT
for (int j = 0; j < nn; j++){
   int col = j*l_mm;
   int pcol = (JPVT[j]-1)*l_mm;
   for (int i = 0; i < rank; i++) TTbuf[pcol+i] = mat_R[col+i];
}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2
cout << "mat_QT perm: " << endl;
for (int i = 0; i < rank; i++){
   for (int j = 0; j < nn; j++) cout << "  " << TTbuf[j*l_mm+i];
   cout << endl;
}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2

// Test with zn ==> res = BT*zn
cblas_dgemv(CblasColMajor,CblasTrans,nn,mm,1.0,BBbuf,nn,zn,1,0.0,res,1);
//for (int i = 0; i < nn; i++) cout << zn[i] << endl;
cout << "NORM res: " << cblas_dnrm2(mm,res,1) << endl;

// Project zn
cblas_dgemv(CblasColMajor,CblasNoTrans,rank,nn,1.0,TTbuf,mm,zn,1,0.0,res,1);
cblas_dgemv(CblasColMajor,CblasTrans,rank,nn,-1.0,TTbuf,mm,res,1,1.0,zn,1);
//for (int i = 0; i < nn; i++) cout << zn[i] << endl;
cblas_dgemv(CblasColMajor,CblasTrans,nn,mm,1.0,BBbuf,nn,zn,1,0.0,res,1);
cout << "NORM res: " << cblas_dnrm2(mm,res,1) << endl;


exit(0);

}
