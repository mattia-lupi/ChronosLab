#include "cpt_resRho.h"
#include <vector>
#include <iostream>
#include "lapack.h"
#include "cblas.h"

void cptRes(int nn_A, int sizeJ, double *A0k, double *AJ, double *mHat, double *res, double &resRelNorm, double &resNorm){
   // Compute the product AJ*mHat and save it in AJ
   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nn_A, 1, sizeJ, 1.0, AJ, sizeJ, mHat, 1, 0.0, res, 1);

   // print_vector("Matrix Aj", nn_A, res);
   // Compute the norm of Aj*mHat
   double normAjMh = cblas_dnrm2(nn_A,res,1);

   // Compute res = Aj*mHat - A0k and save it in AJ
   cblas_daxpy(nn_A, -1.0, A0k, 1, res, 1);

   // Compute the residual norm
   resNorm = cblas_dnrm2(nn_A,res,1);

   // Compute the relative version of the norm
   resRelNorm = 2*resNorm/(normAjMh+resRelNorm);
   // printf("resNorm %f\n", resRelNorm);

   return;
}


// Check when JtildeSize >= 2
void cptRhoJ2(int JtildeSize, double *normColJ, double *AJtilde, int nn_A, double *res, double normRes){
   // Compute the norm for each column. 
   // The matrix is saved in column major so doing the norm is easy
   for (int i = 0; i < JtildeSize; ++i){
      normColJ[i] = cblas_dnrm2(nn_A,&(AJtilde[i*nn_A]),1);
      normColJ[i] *= normColJ[i];

      // If the norm is zero then discard this column
      if (normColJ[i] < 1e-15){
         normColJ[i] = 1e15;
      }
   }

   // print_vector("Vector res", JtildeSize, res);

   // Compute the product AJtilde^T*res and save it in AJtilde
   cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, JtildeSize, 1, nn_A, 1.0, AJtilde, nn_A, res, nn_A, 0.0, AJtilde, nn_A);

   // print_vector("Vector rat", JtildeSize, AJtilde);
   double temp;
   // Compute the rhoJ2 and save it in normColJ
   for (int i = 0; i < JtildeSize; ++i){
      temp = AJtilde[i];
      normColJ[i] = normRes*normRes - temp*temp/normColJ[i];
      // Avoid having possibly negative rhos 
      normColJ[i] = std::max(normColJ[i],0.);
   }

   return;
}

// Compute the index in which rhoJ2 is minimum
int minIdx(double *rhoJ2, int JtildeSize){
   int idx = 0;
   // Initialize the minimum
   double minimum = rhoJ2[0];

   // Loop over the rhoJ2 to get the minimum
   for (int i = 1; i < JtildeSize; ++i){
      if (rhoJ2[i] < minimum){
         idx = i;
         minimum = rhoJ2[i];
      }
   }

   return idx;
}