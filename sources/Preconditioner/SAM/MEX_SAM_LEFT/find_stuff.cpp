#include "find_stuff.h"
#include <vector>
#include <iostream>

void fullA0k(int nn_A, int *iat0, int *ja0, double *coef0, int k, double *A0k){
   // Loop over the rows
   for (int row = 0; row < nn_A; ++row){
      // Loop over the single row entries
      for (int i = iat0[row]; i < iat0[row+1]; ++i){
         // If the column is the same of the current column k
         // Get the nonzero index
         if (ja0[i] == k){
            A0k[row] = coef0[i];
            break;
         }

         // No entry was found, set to zero
         if (i == iat0[i+1] - 1){
            A0k[row] = 0;
         }
      }
   }
   return;
}
void findNonZeroInColJ(int *J, int *iatk, int *jak, int n2, int *I, int &sizeI){

   // Assume the pattern is symmetric
   // the nonzero row entries in column J correspond to 
   // the nonzero column entries in row J

   int initial_count = sizeI;

   // Flag for fast exit in case of repeated index
   bool skip = false;
   // Cycle over J to select the rows to get, need to check for no repeated indices
   for (int i = 0; i < n2; ++i){
   	// printf("J[%d] %d\n", i,J[i]);
      for (int j = iatk[J[i]]; j < iatk[J[i] + 1]; ++j){
      	// printf("from %d to %d\n",iatk[J[i]],iatk[J[i] + 1]);
         // If the sizeI is 0 then avoid checks and just copy all the first row
         if(initial_count != 0){
            // Check for duplicate value
            for (int q = 0; q < sizeI; ++q){
            	// printf("\t\tjak %d, Iq %d, sizeI %d\n", jak[j],I[q],sizeI);
               if(jak[j] == I[q]){
                  skip = true;
                  break;
               }
            }
   
            // This index is repeated, fast exit
            if (skip){
               // Reset flag to false
               skip = false;
            }
            else{
            	// This index is new
         		// Copy this nonzero column entry inside I
         		I[sizeI] = jak[j];
         		sizeI++;
            }
         }
         else{
         	// This index is new
         	// Copy this nonzero column entry inside I
         	I[sizeI] = jak[j];
         	sizeI++;
         }
      }
   }
   return;
}


void getA0k(double *a0k,  int *I, int sizeI, int oldSizeI, int *iat0,  int *ja0, double *coef0,  int k){
   // Cycle over I to get which columns to seach for
   // Cycle over only the new entries of I
   for (int i = oldSizeI; i < sizeI; ++i){
      for (int j = iat0[I[i]]; j < iat0[I[i] + 1]; ++j){
         if (ja0[j] == k){
            // Get the entry of column k
            a0k[i] = coef0[j];
            break;
         }

         if(j == iat0[I[i] + 1] - 1){
            // If entered here there is no column entry equal to k 
            // set to zero
            a0k[i] = 0;
         }
      }
   }
   return;
}

void getAhat(int *I, int sizeI, int *J, int Jstart, int Jend,
             int *iatk, int *jak, double *coefk, double *Ahat, int &Astart){
   // Cycle over the columns in J that have been added
   for (int j = Jstart; j < Jend; ++j){
      // Cycle over the rows in I
      for (int i = 0; i < sizeI; ++i){
         // Cycle over the chosen row
         for(int k = iatk[I[i]]; k < iatk[I[i]+1]; ++k){
            // If the column in the row coincides with the column added then get the nonzero value
            // printf("%d %d %d ", k, jak[k], J[j]);
            if(jak[k] == J[j]){
               Ahat[Astart] = coefk[k];
               // printf("%f\n", coefk[k]);
               Astart++;
               break;
            }

            if(k == iatk[I[i] + 1] - 1){
               // If entered here there is no column entry equal to k 
               // set to zero
               // printf("0\n");
               Ahat[Astart] = 0;
               Astart++;
            }
            // else{
            //    // printf("\n");
            // }

         }
      }
   }
   return;
}



// AJ is saved rowwise then used as transposed in the blas gemm to have better memory access in creating it
void getAJ(int *J, int Jsize, int nn_A, int *iatk, int *jak, double *coefk, double *AJ){
   int Astart = 0;
   // int currJ = 0;

   // printf("%d\n",Jsize);

   // Cycle over the rows in I
   for (int i = 0; i < nn_A; ++i){
      // currJ = 0;

      // Loop over all possible J
      for (int currJ = 0; currJ < Jsize; ++currJ){
      	// Cycle over the chosen row
      	for(int k = iatk[i]; k < iatk[i+1]; ++k){
      		// If column index is equal to the chosen J column then add it
      		if(jak[k] == J[currJ]){
      			// printf("A[%d] = %f\n", Astart,coefk[k]);
      			AJ[Astart] = coefk[k];
         		Astart++;
         		break;
      		}

      		if(k == iatk[i+1] - 1){
				   // If entered here there is no column entry equal to this J
				   // set to zero
				   // printf("A[%d] = %f\n", Astart,0);
				   AJ[Astart] = 0;
				   Astart++;
				   break;
				}
      	}
         // printf("k %d, jak %d, Jcurr %d\n",k,jak[k],J[currJ]);
         // if(jak[k] > J[currJ]){
         //    // If entered here there is no column entry equal to k 
         //    // set to zero
         //    // printf("0\n");
         //    // printf("entered when k %d, jak %d, Jcurr %d\n",k,jak[k],J[currJ]);
         //    AJ[Astart] = 0;
         //    printf("A[%d] = %f\n", Astart,0);
         //    Astart++;
         //    currJ++;
         // }

         // // Check if reached the max size for J
         // if(currJ >= Jsize){
         //    printf("entered break at k = %d\n",k);
         //    break;
         // }

         // // If the column in the row coincides with the column needed then add
         // if(jak[k] == J[currJ]){
         //    AJ[Astart] = coefk[k];
         //    printf("A[%d] = %f\n", Astart,coefk[k]);
         //    Astart++;
         //    currJ++;

         // }

         // // Check if reached the max size for J
         // if(currJ >= Jsize){
         //    printf("entered break at k = %d\n",k);
         //    break;
         // }

         // if(k == iatk[i + 1] - 1){
         //    // If entered here there is no column entry equal to k 
         //    // set to zero
         //    printf("0\n");
         //    for (int j = currJ; j < Jsize; ++j){
         //    	AJ[Astart] = 0;
         //    	Astart++;
         //    }
         // }
      }
   }
   return;
}


void fillL(int *L, double *res, int nn_A, int &usedL){
   usedL = 0;
   // Loop over all residual entries
   for (int i = 0; i < nn_A; ++i){
      // The residual is not numerically zero
      // Add it as a possible column to be computed
      if (std::abs(res[i]) > 1e-13){
         L[usedL] = i;
         usedL++;
      }
   }
   return;
}

void findJtilde(int *Jtilde, int &JtildeSize, int *L, int sizeL, int *iatk, int *jak, int *J, int sizeJ){

   // Sanity check for sizeL == 0
   if (sizeL == 0){
      printf("sizeL == 0, check\n");
      return;
   }

   // Initialize counter to 0
   JtildeSize = 0;

   // Flag for fast exit in case of repeated index
   bool skip = false;
   int maxJsize;
   for (int i = 0; i < sizeL; ++i){
      for (int j = iatk[L[i]]; j < iatk[L[i] + 1]; ++j){
         // Get the max size to search for duplicates
         maxJsize = std::max(sizeJ,JtildeSize);
         // std::cout << "maxJsize = " << maxJsize << std::endl;
         // std::cout << "ja["<< j << "] = " << jak[j] << std::endl;
         // Check for duplicate value
         for (int q = 0; q < maxJsize; ++q){
            // Check if it is already in Jtilde
            if(q < JtildeSize){
               if(jak[j] == Jtilde[q]){
                  skip = true;
                  break;
               }
            }

            // Check if it is already in J
            if (q < sizeJ){
               if(jak[j] == J[q]){
                  skip = true;
                  break;
               }
            }
         }
   
         // This index is repeated, fast exit
         if (skip){
            // Reset flag to false
            skip = false;
         }
         else{
            // Index is not repeated, add it
            Jtilde[JtildeSize] = jak[j];
            JtildeSize++;
         }
      }
   }
   return;
}

// Compute A(:,Jtilde) in colmajor
// check again if Jtilde size >= 2
void fullAJtilde(int nn_A, int *iatk, int *jak, double *coefk, int *Jtilde, int JtildeSize, double *AJtilde){
   
   // Loop over possible columns for Jtilde
   for (int j = 0; j < JtildeSize; ++j){
      // Loop over the rows
      for (int row = 0; row < nn_A; ++row){
         // Loop over the single row entries
         for (int i = iatk[row]; i < iatk[row+1]; ++i){
         
            // If the column is the same of the current column k
            // Get the nonzero index
            if (jak[i] == Jtilde[j]){
               AJtilde[row + j*nn_A] = coefk[i];
               break;
            }
   
            // No entry was found, set to zero
            if (i == iatk[i+1] - 1){
               AJtilde[row + j*nn_A] = 0;
            }
         }
      }
   }
   return;
}
