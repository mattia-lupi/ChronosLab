//----------------------------------------------------------------------------------------

#pragma once

#include <vector>     // to use: vector
#include <cmath>      // to use: abs,sqrt
#include <algorithm>  // to use: min
#include <lapacke.h>  // to use: dpotrf,dpotrs
#include <omp.h>      // omp_get_wtime
#include <iostream>
#include <iomanip>

using namespace std;

//----------------------------------------------------------------------------------------

// Regular integer
typedef int iReg;

// Extended integer
typedef long int iExt;

// Regular real
typedef float rReg;

// Extended real
typedef double rExt;

// Blas integer
typedef int iBlas;

// Lapack integer
typedef int iLapack;

//----------------------------------------------------------------------------------------


rExt cpt_ddot(iReg n, rExt *x1, rExt *x2);

rExt cpt_dnrm2(iReg n, rExt *x1);

void copy_ja(iExt n, iReg *x1, iReg *x2);

void copy_coef(iExt n, rExt *x1, rExt *x2);

void ri_sortsplit(iReg n, iReg ncut, rExt * const R_vec, iReg * const I_vec);

void swapi(iReg &i1, iReg &i2);

void iheapsort(iReg * const x1, iReg n);

void gather_fullsys(bool &nulrhs, iReg irow, iReg mrow, iReg jendbloc, iReg nequ,
                    iExt nterm, iReg mmax, const iExt * const iat, const iReg * const ja,
                    const iReg * const vecinc, const rExt * const mat_A,
                    rExt * const full_A, rExt * const rhs);

void kap_grad(iReg irow, iReg nequ, iExt nterm, iReg mmax, iReg &mrow, iReg jendbloc,
              iReg lfil, const iExt * const iat, const iReg * const ja,
              const rExt * const mat_A, const rExt * const rhs, iReg * const IWN,
              iReg * const JWN, rExt * const WR);

void cpt_afsai_coef(iReg chunk_size, iReg n_step, iReg step_size, rExt tau, rExt eps,
                    iReg shift, iReg nrows, iReg nequ, iExt nterm, iExt &nterm_G,
                    const iExt * const iat, const iReg * const ja, iExt * const istart_G,
                    iExt * const istop_G, iReg * const ja_G, const rExt * const coef_A,
                    rExt * const coef_G);


