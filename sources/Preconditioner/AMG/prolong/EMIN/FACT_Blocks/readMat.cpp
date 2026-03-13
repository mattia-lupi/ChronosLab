#include <stdio.h>
#include <stdlib.h>
#include <libgen.h>
/////////////////////
//#include <iostream>
//using namespace std;
/////////////////////

/* Header structure */
struct header { int nr,nt; };

/* Row structure */
struct row { int ii,jj; double aa; };

void irow2iat(int nr,int nt,int *irow,int *iat)
{
  int j,k,irow_old,irow_new;
  irow_old = 0;
  for (k=0;k<nt;k++) {
    irow_new = irow[k];
    if (irow_new > irow_old) {
      for (j=irow_old;j<irow_new;j++)
        iat[j] = k;
      irow_old = irow_new;
    }
  }
  k = nt;
  for (j=irow_old;j<nr+1;j++)
    iat[j] = k;
}

int readBinaryMatrix(const char *filename, int *nr, int *nc, int *nt,
                     int **iat, int **ja, double **coef)
{
  int ierr = 0;
  FILE *fid;
  fid = fopen(filename, "rb");
  if (!fid) {
    ierr = -1;
    return(ierr);
  }

  struct header hdr_var;
  struct row row_var;

  if (!fread(&hdr_var, sizeof(struct header), 1, fid)) {
    fclose(fid);
    ierr = -2;
    return(ierr);
  }

  *nr = hdr_var.nr;
  *nc = *nr;
  *nt = hdr_var.nt;
  printf("Matrix rows %d columns %d nterm %d\n",*nr,*nc,*nt);

  *iat = (int*) malloc((*nr+1) * sizeof(int));
  int *irow = (int*) malloc(*nt * sizeof(int));
  *ja = (int*) malloc(*nt * sizeof(int));
  *coef = (double*) malloc(*nt * sizeof(double));

  if (irow==NULL) printf("ERROR malloc irow\n");
  if (*iat==NULL) printf("ERROR malloc iat\n");
  if (*ja==NULL) printf("ERROR malloc ja\n");
  if (*coef==NULL) printf("ERROR malloc coef\n");

  int i;
  for (i=0;i<*nt;i++) {
    if (!fread(&row_var, sizeof(struct row), 1, fid)) {
      ierr = i+1;
      return(ierr);
    }
    irow[i] = row_var.ii;
    (*ja)[i] = row_var.jj;
    (*coef)[i] = row_var.aa;
    (*ja)[i]--;
  }
  fclose(fid);

  irow2iat(*nr,*nt,irow,*iat);

  free(irow);

  return(ierr);
}

int readASCIIMatrix(const char *filename, int *nr, int *nc, int *nt,
                    int **iat, int **ja, double **coef)
{
  int ierr = 0;
  FILE *fid;
  fid = fopen(filename,"r");
  if (!fid) {
    ierr = -1;
    return(ierr);
  }
  if (!fscanf(fid, "%d %d", nr, nt)) {
    fclose(fid);
    ierr = -2;
    return(ierr);
  }
  *nc = *nr;
  printf("Matrix rows %d columns %d nterm %d\n",*nr,*nc,*nt);

  *iat = (int*) malloc((*nr+1) * sizeof(int));
  int *irow = (int*) malloc(*nt * sizeof(int));
  *ja = (int*) malloc(*nt * sizeof(int));
  *coef = (double*) malloc(*nt * sizeof(double));

  if (irow==NULL) printf("ERROR malloc irow\n");
  if (*iat==NULL) printf("ERROR malloc iat\n");
  if (*ja==NULL) printf("ERROR malloc ja\n");
  if (*coef==NULL) printf("ERROR malloc coef\n");

  int i;
  for (i=0;i<*nt;i++) {
    if (fscanf(fid, "%d %d %lf", &irow[i], &(*ja)[i], &(*coef)[i]) != 3) {
      fclose(fid);
      ierr = i+1;
      return(ierr);
    }
    (*ja)[i]--;
  }
  fclose(fid);

  irow2iat(*nr,*nt,irow,*iat);

  free(irow);

  return(ierr);
}

int readCSRmat(int *nrows, int *ncols, int *nterm, int **ia, int **ja, double **coef,
               const char *fname, bool BINREAD) {

   printf("Matrix file %s\n",fname);

   int ierr;
   if (BINREAD) {
       ierr = readBinaryMatrix(fname,nrows,ncols,nterm,ia,ja,coef);
       if (ierr != 0) return 1;
   } else {
       ierr = readASCIIMatrix(fname,nrows,ncols,nterm,ia,ja,coef);
       if (ierr != 0) return 1;
   }

   return 0;
}
