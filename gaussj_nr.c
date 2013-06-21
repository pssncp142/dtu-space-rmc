#include "math.h"
#include "common.h"

#define SWAP(a, b) {double temp = (a); (a) = (b); (b) = temp;}	/* function used in gaussj_nr */ 

int gaussj_nr(
     double a[], 
     int n, 
     double b[], 
     int m)

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                             */
/*  Niels Lund, June 2013                                  */        
/*                                                             */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/* -------------------------------------------------------------

docbegin:

*  Parent component:            gaussj_nr from NJW 
*  Programmer:                  Niels Lund (modification to avoid DAL/PIL/RIL dependencies)
*  Affiliation:                 DTU Space 
*  Purpose:                     Fit of shadowgram components for Iterative Removal of Sources
*  Origin date:                 021220
*  Update history:
*  Type of element:             Subroutine
*  Required functions           SWAP

Various remarks:

Gauss-Jordan matrix inversion from "Numerical Recipes"
adapted to follow the storage convention:                   a(i,j) = a[j*N + i]

-------------------------------------------------------------- */
{
   int icol = 0, irow = 0, i, j, k, l, ll, status;
   int indxc[20], indxr[20], ipiv[20];	/* Set up for maximally 20 equations with 20 unknowns */
   double big, dum, pivinv;


   for (j=0;j<n;j++) ipiv[j]=0;
   for (i=0;i<n;i++) {
      big=0.0;
      for (j=0;j<n;j++)
         if (ipiv[j] != 1)
            for (k=0;k<n;k++) {
               if (ipiv[k] == 0) {
                  if (fabs(a[k*n+j]) >= big) {
                     big=fabs(a[k*n+j]);
                     irow=j;
                     icol=k;
                  }
               } else {
                  if (ipiv[k] > 1) {
                     status = -1; 
                     goto exittrace;             /* singular matrix 1 */
                  }
               }
            }
      ++(ipiv[icol]);
      
      if (irow != icol) {
         for (l=0;l<n;l++) SWAP(a[l*n + irow],a[l*n + icol]);
         for (l=0;l<m;l++) SWAP(b[l*n + irow],b[l*n + icol]);
      }
      indxr[i]=irow;
      indxc[i]=icol;


      if (a[icol*n + icol] == 0.0)  {
         status = -2; 
         goto exittrace;	         /* singular matrix 2 */
      }
      
      pivinv=1.0/a[icol*n + icol];
      a[icol*n + icol]=1.0;
      for (l=0;l<n;l++) a[l*n + icol] *= pivinv;
      for (l=0;l<m;l++) b[l*n + icol] *= pivinv;
      for (ll=0;ll<n;ll++)
         if (ll != icol) {
            dum=a[icol*n + ll];
            a[icol*n + ll]=0.0;
            for (l=0;l<n;l++) a[l*n + ll] -= a[l*n + icol]*dum;
            for (l=0;l<m;l++) b[l*n + ll] -= b[l*n + icol]*dum;
         }
   }
   for (l=n-1;l>=0;l--) {
      if (indxr[l] != indxc[l])
         for (k=0;k<n;k++)
            SWAP(a[indxr[l]*n + k],a[indxc[l]*n + k]);
   }
   status = 0;
exittrace:
   return( status );

}
