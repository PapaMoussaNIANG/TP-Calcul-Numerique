/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include <string.h>
void eig_poisson1D(double* eigval, int *la){
  // TODO: Compute all eigenvalues for the 1D Poisson operator
  for (size_t i = 1; i < *la; ++i)
  {
    eigval[i] = 2. - 2. * cos(((double)i * M_PI) / ((double)(*la) + 1.));
  }
}

double eigmax_poisson1D(int *la){
  // TODO: Compute and return the maximum eigenvalue for the 1D Poisson operator
  return 2. - 2. * cos(((double)(*la) * M_PI) / ((double)(*la) + 1.));
}

double eigmin_poisson1D(int *la){
  // TODO: Compute and return the minimum eigenvalue for the 1D Poisson operator
  return 2. - 2. * cos(M_PI / ((double)(*la) + 1.));
}

double richardson_alpha_opt(int *la){
  return 2. / (eigmax_poisson1D(la) + eigmin_poisson1D(la));
}

/**
 * Solve linear system Ax=b using Richardson iteration with fixed relaxation parameter alpha.
 * The iteration is: x^(k+1) = x^(k) + alpha*(b - A*x^(k))
 * Stops when ||b - A*x^(k)||_2  / ||b||_2 < tol or when reaching maxit iterations.
 */
void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  // TODO: Implement Richardson iteration
  // 1. Compute residual r = b - A*x (use dgbmv for matrix-vector product)
  // 2. Update x = x + alpha*r (use daxpy)
  // 3. Check convergence: ||r||_2 < tol (use dnrm2)
  // 4. Store residual norm in resvec and repeat
  double* r = NULL;
  const double TOLERANCE = *tol;
  const int MAX_IT = *maxit;
  const double RHS_norm = cblas_dnrm2(*la, RHS, 1);
  r = (double*)malloc(*la * sizeof(double));
  memcpy(r, RHS, *la  * sizeof(double));
  cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1, AB, *lab, X, 1, 1.0, r, 1);
  resvec[*nbite] = cblas_dnrm2(*la, r, 1) / RHS_norm;
  ++(*nbite);
  while(resvec[(*nbite) - 1] > TOLERANCE)
  {
    cblas_daxpy(*la, *alpha_rich, r, 1, X, 1);
    memcpy(r, RHS, *la * sizeof(double));
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1, AB, *lab, X, 1, 1.0, r, 1);
    resvec[*nbite] = cblas_dnrm2(*la, r, 1) / RHS_norm;
    ++(*nbite);
    if(*nbite ==  MAX_IT)
      break;
  }
}

/**
 * Extract MB for Jacobi method from tridiagonal matrix.
 * Such as the Jacobi iterative process is: x^(k+1) = x^(k) + D^(-1)*(b - A*x^(k))
 */
void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  // TODO: Extract diagonal elements from AB and store in MB
  // MB should contain only the diagonal of A
}

/**
 * Extract MB for Gauss-Seidel method from tridiagonal matrix.
 * Such as the Gauss-Seidel iterative process is: x^(k+1) = x^(k) + (D-E)^(-1)*(b - A*x^(k))
 */
void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  // TODO: Extract diagonal and lower diagonal from AB
  // MB should contain the lower triangular part (including diagonal) of A
}

/**
 * Solve linear system Ax=b using preconditioned Richardson iteration.
 * The iteration is: x^(k+1) = x^(k) + M^(-1)*(b - A*x^(k))
 * where M is either D for Jacobi or (D-E) for Gauss-Seidel.
 * Stops when ||b - A*x^(k)||_2  / ||b||_2 < tol or when reaching maxit iterations.
 */
void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  // TODO: Implement Richardson iterative method
}

