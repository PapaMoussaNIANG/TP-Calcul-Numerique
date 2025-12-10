/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  // TODO: Fill AB with the tridiagonal Poisson operator
  for(int j = 0; j < *kv; ++j)
      AB[j] = 0;
    
  AB[(*kv)] = 0;
  AB[(*kv) + 1] = 2;
  AB[(*kv) + 2] = -1;

  for (int i = 1; i < *la - 1; ++i)
  {
    for(int j = 0; j < *kv; ++j)
      AB[i * (*lab) + j] = 0;
    
    AB[i * (*lab )+ (*kv)] = -1;
    AB[i * (*lab ) + (*kv) + 1] = 2;
    AB[i * (*lab ) + (*kv) + 2] = -1;
  }
  
  for(int j = 0; j < *kv; ++j)
      AB[((*la)-1) * (*lab) + j] = 0;
    
  AB[((*la)-1) * (*lab) + (*kv)] = -1;
  AB[((*la)-1) * (*lab) + (*kv) + 1] = 2;
  AB[((*la)-1) * (*lab) + (*kv) + 2] = 0;
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  // TODO: Fill AB with the identity matrix
  // Only the main diagonal should have 1, all other entries are 0

  for (int i = 1; i < *la; ++i)
  {
    for(int j = 0; j < *kv; ++j)
      AB[i * (*lab) + j] = 0;
    
    AB[i * (*lab )+ (*kv)] = 0;
    AB[i * (*lab ) + (*kv) + 1] = 1;
    AB[i * (*lab ) + (*kv) + 2] = 0;
  }
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  // TODO: Compute RHS vector
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  // TODO: Compute the exact analytical solution at each grid point
  // This depends on the source term f(x) used in set_dense_RHS_DBC_1D
}  

void set_grid_points_1D(double* x, int* la){
  // TODO: Generate uniformly spaced grid points in [0,1]
}

double relative_forward_error(double* x, double* y, int* la){
  // TODO: Compute the relative error using BLAS functions (dnrm2, daxpy or manual loop)
  return 0.0;
}

int indexABCol(int i, int j, int *lab){
  // TODO: Return the correct index formula for column-major band storage
  return 0;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  // TODO: Implement specialized LU factorization for tridiagonal matrices
  return *info;
}
