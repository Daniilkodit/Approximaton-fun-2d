#include "../header.h"

int
minimal_errors_msr_martix_full (int n, const double *A, const int *I,
                                  const double *b, double *x, double *r,
                                  double *u, double *v, double eps,
                                  int maxit, int max_steps, int p, int k, double *sp)
{
  int it = 0, step = 0;
  for (step = 0; step < max_steps; step++)
    {
      int ret = minimal_errors_msr_martix (
          n, A, I, b, x, r, u, v, eps, maxit, p, k, sp);
      if (ret >= 0)
        {
          it += ret;
          break;
        }
      it += maxit; // make restart
    }

  if (step >= max_steps)
    {
      return -1;
    }

  return it;
}

int
minimal_errors_msr_martix (int n, const double *A, const int *I,
                              const double *b, double *x, double *r,
                              double *u, double *v, double eps,
                              int maxit, int p, int k, double *sp)
{
  int it = 0;
  double b_norm = norm_vector (n, b, p, k, sp);
  double prec = eps * eps * b_norm; // accuracy of find
  mult_msr_matrix_vector (A, I, n, x, r, p, k); // r = Ax
  mult_sub_vector (n, r, b, 1, p, k); // r -= 1*b
  for (it = 0; it < maxit; it++)
  {
    apply_preconditioner_msr_matrix (n, A, I, v, u, r, p ,k); // v = M^-1 * r
    mult_msr_matrix_vector (A, I, n, v, u, p, k); // u = Av
    double c1 = scalar_product (n, v, r, p, k, sp);
    double c2 = scalar_product (n, u, v, p, k, sp);
    if (c1 < prec || c2 < prec)
      {
        break;
      }

    double tau = c1 / c2;
    mult_sub_vector (n, x, v, tau, p, k);
    mult_sub_vector (n, r, u, tau, p, k);
  }

  if (it >= maxit)
    {
      return -1;
    }

  return it;
}