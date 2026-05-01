#include "../header.h"

// This is necessary so that the answer does not change from launch to launch.
double
scalar_product (int n, const double *x, const double *y,
                  int p, int k, double *sp)
{
  int i1, i2;
  double s = 0.;
  thread_rows (n, p, k, i1, i2);
  for (int i = i1; i < i2; i++)
    {
      s += x[i] * y[i];
    }

  sp[k] = s;
  reduce_sum (p);
  s = 0.;
  for (int i = 0; i < p; i++)
    {
      s += sp[i];
    }
  reduce_sum (p); // maybe not need?
  return s;
}

void
mult_sub_vector (int n, double *x, const double *y,
                  double t, int p, int k)
{
  int i1, i2;
  thread_rows (n, p, k, i1, i2);
  for (int i = i1; i < i2; i++)
    {
      x[i] -= t * y[i];
    }

  reduce_sum (p);
}

void
mult_msr_matrix_vector (const double *A, const int *I, int n,
                          const double *x, double *b, int p, int k)
{
  int i1, i2;
  thread_rows (n, p, k, i1, i2);
  for (int i = i1; i < i2; i++)
    {
      double s = A[i] * x[i];
      int l = I[i + 1] - I[i];
      int j_start = I[i];
      for (int j = 0; j < l; j++)
        {
          s += A[j_start + j] * x[I[j_start + j]];
        }

      b[i] = s;
    }

  reduce_sum(p); // это может замедлять
}

double
norm_vector (int n, const double *b, int p, int k, double *sp)
{
  int i1, i2;
  thread_rows (n, p, k, i1, i2);
  double s = 0;
  for (int i = i1; i < i2; i++)
    {
      s += b[i] * b[i];
    }

  sp[k] = s;
  reduce_sum (p);
  s = 0.;
  for (int i = 0; i < p; i++)
    {
      s += sp[i];
    }
  reduce_sum (p);
  return sqrt (s);
}

void
thread_rows (int n, int p, int k, int &i1, int &i2)
{
  i1 = n * k;
  i2 = n * (k + 1);
  i1 /= p;
  i2 /= p;
}
