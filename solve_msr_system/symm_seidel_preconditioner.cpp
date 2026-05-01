#include "../header.h"

void
apply_preconditioner_msr_matrix (int n, const double *A, const int *I,
                                 double *v, double *u, double *r, int p, int k)
{
  if (k == 0)
    memset(u, 0, n * sizeof(double));
  reduce_sum(p);

  // first solve (D + L)u = r
  for (int step = 0; step < n && k == 0; step++)
    {
      double s = 0;
      int j_start = I[step];
      int len = I[step + 1] - I[step]; // len max = 6

      for (int i = 0; i < len; i++)
        {
          int col = I[j_start + i];
          if (col > step)
            continue;

          s += u[col] * A[j_start + i];
        }
      u[step] = (r[step] - s) / A[step];
    }
  reduce_sum (p);

  // second solve D^-1u = u;
  int i1, i2;
  thread_rows (n, p, k, i1, i2);
  for (int i = i1; i < i2; i++)
    {
      u[i] *= A[i];
    }
  reduce_sum (p);

  // third solve (D + R)v = u;
  for (int step = n; step > 0 && k == 0; step--)
    {
      double s = 0;
      int idx = step - 1;
      int j_start = I[idx];
      int len = I[idx + 1] - I[idx]; // len max = 6

      for (int i = len - 1; i >= 0; i--)
        {
          int col = I[j_start + i];
          if (col <= idx)
            continue;

          s += v[col] * A[j_start + i];
        }
      v[idx] = (u[idx] - s) / A[idx];
    }
  reduce_sum (p);
}