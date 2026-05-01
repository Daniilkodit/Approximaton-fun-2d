#include "../header.h"

void calculate_residual (int nx, int ny, double hx, double hy,
                          double *x, double *y, double a, double c,
                          double &r1, double &r2, double &r3, double &r4,
                          int func_id, double (*f) (int , double, double),
                          double *solution, int p, int k)
{
  auto Pf = [=, &x, &y, &solution](double x0, double y0) -> double
    {
      return approximation_fun(nx, ny, hx, hy, x, y, x0, y0, a, c, solution);
    };

  r1 = 0.0;
  r2 = 0.0;

  int i1, i2;
  thread_rows (nx + 1, p, k, i1, i2);
  i2 = std::min (i2, nx);
  for (int i = i1; i < i2; i++)
    {
      for (int j = 0; j < ny; j++)
        {
          double mean_x = (x[i] + x[i + 1] + x[i + 1]) / 3.0;
          double mean_y = (y[j] + y[j]     + y[j + 1]) / 3.0;

          double pval = approximation_fun(nx, ny, hx, hy, x, y,
                                           mean_x, mean_y, a, c, solution);
          double diff = fabs(f(func_id, mean_x, mean_y) - pval);
          if (diff > r1)
            r1 = diff;
          r2 += diff;
          mean_x = (x[i] + x[i]     + x[i + 1]) / 3.0;
          mean_y = (y[j] + y[j + 1] + y[j + 1]) / 3.0;

          diff = fabs (f (func_id, mean_x, mean_y) - Pf (mean_x, mean_y));
          if (diff > r1)
            r1 = diff;
          r2 += diff;
        }
    }
  r2 *= hx * hy / 2.0;

  r3 = 0.0;
  r4 = 0.0;
  for (int i = i1; i <= i2; i++)
    {
      for (int j = 0; j <= ny; j++)
        {
          double diff = fabs (f (func_id, x[i], y[j]) - Pf (x[i], y[j]));
          if (diff > r3)
            r3 = diff;
          r4 += diff;
        }
    }
  r4 *= hx * hy;
}
