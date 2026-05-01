#include "../header.h"

double approximation_fun (int nx, int ny, double hx, double hy, double *x, double *y,
                            double x0, double y0, double a, double c, double *solution)
{
  // search need square
  int i = (x0 - a) / hx;
  int j = (y0 - c) / hy;

  if (i < 0)
    i = 0;
  if (i >= nx)
    i = nx - 1;
  if (j < 0)
    j = 0;
  if (j >= ny)
    j = ny - 1;
  // search need triangle
  double linear_equation = (x0 - x[i]) * (y[j + 1] - y[j])
                           - (y0 - y[j]) * (x[i + 1] - x[i]);
  int l1, l2, l3;
  if (linear_equation > 0)
    {
      ij_to_l (nx, ny, i    , j    , l1);
      ij_to_l (nx, ny, i + 1, j    , l2);
      ij_to_l (nx, ny, i + 1, j + 1, l3);
    }
  else
    {
      ij_to_l (nx, ny, i    , j    , l1);
      ij_to_l (nx, ny, i    , j + 1, l2);
      ij_to_l (nx, ny, i + 1, j + 1, l3);
    }

  // create plane
  double x1 = x[i]    , y1 = y[j]    , z1 = solution[l1];
  double x2 = x[i + 1], y2 = y[j]    , z2 = solution[l2];
  double x3 = x[i + 1], y3 = y[j + 1], z3 = solution[l3];

  if (linear_equation < 0)
    {
      x2 = x[i];
      y2 = y[j + 1];
    }

  // calculate plan
  double a12 = x2 - x1, b12 = y2 - y1, c12 = z2 - z1;
  double a13 = x3 - x1, b13 = y3 - y1, c13 = z3 - z1;

  double D  = a12 * b13 - b12 * a13;
  double Dx = b12 * c13 - c12 * b13;
  double Dy = a12 * c13 - c12 * a13;

  return (z1 - ((x0 - x1) * Dx - (y0 - y1) * Dy) / D);
}