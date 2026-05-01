#include "../header.h"

#define F(I, J) (f (func_id, x0 + (I) * hx, y0 + (J) * hy))

void fill_right_side (int nx, int ny, double a, double b, double c, double d,
                        double *B, int p, int k, int func_id,
                        double (*f) (int func_id,double x, double y))
{
  double hx = (b - a) / nx;
  double hy = (d - c) / ny;
  int l1, l2, i, j;
  get_my_rows (nx, ny, p, k, l1, l2);
  for (int l = l1; l < l2; l++)
    {
      l_to_ij (nx, ny, i, j, l);
      B[l] = F_ij (nx, ny, hx, hy, a, c, i, j, func_id, f);
    }

  reduce_sum (p);
}

double F_ij (int nx, int ny, double hx, double hy,
              double x0, double y0, int i, int j,
              int func_id, double (*f) (int func_id,double x, double y))
{
  double w = hx * hy / 192.;
  if (i > 0 && i < nx && j > 0 && j < ny)
    {
      return w * (  36 * F (i      , j      )
                  + 20 * F (i + 0.5, j      )
                  + 20 * F (i      , j - 0.5)
                  + 20 * F (i - 0.5, j - 0.5)
                  + 20 * F (i - 0.5, j      )
                  + 20 * F (i      , j + 0.5)
                  + 20 * F (i + 0.5, j + 0.5)
                  + 4  * F (i + 0.5, j - 0.5)
                  + 4  * F (i - 0.5, j - 1  )
                  + 4  * F (i - 1  , j - 0.5)
                  + 4  * F (i - 0.5, j + 0.5)
                  + 4  * F (i + 0.5, j + 1  )
                  + 4  * F (i + 1  , j + 0.5)
                  + 2  * F (i + 1  , j      )
                  + 2  * F (i      , j - 1  )
                  + 2  * F (i - 1  , j - 1  )
                  + 2  * F (i - 1  , j      )
                  + 2  * F (i      , j + 1  )
                  + 2  * F (i + 1  , j + 1  ));
    }

  if (j == 0 && i > 0 && i < nx)
    {
      return w * (  18 * F (i      , j      )
                  + 10 * F (i + 0.5, j      )
                  + 10 * F (i - 0.5, j      )
                  + 20 * F (i      , j + 0.5)
                  + 20 * F (i + 0.5, j + 0.5)
                  + 4  * F (i - 0.5, j + 0.5)
                  + 4  * F (i + 0.5, j + 1  )
                  + 4  * F (i + 1  , j + 0.5)
                  + 1  * F (i + 1  , j      )
                  + 1  * F (i - 1  , j      )
                  + 2  * F (i      , j + 1  )
                  + 2  * F (i + 1  , j + 1  ));
    }

  if (j == ny && i > 0 && i < nx)
    {
      return w * (  18 * F (i      , j      )
                  + 10 * F (i + 0.5, j      )
                  + 10 * F (i - 0.5, j      )
                  + 20 * F (i      , j - 0.5)
                  + 20 * F (i - 0.5, j - 0.5)
                  + 4  * F (i + 0.5, j - 0.5)
                  + 4  * F (i - 0.5, j - 1  )
                  + 4  * F (i - 1  , j - 0.5)
                  + 2  * F (i      , j - 1  )
                  + 2  * F (i - 1  , j - 1  )
                  + 1  * F (i + 1  , j      )
                  + 1  * F (i - 1  , j      ));
    }

  if (i == 0 && j > 0 && j < ny)
    {
      return w * (  18 * F (i      , j      )
                  + 20 * F (i + 0.5, j      )
                  + 20 * F (i + 0.5, j + 0.5)
                  + 10 * F (i      , j - 0.5)
                  + 10 * F (i      , j + 0.5)
                  + 4  * F (i + 0.5, j - 0.5)
                  + 4  * F (i + 0.5, j + 1  )
                  + 4  * F (i + 1  , j + 0.5)
                  + 2  * F (i + 1  , j      )
                  + 2  * F (i + 1  , j + 1  )
                  + 1  * F (i      , j - 1  )
                  + 1  * F (i      , j + 1  ));
    }

  if (i == nx && j > 0 && j < ny)
    {
      return w * (  18 * F (i      , j      )
                  + 20 * F (i - 0.5, j - 0.5)
                  + 20 * F (i - 0.5, j      )
                  + 10 * F (i      , j - 0.5)
                  + 10 * F (i      , j + 0.5)
                  + 4  * F (i - 0.5, j - 1  )
                  + 4  * F (i - 1  , j - 0.5)
                  + 4  * F (i - 0.5, j + 0.5)
                  + 2  * F (i - 1  , j - 1  )
                  + 2  * F (i - 1  , j      )
                  + 1  * F (i      , j - 1  )
                  + 1  * F (i      , j + 1  ));
    }

  if (i == 0 && j == 0)
    {
      return w * (  12 * F (i      , j      )
                  + 10 * F (i + 0.5, j      )
                  + 10 * F (i      , j + 0.5)
                  + 20 * F (i + 0.5, j + 0.5)
                  + 4  * F (i + 1  , j + 0.5)
                  + 4  * F (i + 0.5, j + 1  )
                  + 1  * F (i + 1  , j      )
                  + 1  * F (i      , j + 1  )
                  + 2  * F (i + 1  , j + 1  ));
    }

  if (i == nx && j == ny)
    {
      return w * (  12 * F (i      , j      )
                  + 10 * F (i      , j - 0.5)
                  + 20 * F (i - 0.5, j - 0.5)
                  + 10 * F (i - 0.5, j      )
                  + 4  * F (i - 0.5, j - 1  )
                  + 4  * F (i - 1  , j - 0.5)
                  + 1  * F (i      , j - 1  )
                  + 2  * F (i - 1  , j - 1  )
                  + 1  * F (i - 1  , j      ));
    }

  if (i == 0 && j == ny)
    {
      return w * (  6  * F (i       , j     )
                  + 10 * F (i + 0.5, j      )
                  + 10 * F (i      , j - 0.5)
                  + 4  * F (i + 0.5, j - 0.5)
                  + 1  * F (i + 1  , j      )
                  + 1  * F (i      , j - 1  ));
    }

  if (i == nx && j == 0)
    {
      return w * (   6 * F (i      , j      )
                  + 10 * F (i - 0.5, j      )
                  + 10 * F (i      , j + 0.5)
                  + 4  * F (i - 0.5, j + 0.5)
                  + 1  * F (i - 1  , j      )
                  + 1  * F (i      , j + 1  ));
    }

  return -99999999;
}
