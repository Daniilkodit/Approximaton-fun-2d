#include "../header.h"

int
check_row_sum (int nx, int ny, double hx, double hy,
               int *I, double *A, int p, int k)
{
  double eps = 1.e-15;
  int l1, l2, err = 0;
  get_my_rows (nx, ny, p, k, l1, l2);

  for (int l = l1; l < l2; l++)
    {
      double s = A[l];
      int r = I[l + 1] - I[l];
      for (int t = 0; t < r; t++)
        {
          s += A[I[l] + t];
        }
      int i, j;
      l_to_ij (nx, ny, i, j, l);
      int count = get_triangles (nx, ny, i, j);
      if (fabs (s - count * hx * hy / 6.0) > eps * hx * hy)
        {
          err = -1;
          break;
        }
    }

  reduce_sum (p, &err, 1);
  return err;
}

int check_symm (int nx, int ny, double hx, double hy,
                int *I, double *A, int p, int k)
{
  int l1, l2;
  int err = 0;
  double w = sqrt (hx * hy);
  double eps = 1.e-15;
  get_my_rows (nx, ny, p, k, l1, l2);

  for (int l = l1; l < l2; l++)
    {
      int r = I[l + 1] - I[l];
      int m = I[l];
      for (int t = 0; t < r; t++)
        {
          int j = I[m + t];
          double alj = A[m + t];
          int rj = I[j + 1] - I[j];
          int mj = I[j], tj;

          for (tj = 0; tj < rj; tj++)
            {
              if (I[mj + tj] == l)
                break;
            }

          if (tj == rj)
            {
              err = -1;
              break;
            }

          if (fabs (alj - A[mj + tj]) > eps * w)
            {
              err = -2;
              break;
            }
        }

      if (err < 0)
        break;
    }

  reduce_sum (p, &err, 1);
  return err;
}

int
get_triangles (int nx, int ny, int i, int j)
{
  if (i > 0 && i < nx && j > 0 && j < ny)
    return 6;

  if ((i == 0 || i == nx) && j > 0 && j < ny)
    return 3;

  if ((j == 0 || j == ny) && i > 0 && i < nx)
    return 3;

  if ((i == 0 && j == 0)  || (i == nx && j == ny))
    return 2;

  if ((i == 0 && j == ny) || (i == nx && j == 0))
    return 1;

  return -999999999;
}